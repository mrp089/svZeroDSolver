from sympy import symbols, Matrix, simplify, pi, Abs, Max, sign, Heaviside, Piecewise, Function, cse
from sympy.printing.c import C99CodePrinter
import re
import pdb
import yaml
import argparse


class warp_signed(Function):
    """Signed phase warp into one cardiac cycle, matching
    ChamberSphere::get_elastance_values. Its derivative with respect to its
    argument is 1 (it shifts by an integer number of periods), so it composes
    transparently with the activation derivatives. The C++ side defines a
    ``warp_signed`` lambda."""

    def fdiff(self, argindex=1):
        return 1

    def _eval_is_real(self):
        return True


class CppPrinter(C99CodePrinter):
    """C printer that emits Heaviside as a single-line ternary (and the
    arithmetic form of sign), and prints the warp_signed helper as a function
    call. Avoids the verbose/NaN-prone default output for activation-function
    derivatives."""

    def _print_Heaviside(self, expr):
        return "(%s > 0 ? 1.0 : 0.0)" % self._print(expr.args[0])

    def _print_warp_signed(self, expr):
        return "warp_signed(%s)" % self._print(expr.args[0])


_printer = CppPrinter()


def ccode(expr):
    return _printer.doprint(expr)

def load_model(filepath):
    with open(filepath, 'r') as file:
        data = yaml.safe_load(file)

    variables = Matrix(symbols(' '.join(data['variables']), real=True))
    derivatives = Matrix(symbols(' '.join(data['derivatives']), real=True))
    constant_syms = symbols(' '.join(data['constants']), real=True)
    constants = Matrix(list(constant_syms) if hasattr(constant_syms, '__iter__')
                       else [constant_syms])

    context = {str(s): s for s in list(variables) + list(derivatives) + list(constants)}
    context['pi'] = pi
    context['abs'] = Abs
    context['warp_signed'] = warp_signed

    # Optional time symbol (e.g. the activation time `t`). It is neither a state
    # variable nor a calibration parameter (it gets no Jacobian column), but it
    # is available to the residual/helper expressions.
    time_symbol = None
    if 'time_symbol' in data:
        time_symbol = symbols(data['time_symbol'], real=True)
        context[str(time_symbol)] = time_symbol

    if 'helper_functions' in data:
        exec(data['helper_functions'], context, context)
    residual_exprs = [eval(res, context, context) for res in data['residuals']]

    residuals = Matrix(residual_exprs)
    time_dependent = {context[s] for s in data.get('time_dependent', [])}
    if time_symbol is not None:
        time_dependent.add(time_symbol)

    assert len(variables) == len(derivatives), f"Number of variables must be equal to number of derivatives"

    return variables, derivatives, constants, residuals, time_dependent, time_symbol

def extract_linear(residuals, y, dy):
    zero_subs = [(v, 0) for v in list(y) + list(dy)]
    nr = residuals.shape[0]
    ny = y.shape[0]
    E = Matrix.zeros(nr, ny)
    F = Matrix.zeros(nr, ny)
    for i in range(nr):
        for j in range(ny):
            for mat, dd in zip([E, F], [dy, y]):
                coeff = residuals[i].coeff(dd[j])
                const_coeff = simplify(coeff.subs(zero_subs))
                if const_coeff:
                    mat[i, j] = const_coeff
    return E, F

def extract_nonlinear(residuals, E, F, y, dy):
    C = simplify(residuals - E * dy - F * y)
    dC_dy = simplify(C.jacobian(y))
    dC_dydot = simplify(C.jacobian(dy))
    return C, dC_dy, dC_dydot

def extract_gradient(residuals, constants):
    # Parameter Jacobian dr_i/dalpha_j used by the calibrator's update_gradient.
    # Simplify only smooth entries; entries from differentiating Abs/Max (the
    # activation function) are left unsimplified so they stay compact and avoid
    # the sign = x/|x| form that NaNs at zero.
    def smooth(e):
        return not e.has(Abs, Max, sign, Heaviside, Piecewise)
    jac = residuals.jacobian(constants)
    return jac.applyfunc(lambda e: simplify(e) if smooth(e) else e)

def depends_on(expr, symbols_set):
    return any(sym in expr.free_symbols for sym in symbols_set)

def partition_terms(E, F, C, dC_dy, dC_dydot, variables, time_dependent_symbols):
    parts = {"constant": [], "time": [], "solution": []}
    variable_symbols = variables.free_symbols
    for label, mat in [("E", E), ("F", F), ("C", C), ("dC_dy", dC_dy), ("dC_dydot", dC_dydot)]:
        for i in range(mat.rows):
            for j in range(mat.cols):
                if mat.cols == 1:
                    expr = mat[i]
                else:
                    expr = mat[i, j]
                if expr != 0:
                    if depends_on(expr, variable_symbols):
                        target = "solution"
                    elif depends_on(expr, time_dependent_symbols):
                        target = "time"
                    else:
                        target = "constant"
                    if mat.cols == 1:
                        parts[target].append((label, i, expr))
                    else:
                        parts[target].append((label, i, j, expr))
    return parts

def get_expressions(parts, type):
    expressions = [part[-1] for part in parts]
    deps = set().union(*(expr.free_symbols for expr in expressions))
    return deps.intersection(type)

def replace_symbolic_indices(expr, base):
    pattern = re.compile(rf'\b{base}(\d+)\b')
    return pattern.sub(rf'{base}[global_var_ids[\1]]', expr)

def format_cpp_expr(expr):
    try:
        cpp_expr = ccode(expr).replace('3.141592653589793', 'M_PI')
    except:
        pdb.set_trace()
    cpp_expr = replace_symbolic_indices(cpp_expr, 'y')
    cpp_expr = replace_symbolic_indices(cpp_expr, 'dy')
    return cpp_expr

def print_index(i, j=None):
    print(f".coeffRef(global_eqn_ids[{i}]", end='')
    if j is not None:
        print(f", global_var_ids[{j}]", end='')
    print(f")", end='')

def print_system(parts):
    for part in parts:
        if part[-1] != 0:
            print(f"  system.{part[0]}", end='')
            if len(part) == 3:
                print_index(part[1])
            elif len(part) == 4:
                print_index(part[1], part[2])
            print(f" = {format_cpp_expr(part[-1])};")

def print_constants(parts, constants, time_dependent):
    for out in get_expressions(parts, constants):
        if out in time_dependent:
            print(f"  // compute time dependent constant {out}")
        else:
            print(f"  const double {out} = parameters[global_param_ids[ParamId::{out}]];")

def print_variables(parts, y, dy):
    for out in get_expressions(parts, y.free_symbols | dy.free_symbols):
        for name, vec in zip(['y', 'dy'], [y, dy]):
            if out in vec:
                index = next(i for i, sym in enumerate(vec) if sym == out)
                print(f"  const double {out} = {name}[global_var_ids[{index}]];")

def free_symbols_of(exprs):
    nonzero = [e for e in exprs if e != 0]
    if not nonzero:
        return set()
    return set().union(*(e.free_symbols for e in nonzero))

def uses_warp(exprs):
    return any(e.has(warp_signed) for e in exprs if e != 0)

def print_alpha_constants(exprs, constants):
    """Declare calibration parameters read from the ``alpha`` vector."""
    deps = free_symbols_of(exprs)
    for c in constants:
        if c in deps:
            print(f"  const double {c} = alpha[global_param_ids[ParamId::{c}]];")

def print_gradient_variables(exprs, y, dy):
    """Declare state variables (``y``) and derivatives (``dy``)."""
    deps = free_symbols_of(exprs)
    for name, vec in zip(['y', 'dy'], [y, dy]):
        for index, sym in enumerate(vec):
            if sym in deps:
                print(f"  const double {sym} = {name}[global_var_ids[{index}]];")

def print_time(exprs, time_symbol):
    """Declare the in-cycle time and (if used) the periodic phase-warp helper."""
    if time_symbol is None or time_symbol not in free_symbols_of(exprs):
        return
    print(f"  const double {time_symbol} = "
          "model->cardiac_cycle_period > 0.0 ? "
          "fmod(model->time, model->cardiac_cycle_period) : model->time;")
    if uses_warp(exprs):
        print("  const double T_cardiac = model->cardiac_cycle_period;")
        print("  auto warp_signed = [T_cardiac](double dt) {")
        print("    return T_cardiac > 0.0 ? "
              "fmod(dt + 1.5 * T_cardiac, T_cardiac) - 0.5 * T_cardiac : dt;")
        print("  };")

def print_gradient(residuals, jacobian, constants, y, dy, time_symbol):
    names = [str(c) for c in constants]

    # Collect the non-zero residual and Jacobian entries. Columns are grouped by
    # parameter for readability.
    items = []  # (line_prefix, expr)
    for i in range(residuals.shape[0]):
        if residuals[i] != 0:
            items.append((f"  residual(global_eqn_ids[{i}])", residuals[i]))
    for j in range(jacobian.cols):
        for i in range(jacobian.rows):
            if jacobian[i, j] != 0:
                items.append((f"  jacobian.coeffRef(global_eqn_ids[{i}], "
                              f"global_param_ids[ParamId::{names[j]}])",
                              jacobian[i, j]))

    all_exprs = [expr for _, expr in items]

    # Common-subexpression elimination: the activation-function derivatives share
    # large repeated subexpressions (the tanh terms, S_plus/S_minus, ...). cse
    # hoists them into temporaries so the columns stay compact.
    replacements, reduced = cse(all_exprs)

    print('update_gradient')
    print_alpha_constants(all_exprs, constants)
    print_gradient_variables(all_exprs, y, dy)
    print_time(all_exprs, time_symbol)
    print()
    for sym, sub in replacements:
        print(f"  const double {sym} = {format_cpp_expr(sub)};")
    if replacements:
        print()
    for (prefix, _), expr in zip(items, reduced):
        print(f"{prefix} = {format_cpp_expr(expr)};")
    print()

def main(yaml_path):
    # read model from yaml file
    y, dy, constants, residuals, time_dependent, time_symbol = load_model(yaml_path)

    # extract linear and nonlinear terms
    E, F = extract_linear(residuals, y, dy)
    C, dC_dy, dC_dydot = extract_nonlinear(residuals, E, F, y, dy)

    # split into constant, time dependent and solution parts
    parts = partition_terms(E, F, C, dC_dy, dC_dydot, Matrix(list(y) + list(dy)), time_dependent)

    # print C++ code
    for section in ['constant', 'time', 'solution']:
        print(f'update_{section}')
        print_constants(parts[section], constants, time_dependent)
        print_variables(parts[section], y, dy)
        print_system(parts[section])
        print()

    # parameter Jacobian for the calibrator (update_gradient)
    jac = extract_gradient(residuals, constants)
    print_gradient(residuals, jac, constants, y, dy, time_symbol)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process YAML model file to generate C++ code for svZeroDSolver')
    parser.add_argument('yaml_file', help='Path to YAML model file')
    args = parser.parse_args()
    main(args.yaml_file)