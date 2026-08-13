// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberCylinder.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <utility>

#include "Model.h"

namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

// The continuum + active-law kernel is templated on the scalar type so it can
// be evaluated in double precision (residual) or with a complex perturbation
// (complex-step differentiation) to obtain machine-accuracy tangents without
// finite-difference cancellation.
template <typename T>
using Mat3 = Eigen::Matrix<T, 3, 3>;
template <typename T>
using Vec3 = Eigen::Matrix<T, 3, 1>;
using cplx = std::complex<double>;

/// Real part of a scalar (identity for double, Re(.) for complex) so ordering
/// comparisons can be made on a complex-step argument.
inline double re(double x) { return x; }
inline double re(const cplx& x) { return x.real(); }

/// Upper clamp that stays holomorphic for complex-step: below the cap the value
/// (and its imaginary perturbation) pass through; at/above the cap it is the
/// real constant `m` (zero derivative, matching the clamped plateau).
template <typename T>
inline T clamp_max(T x, double m) {
  return re(x) < m ? x : T(m);
}

/// Material parameters passed to the continuum kernel.
struct MatParams {
  double C1, C2, C3, C4, C5, C6, kappa, gamma;
  bool mixed = false;    ///< true: bulk term uses the mixed pressure p_mix
  double p_mix = 0.0;    ///< mixed u/p element hydrostatic pressure (if mixed)
};

/// Active-stress input at an integration point. Simple model: active fiber
/// stress magnitude `tau`. Genet BCS: bond stress `tau_c` and contractile-
/// element strain rate `ecdot`, giving sigma_1D = (tau_c + mu*ecdot)/sqrt(I4).
struct ActiveInput {
  bool bcs = false;
  double tau = 0.0;    ///< simple model active fiber stress magnitude
  double tauc = 0.0;   ///< BCS active bond stress tau_c
  double ecdot = 0.0;  ///< BCS contractile-element strain rate
  double mu = 0.0;     ///< BCS active dissipation
  double i4pow = 0.5;  ///< sigma_1D = T_fib / I4^i4pow; 0.5 = Genet Eq. 30
                       ///< (T_fib/(1+e_fib)), 1.0 = Eq. 59 limit (T_fib/I4)
};

/// Full contraction A:B of two symmetric 3x3 tensors (no conjugation, so it is
/// the holomorphic bilinear form required by complex-step).
template <typename T>
inline T contract(const Mat3<T>& A, const Mat3<T>& B) {
  return (A.array() * B.array()).sum();
}

/**
 * @brief Kinematics at a through-wall integration point.
 *
 * Holds the right Cauchy-Green tensor C and the six derivatives of the
 * Green-Lagrange strain with respect to the local kinematic quantities
 * xi = [rho, rho', phi', beta, eps, eta'] (\cite genet23 Eqs. 15, 16, 21).
 */
template <typename T>
struct Kinematics {
  Mat3<T> C;      ///< Right Cauchy-Green deformation tensor
  Mat3<T> dE[6];  ///< dE/dxi_k for xi=[rho,rho',phi',beta,eps,eta']
  T J;            ///< Volume ratio det(F)
};

/// Build C, J and the strain derivatives from local kinematics at radius R.
template <typename T>
Kinematics<T> compute_kinematics(const T xi[6], double R) {
  const T rho = xi[0], drho = xi[1], dphi = xi[2];
  const T beta = xi[3], eps = xi[4], deta = xi[5];
  const T a = 1.0 + drho;      // radial stretch F_RR
  const T b = 1.0 + rho / R;   // circumferential stretch F_ThTh
  const T e = 1.0 + eps;       // longitudinal stretch F_ZZ

  Kinematics<T> k;
  // Right Cauchy-Green tensor in the orthonormal (e_R, e_Th, e_Z) basis.
  const T C00 = a * a + (R * b * dphi) * (R * b * dphi) + deta * deta;
  const T C11 = b * b;
  const T C22 = (R * b * beta) * (R * b * beta) + e * e;
  const T C01 = R * b * b * dphi;
  const T C02 = R * R * b * b * dphi * beta + deta * e;
  const T C12 = R * b * b * beta;
  k.C << C00, C01, C02, C01, C11, C12, C02, C12, C22;
  k.J = a * b * e;

  for (int i = 0; i < 6; i++) k.dE[i].setZero();

  // dE/drho (with d b / d rho = 1/R)
  {
    Mat3<T>& m = k.dE[0];
    m(0, 0) = R * b * dphi * dphi;
    m(1, 1) = b / R;
    m(2, 2) = R * b * beta * beta;
    m(0, 1) = m(1, 0) = b * dphi;
    m(0, 2) = m(2, 0) = R * b * dphi * beta;
    m(1, 2) = m(2, 1) = b * beta;
  }
  // dE/drho'
  k.dE[1](0, 0) = a;
  // dE/dphi'
  {
    Mat3<T>& m = k.dE[2];
    m(0, 0) = R * R * b * b * dphi;
    m(0, 1) = m(1, 0) = 0.5 * R * b * b;
    m(0, 2) = m(2, 0) = 0.5 * R * R * b * b * beta;
  }
  // dE/dbeta
  {
    Mat3<T>& m = k.dE[3];
    m(2, 2) = R * R * b * b * beta;
    m(0, 2) = m(2, 0) = 0.5 * R * R * b * b * dphi;
    m(1, 2) = m(2, 1) = 0.5 * R * b * b;
  }
  // dE/deps
  {
    Mat3<T>& m = k.dE[4];
    m(2, 2) = e;
    m(0, 2) = m(2, 0) = 0.5 * deta;
  }
  // dE/deta'
  {
    Mat3<T>& m = k.dE[5];
    m(0, 0) = deta;
    m(0, 2) = m(2, 0) = 0.5 * e;
  }
  return k;
}

/// Myofiber structural tensor e_F (x) e_F at helix angle alpha (\cite genet23
/// Eq. 22). e_F lies in the circumferential-longitudinal plane. The angle is a
/// fixed geometric quantity, so it is evaluated in double and cast to T.
template <typename T>
inline Mat3<T> fiber_tensor(double alpha) {
  const double ca = std::cos(alpha), sa = std::sin(alpha);
  Vec3<T> ef(T(0.0), T(ca), T(sa));
  return ef * ef.transpose();
}

/**
 * @brief Second Piola-Kirchhoff stress at an integration point.
 *
 * Sigma = Sigma^d (transversely isotropic, \cite genet23 Eq. 25) + Sigma^b
 * (incompressibility penalty) + Sigma^v (viscous, Eq. 28) + Sigma^a (active
 * fiber stress, Eq. 29).
 */
template <typename T>
Mat3<T> compute_stress(const Kinematics<T>& k, const T xidot[6],
                       const ActiveInput& act, double alpha,
                       const MatParams& p) {
  const Mat3<T>& C = k.C;
  const T J = k.J;
  const Mat3<T> Cinv = C.inverse();
  const Mat3<T> I = Mat3<T>::Identity();
  const Mat3<T> M = fiber_tensor<T>(alpha);

  const T I1 = C.trace();
  const T trC2 = (C.array() * C.array()).sum();
  const T I2 = 0.5 * (I1 * I1 - trC2);
  const T I4 = (C * M).trace();  // e_F . C . e_F

  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const T Jm43 = Jm23 * Jm23;
  const T I1b = Jm23 * I1;
  const T I2b = Jm43 * I2;
  const T I4b = Jm23 * I4;
  (void)I2b;  // reduced invariant I2b enters only through its derivative below

  // Derivatives of the strain energy w.r.t. the isochoric invariants. The
  // exponent is clamped to keep values finite for non-physical trial states
  // (e.g. during sparsity-pattern setup); physical arguments are small.
  const T a1 = p.C4 * (I1b - 3.0) * (I1b - 3.0);
  const T a4 = p.C6 * (I4b - 1.0) * (I4b - 1.0);
  const T w1 =
      p.C1 + 2.0 * p.C3 * p.C4 * (I1b - 3.0) * std::exp(clamp_max(a1, 300.0));
  const T w2 = T(p.C2);
  const T w4 =
      2.0 * p.C5 * p.C6 * (I4b - 1.0) * std::exp(clamp_max(a4, 300.0));

  // d(isochoric invariant)/dC
  const Mat3<T> dI1b = Jm23 * (I - (I1 / 3.0) * Cinv);
  const Mat3<T> dI2b = Jm43 * (I1 * I - C - (2.0 / 3.0) * I2 * Cinv);
  const Mat3<T> dI4b = Jm23 * (M - (I4 / 3.0) * Cinv);

  Mat3<T> Sd = 2.0 * (w1 * dI1b + w2 * dI2b + w4 * dI4b);   // deviatoric
  // Bulk (incompressibility) 2nd-PK stress Sigma^b = Pi J C^{-1}. Penalty:
  // Pi = kappa(J-1); mixed u/p: Pi = p_mix, the element hydrostatic pressure
  // (held fixed under the xi complex-step; its column is added analytically).
  const T bulk = p.mixed ? T(p.p_mix) : p.kappa * (J - 1.0);
  Mat3<T> Sb = bulk * J * Cinv;                            // bulk

  Mat3<T> Edot = Mat3<T>::Zero();                           // viscous
  for (int m = 0; m < 6; m++) Edot += k.dE[m] * xidot[m];
  Mat3<T> Sv = p.gamma * Edot;

  // active fiber stress magnitude
  T sigma_act;
  if (act.bcs) {
    // sigma_1D = (tau_c + mu * e_c_dot) / (1 + e_fib), 1 + e_fib = sqrt(I4).
    // Written via the fiber tension tau_c + mu*e_c_dot (= T_fib) so the stiff
    // series stiffness k_s stays out of the mechanical residual.
    sigma_act = (act.tauc + act.mu * act.ecdot) / std::pow(I4, act.i4pow);
  } else {
    sigma_act = T(act.tau);
  }
  Mat3<T> Sa = sigma_act * M;

  return Sd + Sb + Sv + Sa;
}

/// Generalized internal forces g_k = Sigma : dE/dxi_k at an integration point.
template <typename T>
void compute_g(const T xi[6], const T xidot[6], const ActiveInput& act,
               double R, double alpha, const MatParams& p, T g[6]) {
  Kinematics<T> k = compute_kinematics<T>(xi, R);
  Mat3<T> S = compute_stress<T>(k, xidot, act, alpha, p);
  for (int i = 0; i < 6; i++) g[i] = contract<T>(S, k.dE[i]);
}

/// BCS Frank-Starling force-length function n_0(e_c): the calibrated
/// piecewise-linear length-dependence curve from the MEDISIM PhysioBlocks
/// reference implementation (n_0 = interp(e_c, abscissas, ordinates); Caruel et
/// al. 2013, \cite genet23 Ref. 13). It represents the actin-myosin overlap
/// (length-tension) relationship: n_0 rises from 0, plateaus at 1 over the
/// physiological contractile strain e_c in [0.20, 0.47], and falls back to 0.
/// NOTE: the plateau sits at e_c ~ 0.2, not ~1 as Caruel Fig. 7(a) suggests -
/// that figure is drawn over the wide isotonic papillary-muscle strain range.
/// (The center/width of the earlier Gaussian placeholder are unused.)
template <typename T>
inline T frank_starling(T e_c, double /*center*/, double /*width*/) {
  // Breakpoints (e_c, n_0) from PhysioBlocks (physioblocks/physioblocks).
  constexpr int NP = 9;
  static const double X[NP] = {-0.1668, -0.0073, 0.0534, 0.0969, 0.1326,
                               0.2016,  0.4663,  0.9187, 1.1762};
  static const double Y[NP] = {0.0,    0.5614, 0.7748, 0.8933, 0.9618,
                               1.0,    1.0,    0.1075, 0.0};
  const double x = re(e_c);
  if (x <= X[0]) return T(Y[0]);
  if (x >= X[NP - 1]) return T(Y[NP - 1]);
  int k = 0;
  while (k < NP - 1 && x >= X[k + 1]) ++k;
  const double slope = (Y[k + 1] - Y[k]) / (X[k + 1] - X[k]);
  // e_c may carry a complex-step perturbation; keep it in the linear term so
  // Im(n_0)/h recovers the exact segment slope.
  return Y[k] + slope * (e_c - X[k]);
}

/// BCS active-law parameters.
struct BcsParams {
  double k_s, k_0, mu, alpha, sigma0, n0_center, n0_width;
};

/// Residuals of the three BCS internal-variable ODEs (Genet Eqs. 32-33 /
/// Chapelle Eq. 9) at a quadrature point. a = [e_c, tau_c, k_c], ad = its rate;
/// the contractile-element strain rate is e_c_dot = ad[0].
template <typename T>
void bcs_residual(const T a[3], const T ad[3], T e_fib,
                  double nu_abs, double nu_plus, const BcsParams& bp,
                  T res[3]) {
  const T ecdot = ad[0];
  // Smoothed |e_c dot| (eps ~ 1e-2 /s regularizes the force-velocity kink so
  // Newton and its complex-step tangent behave; negligible vs |nu| physically).
  const T abs_ecd = std::sqrt(ecdot * ecdot + 1e-4);
  const T n0 = frank_starling<T>(a[0], bp.n0_center, bp.n0_width);
  const T decay = nu_abs + bp.alpha * abs_ecd;
  // series-spring force balance: mu * e_c_dot = k_s (e_fib - e_c) - tau_c
  res[0] = ecdot - (bp.k_s * (e_fib - a[0]) - a[1]) / bp.mu;
  // tau_c dot = -decay*tau_c + e_c_dot*k_c + n0*sigma0*nu_+
  res[1] = ad[1] - (-decay * a[1] + ecdot * a[2] + n0 * bp.sigma0 * nu_plus);
  // k_c dot = -decay*k_c + n0*k_0*nu_+
  res[2] = ad[2] - (-decay * a[2] + n0 * bp.k_0 * nu_plus);
}

}  // namespace

void ChamberCylinder::setup_dofs(DOFHandler& dofhandler) {
  // Read geometry and discretization (parameter values are available here).
  const double ri = model->get_parameter_value(global_param_ids[ParamId::Ri]);
  const double re = model->get_parameter_value(global_param_ids[ParamId::Re]);
  const double a_endo =
      model->get_parameter_value(global_param_ids[ParamId::alpha_endo]);
  const double a_epi =
      model->get_parameter_value(global_param_ids[ParamId::alpha_epi]);
  n_ele = std::max(
      1, (int)std::lround(
             model->get_parameter_value(global_param_ids[ParamId::num_elements])));
  n_node = n_ele + 1;
  use_bcs =
      model->get_parameter_value(global_param_ids[ParamId::active_model]) > 0.5;
  is_dynamic =
      model->get_parameter_value(global_param_ids[ParamId::use_inertia]) > 0.5;
  use_mixed =
      model->get_parameter_value(global_param_ids[ParamId::mixed]) > 0.5;

  // Uniform nodes from endocardium (index 0, R_i) to epicardium (R_e).
  node_R.resize(n_node);
  const double h = (re - ri) / n_ele;
  for (int i = 0; i < n_node; i++) node_R[i] = ri + i * h;

  // Three-point Gauss quadrature per element. The reference cylindrical volume
  // element R dR dTheta dZ is used, so the radius R is folded into the weight.
  const double gp[3] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
  const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  const double deg2rad = M_PI / 180.0;
  quad.clear();
  for (int e = 0; e < n_ele; e++) {
    const double Ra = node_R[e], Rb = node_R[e + 1];
    for (int q = 0; q < 3; q++) {
      QuadPoint pt;
      pt.elem = e;
      pt.R = 0.5 * (Ra + Rb) + 0.5 * (Rb - Ra) * gp[q];
      pt.w = gw[q] * 0.5 * (Rb - Ra) * pt.R;  // Gauss weight * dR/dxi * R
      const double frac = (pt.R - ri) / (re - ri);
      pt.alpha = deg2rad * (a_endo * (1.0 - frac) + a_epi * frac);
      pt.N[0] = (Rb - pt.R) / (Rb - Ra);
      pt.N[1] = (pt.R - Ra) / (Rb - Ra);
      pt.dN[0] = -1.0 / (Rb - Ra);
      pt.dN[1] = 1.0 / (Rb - Ra);
      quad.push_back(pt);
    }
  }
  n_quad = (int)quad.size();

  // Internal variables: rho/phi/eta fields, beta, eps, active block, volume.
  std::list<std::string> int_names;
  for (int a = 0; a < n_node; a++) int_names.push_back("rho_" + std::to_string(a));
  for (int a = 0; a < n_node; a++) int_names.push_back("phi_" + std::to_string(a));
  for (int a = 0; a < n_node; a++) int_names.push_back("eta_" + std::to_string(a));
  int_names.push_back("beta");
  int_names.push_back("eps");
  if (use_bcs) {
    for (int q = 0; q < n_quad; q++) {
      int_names.push_back("ec_" + std::to_string(q));
      int_names.push_back("tauc_" + std::to_string(q));
      int_names.push_back("kc_" + std::to_string(q));
    }
  } else {
    int_names.push_back("tau");
  }
  int_names.push_back("volume");

  // Full dynamics: velocity companion fields w = d(field)/dt (genet23 Eqs. 8,
  // 18, 45, A5-A6). Appended after the volume DOF; absent when quasi-static.
  if (is_dynamic) {
    for (int a = 0; a < n_node; a++) int_names.push_back("vrho_" + std::to_string(a));
    for (int a = 0; a < n_node; a++) int_names.push_back("vphi_" + std::to_string(a));
    for (int a = 0; a < n_node; a++) int_names.push_back("veta_" + std::to_string(a));
    int_names.push_back("vbeta");
    int_names.push_back("veps");
  }

  // Mixed u/p: one element-wise-constant (P0) hydrostatic pressure per element,
  // appended last so the penalty/dynamics layouts above are unchanged.
  if (use_mixed) {
    for (int e = 0; e < n_ele; e++)
      int_names.push_back("pmix_" + std::to_string(e));
  }

  // 3*n_node field eqns + beta + eps + active eqns + volume + mass + pressure
  const int n_vel = is_dynamic ? 3 * n_node + 2 : 0;
  const int n_pmix = use_mixed ? n_ele : 0;
  n_eqn = 3 * n_node + 2 + n_active_var() + 3 + n_vel + n_pmix;
  Block::setup_dofs_(dofhandler, n_eqn, int_names);
  n_var = 4 + 3 * n_node + 2 + n_active_var() + 1 + n_vel + n_pmix;

  // Dense block for dC_dy and dC_dydot; a few linear entries for F/E.
  num_triplets.D = n_eqn * n_var;
  num_triplets.F = 3 * n_node + 16 + (is_dynamic ? 3 * n_node + 4 : 0);
  num_triplets.E = (use_bcs ? 1 : 4) + 1 +  // +1 for the C_valve dPv/dt term
                   (is_dynamic ? 21 * n_node + 8 : 0);  // mass + companion
}

void ChamberCylinder::update_constant(SparseSystem& system,
                                      std::vector<double>& parameters) {
  // Simple active stress ODE: d(tau)/dt term. (The BCS model puts all of its
  // rate terms in the nonlinear residual C, so it needs no E/F entries here.)
  if (!use_bcs) {
    system.E.coeffRef(global_eqn_ids[e_active()], global_var_ids[i_tau()]) = 1.0;
  }

  // Mass conservation: Qin - Qout - Vdot - C_valve * dPv/dt = 0
  // (Vdot handled in update_solution; the compliance term is linear in dPin/dt).
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[1]) = 1.0;   // Qin
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[3]) = -1.0;  // Qout
  const double cvalve = parameters[global_param_ids[ParamId::c_valve]];
  system.E.coeffRef(global_eqn_ids[e_mass()], global_var_ids[0]) = -cvalve;

  // Pressure equality: Pin - Pout = 0.
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[2]) = -1.0;

  // Cavity volume definition: V - pi (R_i + rho_i)^2 (1+eps) L = 0 (linear V).
  system.F.coeffRef(global_eqn_ids[e_volume()], global_var_ids[i_vol()]) = 1.0;

  // Rigid-body pins: phi(R_i) = 0 and eta(R_i) = 0 (inner node, index 0).
  system.F.coeffRef(global_eqn_ids[e_phi(0)], global_var_ids[i_phi(0)]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_eta(0)], global_var_ids[i_eta(0)]) = 1.0;

  // --- Full dynamics: consistent mass matrix + velocity companion equations ---
  // (genet23 Eqs. 8,18,45; App. A.2 Eqs. A5-A6). The velocity field is the
  // push-forward v = Du[zeta] zeta_dot; introducing companion DOFs w = zeta_dot
  // turns the 2nd-order system into 1st order. The inertia force on field k is
  // sum_j M_kj w_dot_j with the consistent mass M = int rho0 (Du)^T Du dOmega.
  // Du (Eq. A5) is orthogonal between the radial/azimuthal (rho,phi,beta) and
  // axial (eta,eps) groups, so after integrating Theta (2 pi) and Z the only
  // couplings are rho-rho, eta-eta, phi-phi, eps-eps, beta-beta, phi-beta and
  // eta-eps. The mass is evaluated at the reference configuration (R), a
  // standard approximation - inertia is ~1e-4 of the internal/pressure forces
  // for cardiac parameters. The O(zeta_dot^2) centrifugal term D2u(zd,zd) is
  // likewise negligible and omitted.
  if (is_dynamic) {
    const double rho0 = parameters[global_param_ids[ParamId::density]];
    const double Lp = parameters[global_param_ids[ParamId::length]];
    const double cL1 = rho0 * 2.0 * M_PI * Lp;                    // int_Z 1  = L
    const double cL2 = rho0 * 2.0 * M_PI * Lp * Lp / 2.0;         // int_Z Z  = L^2/2
    const double cL3 = rho0 * 2.0 * M_PI * Lp * Lp * Lp / 3.0;    // int_Z Z^2= L^3/3
    // Accumulate the mass over quadrature points into a local map, then ASSIGN
    // (=) to the sparse system: update_constant may be called multiple times
    // (e.g. per step-size change), so += into the shared matrix would double it.
    std::map<std::pair<int, int>, double> M;
    auto Madd = [&](int r, int c, double v) { M[{r, c}] += v; };
    for (int q = 0; q < n_quad; q++) {
      const QuadPoint& pt = quad[q];
      const int A[2] = {pt.elem, pt.elem + 1};
      const double N[2] = {pt.N[0], pt.N[1]};
      const double R2 = pt.R * pt.R;
      Madd(e_beta(), i_vbeta(), cL3 * pt.w * R2);   // beta-beta
      Madd(e_eps(), i_veps(), cL3 * pt.w);          // eps-eps
      for (int i = 0; i < 2; i++) {
        Madd(e_beta(), i_vphi(A[i]), cL2 * pt.w * R2 * N[i]);  // beta-phi
        Madd(e_eps(), i_veta(A[i]), cL2 * pt.w * N[i]);        // eps-eta
        if (A[i] != 0) {  // skip the pinned phi_0/eta_0 momentum rows
          Madd(e_phi(A[i]), i_vbeta(), cL2 * pt.w * R2 * N[i]);
          Madd(e_eta(A[i]), i_veps(), cL2 * pt.w * N[i]);
        }
        for (int j = 0; j < 2; j++) {
          const double mrr = cL1 * pt.w * N[i] * N[j];       // rho-rho / eta-eta
          Madd(e_rho(A[i]), i_vrho(A[j]), mrr);
          if (A[i] != 0) {
            Madd(e_eta(A[i]), i_veta(A[j]), mrr);
            Madd(e_phi(A[i]), i_vphi(A[j]), cL1 * pt.w * R2 * N[i] * N[j]);
          }
        }
      }
    }
    for (const auto& kv : M)
      system.E.coeffRef(global_eqn_ids[kv.first.first],
                        global_var_ids[kv.first.second]) = kv.second;
    // Velocity companion equations: w_k - d(field_k)/dt = 0.
    auto companion = [&](int erow, int ivel, int ifield) {
      system.F.coeffRef(global_eqn_ids[erow], global_var_ids[ivel]) = 1.0;
      system.E.coeffRef(global_eqn_ids[erow], global_var_ids[ifield]) = -1.0;
    };
    for (int a = 0; a < n_node; a++) {
      companion(e_vrho(a), i_vrho(a), i_rho(a));
      companion(e_vphi(a), i_vphi(a), i_phi(a));
      companion(e_veta(a), i_veta(a), i_eta(a));
    }
    companion(e_vbeta(), i_vbeta(), i_beta());
    companion(e_veps(), i_veps(), i_eps());
  }
}

void ChamberCylinder::update_time(SparseSystem& system,
                                  std::vector<double>& parameters) {
  get_activation(parameters);  // sets act = |nu|, act_plus = |nu|_+
  // Simple active stress ODE: a(t) * tau term. (BCS uses nu inside C.)
  if (!use_bcs) {
    system.F.coeffRef(global_eqn_ids[e_active()], global_var_ids[i_tau()]) = act;
  }
}

void ChamberCylinder::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  MatParams p;
  p.C1 = parameters[global_param_ids[ParamId::C1]];
  p.C2 = parameters[global_param_ids[ParamId::C2]];
  p.C3 = parameters[global_param_ids[ParamId::C3]];
  p.C4 = parameters[global_param_ids[ParamId::C4]];
  p.C5 = parameters[global_param_ids[ParamId::C5]];
  p.C6 = parameters[global_param_ids[ParamId::C6]];
  p.kappa = parameters[global_param_ids[ParamId::kappa]];
  p.gamma = parameters[global_param_ids[ParamId::gamma]];
  p.mixed = use_mixed;  // bulk term uses the per-element pressure DOF below
  const double sigma_max = parameters[global_param_ids[ParamId::sigma_max]];
  const double Rip = parameters[global_param_ids[ParamId::Ri]];
  const double Lp = parameters[global_param_ids[ParamId::length]];

  BcsParams bp;
  bp.k_s = parameters[global_param_ids[ParamId::k_s]];
  bp.k_0 = parameters[global_param_ids[ParamId::k_0]];
  bp.mu = parameters[global_param_ids[ParamId::mu]];
  bp.alpha = parameters[global_param_ids[ParamId::bcs_alpha]];
  bp.sigma0 = sigma_max;
  bp.n0_center = parameters[global_param_ids[ParamId::n0_center]];
  bp.n0_width = parameters[global_param_ids[ParamId::n0_width]];
  const double i4pow = parameters[global_param_ids[ParamId::active_i4pow]];

  // Gather block DOFs and their rates.
  auto Y = [&](int lv) { return y[global_var_ids[lv]]; };
  auto Yd = [&](int lv) { return dy[global_var_ids[lv]]; };
  const double beta = Y(i_beta()), eps = Y(i_eps()), tau = Y(i_tau());
  const double beta_d = Yd(i_beta()), eps_d = Yd(i_eps());
  const double Pin = Y(0);

  // Local dense residual and tangents (assigned to the sparse system at end).
  std::vector<double> Cloc(n_eqn, 0.0);
  std::vector<double> Kloc(n_eqn * n_var, 0.0);
  std::vector<double> Kdloc(n_eqn * n_var, 0.0);
  auto Kat = [&](int r, int c) -> double& { return Kloc[r * n_var + c]; };
  auto Kdat = [&](int r, int c) -> double& { return Kdloc[r * n_var + c]; };

  // --- Finite element assembly of the internal virtual work over thickness ---
  for (int q = 0; q < n_quad; q++) {
    const QuadPoint& pt = quad[q];
    const int A0 = pt.elem, A1 = pt.elem + 1;
    const double N0 = pt.N[0], N1 = pt.N[1], dN0 = pt.dN[0], dN1 = pt.dN[1];

    // Interpolate local kinematics xi = [rho, rho', phi', beta, eps, eta'].
    const double rho = N0 * Y(i_rho(A0)) + N1 * Y(i_rho(A1));
    const double drho = dN0 * Y(i_rho(A0)) + dN1 * Y(i_rho(A1));
    const double dphi = dN0 * Y(i_phi(A0)) + dN1 * Y(i_phi(A1));
    const double deta = dN0 * Y(i_eta(A0)) + dN1 * Y(i_eta(A1));
    const double xi[6] = {rho, drho, dphi, beta, eps, deta};

    // Interpolate rates (for the viscous term).
    const double rho_d = N0 * Yd(i_rho(A0)) + N1 * Yd(i_rho(A1));
    const double drho_d = dN0 * Yd(i_rho(A0)) + dN1 * Yd(i_rho(A1));
    const double dphi_d = dN0 * Yd(i_phi(A0)) + dN1 * Yd(i_phi(A1));
    const double deta_d = dN0 * Yd(i_eta(A0)) + dN1 * Yd(i_eta(A1));
    const double xidot[6] = {rho_d, drho_d, dphi_d, beta_d, eps_d, deta_d};

    // Active-stress input: simple model uses the global tau; BCS uses this
    // quadrature point's contractile-element strain e_c and computes the fiber
    // stress from the fiber stretch internally.
    ActiveInput ai;
    ai.bcs = use_bcs;
    ai.i4pow = i4pow;
    if (use_bcs) {
      ai.tauc = Y(i_tauc(q));
      ai.ecdot = Yd(i_ec(q));
      ai.mu = bp.mu;
    } else {
      ai.tau = tau;
    }

    // Mixed u/p: this element's hydrostatic pressure enters the bulk stress.
    if (use_mixed) p.p_mix = Y(i_pmix(pt.elem));

    // Base (double) forces plus the analytic active/viscous rate derivatives.
    Kinematics<double> kin = compute_kinematics<double>(xi, pt.R);
    Matrix3d M = fiber_tensor<double>(pt.alpha);
    Matrix3d S = compute_stress<double>(kin, xidot, ai, pt.alpha, p);
    const double I4q = (kin.C * M).trace();  // fiber invariant e_F . C . e_F
    double g[6], dgdtau[6], dgdp[6], Kd[6][6], K[6][6];
    // Mixed u/p: d(g[k])/d(p_mix) = (J C^{-1}) : dE[k], since Sigma^b = p_mix J C^{-1}.
    const Matrix3d JCinv =
        use_mixed ? (kin.J * kin.C.inverse()).eval() : Matrix3d::Zero();
    for (int i = 0; i < 6; i++) {
      g[i] = contract<double>(S, kin.dE[i]);
      dgdtau[i] = contract<double>(M, kin.dE[i]);
      dgdp[i] = use_mixed ? contract<double>(JCinv, kin.dE[i]) : 0.0;
      for (int j = 0; j < 6; j++)
        Kd[i][j] = p.gamma * contract<double>(kin.dE[i], kin.dE[j]);
    }
    // Complex-step material + geometric tangent K[k][m] = d g[k] / d xi[m].
    // A tiny imaginary perturbation gives the derivative from Im(g)/h with no
    // subtractive cancellation (rates held fixed; their tangent is Kd above).
    {
      const double h = 1e-30;
      cplx xic[6], xidotc[6], gc[6];
      for (int i = 0; i < 6; i++) {
        xic[i] = xi[i];
        xidotc[i] = xidot[i];
      }
      for (int m = 0; m < 6; m++) {
        xic[m] += cplx(0.0, h);
        compute_g<cplx>(xic, xidotc, ai, pt.R, pt.alpha, p, gc);
        for (int k = 0; k < 6; k++) K[k][m] = std::imag(gc[k]) / h;
        xic[m] = xi[m];
      }
    }

    // Active-DOF coupling of the mechanical forces: simple model couples to the
    // global tau (d sigma_a/d tau = 1); BCS couples to this point's tau_c via
    // d sigma_1D/d tau_c = 1/(1 + e_fib) = 1/sqrt(I4).
    const int active_col = use_bcs ? i_tauc(q) : i_tau();
    const double i4denom = std::pow(I4q, i4pow);  // sigma_1D = T_fib / I4^i4pow
    const double active_coef = use_bcs ? (1.0 / i4denom) : 1.0;

    // Row (equation) and column (variable) contributions of each local force.
    // xi index: 0=rho(N), 1=rho'(dN), 2=phi'(dN), 3=beta, 4=eps, 5=eta'(dN).
    int req[6][2], cvr[6][2], rn[6], cn[6];
    double rw[6][2], cw[6][2];
    req[0][0] = e_rho(A0); rw[0][0] = N0; req[0][1] = e_rho(A1); rw[0][1] = N1; rn[0] = 2;
    cvr[0][0] = i_rho(A0); cw[0][0] = N0; cvr[0][1] = i_rho(A1); cw[0][1] = N1; cn[0] = 2;
    req[1][0] = e_rho(A0); rw[1][0] = dN0; req[1][1] = e_rho(A1); rw[1][1] = dN1; rn[1] = 2;
    cvr[1][0] = i_rho(A0); cw[1][0] = dN0; cvr[1][1] = i_rho(A1); cw[1][1] = dN1; cn[1] = 2;
    {  // phi': rows skip the pinned inner node
      int c = 0;
      if (A0 != 0) { req[2][c] = e_phi(A0); rw[2][c] = dN0; c++; }
      req[2][c] = e_phi(A1); rw[2][c] = dN1; c++; rn[2] = c;
    }
    cvr[2][0] = i_phi(A0); cw[2][0] = dN0; cvr[2][1] = i_phi(A1); cw[2][1] = dN1; cn[2] = 2;
    req[3][0] = e_beta(); rw[3][0] = 1.0; rn[3] = 1;
    cvr[3][0] = i_beta(); cw[3][0] = 1.0; cn[3] = 1;
    req[4][0] = e_eps(); rw[4][0] = 1.0; rn[4] = 1;
    cvr[4][0] = i_eps(); cw[4][0] = 1.0; cn[4] = 1;
    {  // eta': rows skip the pinned inner node
      int c = 0;
      if (A0 != 0) { req[5][c] = e_eta(A0); rw[5][c] = dN0; c++; }
      req[5][c] = e_eta(A1); rw[5][c] = dN1; c++; rn[5] = c;
    }
    cvr[5][0] = i_eta(A0); cw[5][0] = dN0; cvr[5][1] = i_eta(A1); cw[5][1] = dN1; cn[5] = 2;

    const double pref = 2.0 * M_PI * Lp * pt.w;
    for (int k = 0; k < 6; k++) {
      for (int ri = 0; ri < rn[k]; ri++) {
        const int re = req[k][ri];
        const double rwt = pref * rw[k][ri];
        Cloc[re] += rwt * g[k];
        Kat(re, active_col) += rwt * active_coef * dgdtau[k];
        // Mixed u/p: coupling of the equilibrium to this element's pressure DOF.
        if (use_mixed) Kat(re, i_pmix(pt.elem)) += rwt * dgdp[k];
        // BCS: sigma_1D also depends on e_c_dot via the mu*e_c_dot term.
        if (use_bcs)
          Kdat(re, i_ec(q)) += rwt * (bp.mu / i4denom) * dgdtau[k];
        for (int m = 0; m < 6; m++) {
          for (int ci = 0; ci < cn[m]; ci++) {
            const int cv = cvr[m][ci];
            Kat(re, cv) += rwt * K[k][m] * cw[m][ci];
            Kdat(re, cv) += rwt * Kd[k][m] * cw[m][ci];
          }
        }
      }
    }

    // --- Mixed u/p: weak incompressibility constraint for this element (P0) ---
    // Residual C[e_pmix] = sum_{q in e} pref (J-1) = 0. J depends only on
    // rho, rho', eps; map the analytic dJ/dxi through the same columns as the
    // force assembly. The pressure-pressure block is zero (saddle-point form).
    if (use_mixed) {
      const int erow = e_pmix(pt.elem);
      Cloc[erow] += pref * (kin.J - 1.0);
      const double Rq = pt.R;
      double dJ[6] = {0, 0, 0, 0, 0, 0};
      dJ[0] = (1.0 + xi[1]) * (1.0 / Rq) * (1.0 + xi[4]);  // dJ/drho
      dJ[1] = (1.0 + xi[0] / Rq) * (1.0 + xi[4]);          // dJ/drho'
      dJ[4] = (1.0 + xi[1]) * (1.0 + xi[0] / Rq);          // dJ/deps
      for (int m = 0; m < 6; m++) {
        if (dJ[m] == 0.0) continue;
        for (int ci = 0; ci < cn[m]; ci++)
          Kat(erow, cvr[m][ci]) += pref * dJ[m] * cw[m][ci];
      }
    }

    // --- BCS internal-variable ODEs at this quadrature point (Eqs. 32-33) ---
    if (use_bcs) {
      const double sfib = std::sqrt(I4q);
      const double e_fib = sfib - 1.0;
      const int erow[3] = {e_ec(q), e_tauc(q), e_kc(q)};
      const int avar[3] = {i_ec(q), i_tauc(q), i_kc(q)};
      double a[3] = {Y(i_ec(q)), Y(i_tauc(q)), Y(i_kc(q))};
      double ad[3] = {Yd(i_ec(q)), Yd(i_tauc(q)), Yd(i_kc(q))};
      double Ra[3];
      bcs_residual<double>(a, ad, e_fib, act, act_plus, bp, Ra);
      for (int i = 0; i < 3; i++) Cloc[erow[i]] += Ra[i];
      // Complex-step tangents w.r.t. the local active state (dC_dy) and its
      // rates (dC_dydot). The regularized |e_c dot| makes the force-velocity
      // decay smooth, so Im(R)/h is the exact tangent even for large bp.alpha.
      {
        const double h = 1e-30;
        cplx ac[3], adc[3], Rc[3];
        for (int t = 0; t < 3; t++) {
          ac[t] = a[t];
          adc[t] = ad[t];
        }
        for (int jj = 0; jj < 3; jj++) {
          ac[jj] += cplx(0.0, h);
          bcs_residual<cplx>(ac, adc, cplx(e_fib), act, act_plus, bp, Rc);
          for (int ii = 0; ii < 3; ii++)
            Kat(erow[ii], avar[jj]) += std::imag(Rc[ii]) / h;
          ac[jj] = a[jj];
          adc[jj] += cplx(0.0, h);
          bcs_residual<cplx>(ac, adc, cplx(e_fib), act, act_plus, bp, Rc);
          for (int ii = 0; ii < 3; ii++)
            Kdat(erow[ii], avar[jj]) += std::imag(Rc[ii]) / h;
          adc[jj] = ad[jj];
        }
      }
      // Only the e_c ODE couples to the wall kinematics, through the fiber
      // strain: dR_ec/de_fib = -k_s/mu, de_fib/dxi_m = dgdtau[m]/sqrt(I4). The
      // tau_c/k_c ODEs couple to the mechanics only via e_c (handled above).
      const double dRec_defib = -bp.k_s / bp.mu;
      for (int m = 0; m < 6; m++) {
        const double coef = dRec_defib * dgdtau[m] / sfib;
        for (int ci = 0; ci < cn[m]; ci++)
          Kat(e_ec(q), cvr[m][ci]) += coef * cw[m][ci];
      }
    }
  }

  // --- External pressure loading on the inner surface and lids (Eq. 63) ---
  const double rho_i = Y(i_rho(0));
  const double Rinner = Rip + rho_i;
  const double twoPiL = 2.0 * M_PI * Lp;
  const double piL = M_PI * Lp;
  // rho equation at inner node: - Pv * 2 pi L * Rinner * (1+eps)
  {
    const int re = e_rho(0);
    Cloc[re] -= Pin * twoPiL * Rinner * (1.0 + eps);
    Kat(re, 0) -= twoPiL * Rinner * (1.0 + eps);          // d/dPin
    Kat(re, i_rho(0)) -= Pin * twoPiL * (1.0 + eps);      // d/drho_i
    Kat(re, i_eps()) -= Pin * twoPiL * Rinner;            // d/deps
  }
  // eps equation: - Pv * pi L * Rinner^2
  {
    const int re = e_eps();
    Cloc[re] -= Pin * piL * Rinner * Rinner;
    Kat(re, 0) -= piL * Rinner * Rinner;                  // d/dPin
    Kat(re, i_rho(0)) -= Pin * piL * 2.0 * Rinner;        // d/drho_i
  }

  // --- Cavity volume definition: V - pi Rinner^2 (1+eps) L = 0 ---
  {
    const int re = e_volume();
    Cloc[re] -= piL * Rinner * Rinner * (1.0 + eps);
    Kat(re, i_rho(0)) -= piL * 2.0 * Rinner * (1.0 + eps);
    Kat(re, i_eps()) -= piL * Rinner * Rinner;
  }

  // --- Mass conservation: Qin - Qout - Vdot = 0, with the geometric rate ---
  {
    const int re = e_mass();
    const double rho_i_d = Yd(i_rho(0));
    const double Vdot = piL * (2.0 * Rinner * (1.0 + eps) * rho_i_d +
                               Rinner * Rinner * eps_d);
    Cloc[re] -= Vdot;
    // d/dy
    Kat(re, i_rho(0)) -=
        piL * (2.0 * (1.0 + eps) * rho_i_d + 2.0 * Rinner * eps_d);
    Kat(re, i_eps()) -= piL * 2.0 * Rinner * rho_i_d;
    // d/dydot
    Kdat(re, i_rho(0)) -= piL * 2.0 * Rinner * (1.0 + eps);
    Kdat(re, i_eps()) -= piL * Rinner * Rinner;
  }

  // --- Simple active stress ODE: -sigma_max * a_+ (linear parts in E/F). The
  // BCS ODEs are assembled per quadrature point in the loop above. ---
  if (!use_bcs) Cloc[e_active()] -= sigma_max * act_plus;

  // --- Scatter local contributions into the sparse system (assign) ---
  for (int i = 0; i < n_eqn; i++) {
    system.C(global_eqn_ids[i]) = Cloc[i];
    for (int j = 0; j < n_var; j++) {
      system.dC_dy.coeffRef(global_eqn_ids[i], global_var_ids[j]) = Kat(i, j);
      system.dC_dydot.coeffRef(global_eqn_ids[i], global_var_ids[j]) = Kdat(i, j);
    }
  }
}

void ChamberCylinder::get_activation(std::vector<double>& parameters) {
  const double alpha_max = parameters[global_param_ids[ParamId::alpha_max]];
  const double alpha_min = parameters[global_param_ids[ParamId::alpha_min]];
  const double tsys = parameters[global_param_ids[ParamId::tsys]];
  const double tdias = parameters[global_param_ids[ParamId::tdias]];
  const double steepness = parameters[global_param_ids[ParamId::steepness]];

  const auto T_cardiac = model->cardiac_cycle_period;
  const auto t_in_cycle = fmod(model->time, T_cardiac);

  auto warp_signed = [T_cardiac](double dt) {
    return fmod(dt + 1.5 * T_cardiac, T_cardiac) - 0.5 * T_cardiac;
  };

  const double phase_tsys = warp_signed(t_in_cycle - tsys);
  const double phase_tdias = warp_signed(t_in_cycle - tdias);

  const double S_plus = 0.5 * (1.0 + tanh(phase_tsys / steepness));
  const double S_minus = 0.5 * (1.0 - tanh(phase_tdias / steepness));
  const double f = S_plus * S_minus;

  const double act_t = alpha_max * f + alpha_min * (1.0 - f);
  act = std::abs(act_t);
  act_plus = std::max(act_t, 0.0);
}
