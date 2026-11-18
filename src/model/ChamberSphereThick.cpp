// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberSphereThick.h"

#include <algorithm>
#include <cmath>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>

#include "Model.h"

namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

/// Double contraction A:B.
inline double contract(const Matrix3d& A, const Matrix3d& B) {
  return (A.array() * B.array()).sum();
}

/**
 * @brief Kinematics of the incompressible spherically symmetric wall.
 *
 * At the reference radius R, the exact incompressible motion
 * r^3 = R^3 + r_i^3 - R_i^3 gives the circumferential stretch lambda = r/R and
 * the derivatives of r (and hence of lambda = r/R) with respect to the single
 * unknown rho = r_i - R_i.
 */
struct Kinematics {
  double lam;     ///< Circumferential stretch r/R
  double lam_p;   ///< d(lambda)/d(rho)
  double lam_pp;  ///< d2(lambda)/d(rho)2
  double r_p;     ///< d(r)/d(rho)
  double r_pp;    ///< d2(r)/d(rho)2
  double r_ppp;   ///< d3(r)/d(rho)3
  Matrix3d C;     ///< Right Cauchy-Green tensor diag(l^-4, l^2, l^2)
  Matrix3d A;     ///< dE/d(lambda) = diag(-2 l^-5, l, l)
  Matrix3d B;     ///< d2E/d(lambda)2 = diag(10 l^-6, 1, 1)
};

/// Build the kinematics at reference radius R for the endocardial radius ri
/// of a shell whose reference endocardial radius is Ri.
Kinematics compute_kinematics(double R, double ri, double Ri) {
  const double u = ri * ri * ri - Ri * Ri * Ri;
  const double r = std::cbrt(R * R * R + u);
  const double r2 = r * r;
  const double ri2 = ri * ri;

  Kinematics k;
  k.lam = r / R;
  // r_,rho = ri^2 / r^2 from 3 r^2 r_,rho = du/drho = 3 ri^2
  k.r_p = ri2 / r2;
  k.r_pp = 2.0 * ri / r2 - 2.0 * ri2 * ri2 / (r2 * r2 * r);
  k.r_ppp = 10.0 * ri2 * ri2 * ri2 / (r2 * r2 * r2 * r2) -
            12.0 * ri2 * ri / (r2 * r2 * r) + 2.0 / r2;
  k.lam_p = k.r_p / R;
  k.lam_pp = k.r_pp / R;

  const double l = k.lam;
  const double l2 = l * l;
  k.C = Vector3d(1.0 / (l2 * l2), l2, l2).asDiagonal();
  k.A = Vector3d(-2.0 / (l2 * l2 * l), l, l).asDiagonal();
  k.B = Vector3d(10.0 / (l2 * l2 * l2), 1.0, 1.0).asDiagonal();
  return k;
}

}  // namespace

void ChamberSphereThick::setup_dofs(DOFHandler& dofhandler) {
  check_inputs();

  const double ri = model->get_parameter_value(global_param_ids[ParamId::Ri]);
  const double re = model->get_parameter_value(global_param_ids[ParamId::Re]);
  const int n_ele = std::max(1, (int)std::lround(model->get_parameter_value(
                                    global_param_ids[ParamId::num_elements])));
  n_law = active_stress_->num_states();

  // Composite three-point Gauss quadrature over [Ri, Re], from the
  // endocardium (q = 0) to the epicardium. The weight carries the R^2 of the
  // spherical volume element 4 pi R^2 dR; the 4 pi is a prefactor below.
  const double gp[3] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
  const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  const double h = (re - ri) / n_ele;
  quad.clear();
  for (int e = 0; e < n_ele; e++) {
    const double Ra = ri + e * h, Rb = ri + (e + 1) * h;
    for (int q = 0; q < 3; q++) {
      QuadPoint pt;
      pt.R = 0.5 * (Ra + Rb) + 0.5 * (Rb - Ra) * gp[q];
      pt.w = gw[q] * 0.5 * (Rb - Ra) * pt.R * pt.R;
      quad.push_back(pt);
    }
  }
  n_quad = (int)quad.size();

  // Internal variables: endocardial displacement, the active stress states of
  // every quadrature point, the cavity volume and the velocity companion.
  std::list<std::string> int_names;
  int_names.push_back("radius");
  const std::vector<std::string> law_states = active_stress_->state_names();
  for (int q = 0; q < n_quad; q++) {
    for (const auto& name : law_states)
      int_names.push_back(name + "_" + std::to_string(q));
  }
  int_names.push_back("volume");
  int_names.push_back("velo");

  // wall + active stress + volume + mass + pressure + velocity companion
  n_eqn = 1 + n_active_var() + 4;
  Block::setup_dofs_(dofhandler, n_eqn, int_names);

  // Reservation estimates
  num_triplets.F = 6;
  num_triplets.E = 1;
  num_triplets.D = 5 + n_quad + 2 * n_law * n_quad * (n_law + 1);
}

void ChamberSphereThick::update_constant(SparseSystem& system,
                                         std::vector<double>& parameters) {
  (void)parameters;

  // Cavity volume V - 4/3 pi r_i^3 = 0 (linear part)
  system.F.coeffRef(global_eqn_ids[e_volume()], global_var_ids[i_vol()]) = 1.0;

  // Mass conservation Qin - Qout - Vdot = 0, with Vdot = 4 pi r_i^2 v from
  // the velocity companion (update_solution). The volume unknown is a pure
  // output and does not feed back, so an inconsistent initial volume cannot
  // drive a spurious flow.
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[3]) = -1.0;

  // Pressure equality Pin - Pout = 0
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[2]) = -1.0;

  // Velocity companion v - d(rho)/dt = 0
  system.F.coeffRef(global_eqn_ids[e_velo()], global_var_ids[i_velo()]) = 1.0;
  system.E.coeffRef(global_eqn_ids[e_velo()], global_var_ids[i_rho()]) = -1.0;
}

void ChamberSphereThick::update_time(SparseSystem& system,
                                     std::vector<double>& parameters) {
  (void)system;
  (void)parameters;
  active_stress_->update_time(activation_function_->compute(model->time));
}

void ChamberSphereThick::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  const double Ri_ = parameters[global_param_ids[ParamId::Ri]];
  const double gam = parameters[global_param_ids[ParamId::gamma]];
  const double dens = parameters[global_param_ids[ParamId::density]];

  auto Y = [&](int lv) { return y[global_var_ids[lv]]; };
  auto Yd = [&](int lv) { return dy[global_var_ids[lv]]; };

  const double rho = Y(i_rho());
  const double v = Y(i_velo());
  const double vdot = Yd(i_velo());
  const double Pin = Y(0);
  const double ri = Ri_ + rho;

  // The fiber direction is immaterial as long as it lies in the wall plane:
  // C has equal in-plane eigenvalues, so I4 = lambda^2 for any such direction.
  const Material::Vec3 f(0.0, 1.0, 0.0);
  // In-plane isotropic average <e_F (x) e_F> of the active fiber stress
  const Matrix3d Pact = Vector3d(0.0, 0.5, 0.5).asDiagonal();

  // --- Wall equation: inertia + internal virtual work over the thickness ---
  double C_wall = 0.0, dCw_drho = 0.0, dCw_dv = 0.0, dCw_dvdot = 0.0;
  std::vector<double> a(n_law), ad(n_law);
  for (int q = 0; q < n_quad; q++) {
    const QuadPoint& pt = quad[q];
    const Kinematics k = compute_kinematics(pt.R, ri, Ri_);
    const double tau = Y(i_law(q, 0));

    // Sigma = Sigma^d + gamma Edot + (tau/2)(e_Th (x) e_Th + e_Ph (x) e_Ph)
    const Matrix3d dC = 2.0 * k.A;  // dC/d(lambda)
    const double dJ = 0.0;          // the motion is isochoric
    Matrix3d dSd;
    const Matrix3d Sd =
        material_->stress_derivatives(k.C, 1.0, f, 1, &dC, &dJ, &dSd);
    const double lamdot = k.lam_p * v;
    const Matrix3d S = Sd + (gam * lamdot) * k.A + tau * Pact;

    // Internal force g = Sigma : dE/d(lambda) d(lambda)/d(rho)
    const double G = contract(S, k.A);
    const double g = G * k.lam_p;
    const Matrix3d dS_drho =
        dSd * k.lam_p +
        (gam * v) * (k.B * (k.lam_p * k.lam_p) + k.A * k.lam_pp);
    const double dG_drho = contract(dS_drho, k.A) + contract(S, k.B) * k.lam_p;
    const double dg_drho = dG_drho * k.lam_p + G * k.lam_pp;
    const double dg_dv = gam * k.lam_p * k.lam_p * contract(k.A, k.A);
    const double dg_dtau = k.lam * k.lam_p;  // Pact : A = lambda

    // Inertia M(rho) vdot + 1/2 M'(rho) v^2 of the radial motion
    const double iner = dens * (k.r_p * k.r_p * vdot + k.r_p * k.r_pp * v * v);
    const double diner_drho =
        dens * (2.0 * k.r_p * k.r_pp * vdot +
                (k.r_pp * k.r_pp + k.r_p * k.r_ppp) * v * v);
    const double diner_dv = dens * 2.0 * k.r_p * k.r_pp * v;
    const double diner_dvdot = dens * k.r_p * k.r_p;

    const double pref = 4.0 * M_PI * pt.w;
    C_wall += pref * (iner + g);
    dCw_drho += pref * (diner_drho + dg_drho);
    dCw_dv += pref * (diner_dv + dg_dv);
    dCw_dvdot += pref * diner_dvdot;
    system.dC_dy.coeffRef(global_eqn_ids[e_wall()],
                          global_var_ids[i_law(q, 0)]) = pref * dg_dtau;

    // --- Active stress law at this quadrature point ---
    for (int i = 0; i < n_law; i++) {
      a[i] = Y(i_law(q, i));
      ad[i] = Yd(i_law(q, i));
    }
    ActiveStress::Linearization lin;
    active_stress_->linearize(a.data(), ad.data(), k.lam, lin);
    for (int i = 0; i < n_law; i++) {
      const int erow = global_eqn_ids[e_law(q, i)];
      system.C.coeffRef(erow) = lin.R[i];
      for (int j = 0; j < n_law; j++) {
        system.dC_dy.coeffRef(erow, global_var_ids[i_law(q, j)]) =
            lin.dR_da[i * n_law + j];
        system.dC_dydot.coeffRef(erow, global_var_ids[i_law(q, j)]) =
            lin.dR_dadot[i * n_law + j];
      }
      // The law depends on the wall through the fiber stretch lambda only.
      system.dC_dy.coeffRef(erow, global_var_ids[i_rho()]) =
          lin.dR_dlambda[i] * k.lam_p;
    }
  }

  // --- Cavity pressure on the endocardium: - P_v 4 pi r_i^2 ---
  C_wall -= Pin * 4.0 * M_PI * ri * ri;
  dCw_drho -= Pin * 8.0 * M_PI * ri;
  system.dC_dy.coeffRef(global_eqn_ids[e_wall()], global_var_ids[0]) =
      -4.0 * M_PI * ri * ri;

  system.C.coeffRef(global_eqn_ids[e_wall()]) = C_wall;
  system.dC_dy.coeffRef(global_eqn_ids[e_wall()], global_var_ids[i_rho()]) =
      dCw_drho;
  system.dC_dy.coeffRef(global_eqn_ids[e_wall()], global_var_ids[i_velo()]) =
      dCw_dv;
  system.dC_dydot.coeffRef(global_eqn_ids[e_wall()], global_var_ids[i_velo()]) =
      dCw_dvdot;

  // --- Cavity volume V - 4/3 pi r_i^3 = 0 ---
  system.C.coeffRef(global_eqn_ids[e_volume()]) =
      -4.0 / 3.0 * M_PI * ri * ri * ri;
  system.dC_dy.coeffRef(global_eqn_ids[e_volume()], global_var_ids[i_rho()]) =
      -4.0 * M_PI * ri * ri;

  // --- Mass conservation Qin - Qout - Vdot = 0 with Vdot = 4 pi r_i^2 v ---
  system.C.coeffRef(global_eqn_ids[e_mass()]) = -4.0 * M_PI * ri * ri * v;
  system.dC_dy.coeffRef(global_eqn_ids[e_mass()], global_var_ids[i_rho()]) =
      -8.0 * M_PI * ri * v;
  system.dC_dy.coeffRef(global_eqn_ids[e_mass()], global_var_ids[i_velo()]) =
      -4.0 * M_PI * ri * ri;
}

void ChamberSphereThick::check_inputs() const {
  const std::string where = "ChamberSphereThick '" + get_name() + "': ";
  const double ri = model->get_parameter_value(global_param_ids[ParamId::Ri]);
  const double re = model->get_parameter_value(global_param_ids[ParamId::Re]);
  if (ri <= 0.0) throw std::runtime_error(where + "Ri must be positive.");
  if (re <= ri) throw std::runtime_error(where + "Re must be larger than Ri.");
  if (model->get_parameter_value(global_param_ids[ParamId::density]) < 0.0)
    throw std::runtime_error(where + "density must not be negative.");
  if (!material_)
    throw std::runtime_error(
        where + "material not set. Provide a \"material\" JSON block.");
  if (!active_stress_)
    throw std::runtime_error(
        where +
        "active_stress not set. Provide an \"active_stress\" JSON block.");
  if (!activation_function_)
    throw std::runtime_error(where +
                             "activation_function not set. Provide an "
                             "\"activation_function\" JSON block.");
  if (activation_function_->signal_type() !=
      active_stress_->required_signal_type())
    throw std::runtime_error(
        where + "active_stress law requires a " +
        (active_stress_->required_signal_type() == ActivationSignalType::rate
             ? "rate-type (e.g. piecewise_rate)"
             : "normalized (e.g. half_cosine)") +
        " activation function.");
}
