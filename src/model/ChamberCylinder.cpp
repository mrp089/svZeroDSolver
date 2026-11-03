// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ChamberCylinder.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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
  double C1, C2, C3, C4, C5, C6, gamma;
  double p_mix = 0.0;    ///< mixed u/p element hydrostatic pressure
  /// sigma_1D = T_fib / (1+e_fib)^stretch_pow: 1 = genet23 Eq 30; 2 = the
  /// dt -> 0 limit of Eq 59 as printed (T_fib / (1 + 2 e_F.E.e_F)).
  double stretch_pow = 1.0;
  /// Viscous law: 0 = W_v = gamma/2 |E_dot|^2 (Eq 28); 1 = W_v = gamma/2
  /// (tr E_dot)^2, the literal reading of Chapelle 2012 / Caruel 2013; 2 =
  /// W_v = 1/2 sum_ab gam_ab (E_dot_ab)^2 in the fiber basis (f, r, n)
  /// (Kelvin-Voigt surrogate of Nordsletten 2021).
  int visc_mode = 0;
  Eigen::Matrix3d gam = Eigen::Matrix3d::Zero();
  /// Passive law: 0 = genet23 Eq 25; 1 = Fung W = C/2 (exp Q - 1),
  /// Q = 2 blin tr(Eb) + sum_ab b_ab Eb_ab^2 in the fiber basis (f, r, n), Eb
  /// the isochoric Green-Lagrange strain; 2 = Holzapfel-Ogden (HO 2009 with
  /// tension-only f and s terms; I1 isochoric, I4/I8 not).
  int passive_law = 0;
  double fung_C = 0.0, fung_blin = 0.0;
  Eigen::Matrix3d fung_b = Eigen::Matrix3d::Zero();
  double ho_a = 0.0, ho_b = 1.0, ho_af = 0.0, ho_bf = 1.0, ho_as = 0.0,
         ho_bs = 1.0, ho_afs = 0.0, ho_bfs = 1.0;
  int ho_sheet = 0;  ///< 0 = transmural e_R, 1 = in-plane e_F x e_R
};

/// Orthonormal fiber basis as columns (e_F, e_R, e_N), e_N = e_F x e_R the
/// in-plane cross-fiber direction, in the (e_R, e_Th, e_Z) frame.
inline Eigen::Matrix3d fiber_basis(double alpha) {
  const double ca = std::cos(alpha), sa = std::sin(alpha);
  Eigen::Matrix3d Q;
  Q << 0.0, 1.0, 0.0,
       ca, 0.0, sa,
       sa, 0.0, -ca;
  return Q;
}

/// Active-stress input at an integration point (Genet BCS): bond stress `tauc`
/// and contractile-element strain rate `ecdot`, giving
/// sigma_1D = (tscale*tauc + mu*ecdot)/sqrt(I4).
struct ActiveInput {
  double tauc = 0.0;   ///< BCS active bond stress tau_c
  double ecdot = 0.0;  ///< BCS contractile-element strain rate
  double mu = 0.0;     ///< BCS active dissipation
  /// Factor on tau_c in the fiber tension (Sainte-Marie 2006: the Frank-
  /// Starling factor d(e_c) scales the output stress); 1 = genet23.
  double tscale = 1.0;
  /// Series law: 0 = additive in engineering strain (genet23 Eqs 30-32,
  /// Kimmig 2019); 1 = multiplicative in Green-Lagrange strain (Sainte-Marie
  /// 2006, Chapelle 2012, Caruel 2013): sigma_1D = (1+2e_c)/(1+2e_1D) T_fib.
  int series = 0;
  double ec = 0.0;     ///< contractile strain e_c (needed for series = 1)
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

/// Fung-type deviatoric 2nd-PK stress (MatParams::passive_law = 1):
/// W = C/2 (exp Q - 1), Q = 2 blin tr(Eb) + sum_ab b_ab Eb_ab^2 with Eb =
/// (J^-2/3 C - I)/2 in the fiber basis (off-diagonal pairs counted twice, i.e.
/// 2 b_ab Eb_ab^2). Sigma = J^-2/3 (Sb - (Sb:C)/3 C^-1), Sb = dW/dEb.
template <typename T>
Mat3<T> passive_fung(const Mat3<T>& C, const T& J, const Mat3<T>& Cinv,
                     double alpha, const MatParams& p) {
  const Mat3<T> Qb = fiber_basis(alpha).template cast<T>();
  const Mat3<T> I = Mat3<T>::Identity();
  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const Mat3<T> Eb = 0.5 * (Jm23 * C - I);
  const Mat3<T> Ep = Qb.transpose() * Eb * Qb;
  const Mat3<T> Bp = p.fung_b.template cast<T>();
  const T Q = 2.0 * p.fung_blin * Ep.trace() + (Bp.array() * Ep.array() * Ep.array()).sum();
  const Mat3<T> dQ = 2.0 * p.fung_blin * I + 2.0 * (Bp.array() * Ep.array()).matrix();
  const Mat3<T> Sbar = (0.5 * p.fung_C * std::exp(clamp_max(Q, 300.0))) *
                       (Qb * dQ * Qb.transpose());
  return Jm23 * (Sbar - (contract<T>(Sbar, C) / 3.0) * Cinv);
}

/// Holzapfel-Ogden deviatoric 2nd-PK stress (MatParams::passive_law = 2):
/// W = a/2b (e^{b(I1b-3)} - 1) + sum_{i=f,s} a_i/2b_i (e^{b_i(I4i-1)^2} - 1)
/// [I4i > 1] + a_fs/2b_fs (e^{b_fs I8^2} - 1).
template <typename T>
Mat3<T> passive_ho(const Mat3<T>& C, const T& J, const Mat3<T>& Cinv,
                   double alpha, const MatParams& p) {
  const Eigen::Matrix3d Qb = fiber_basis(alpha);
  const Vec3<T> f = Qb.col(0).template cast<T>();
  const Vec3<T> s = Qb.col(p.ho_sheet == 1 ? 2 : 1).template cast<T>();
  const Mat3<T> I = Mat3<T>::Identity();
  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const T I1 = C.trace();
  Mat3<T> S = (p.ho_a * std::exp(clamp_max(p.ho_b * (Jm23 * I1 - 3.0), 300.0))) *
              (Jm23 * (I - (I1 / 3.0) * Cinv));
  const T I4f = f.dot(C * f), I4s = s.dot(C * s), I8 = f.dot(C * s);
  if (re(I4f) > 1.0)
    S += (2.0 * p.ho_af * (I4f - 1.0) *
          std::exp(clamp_max(p.ho_bf * (I4f - 1.0) * (I4f - 1.0), 300.0))) * (f * f.transpose());
  if (re(I4s) > 1.0)
    S += (2.0 * p.ho_as * (I4s - 1.0) *
          std::exp(clamp_max(p.ho_bs * (I4s - 1.0) * (I4s - 1.0), 300.0))) * (s * s.transpose());
  S += (p.ho_afs * I8 * std::exp(clamp_max(p.ho_bfs * I8 * I8, 300.0))) *
       (f * s.transpose() + s * f.transpose());
  return S;
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
  if (p.passive_law == 1) Sd = passive_fung<T>(C, J, Cinv, alpha, p);
  else if (p.passive_law == 2) Sd = passive_ho<T>(C, J, Cinv, alpha, p);
  // Bulk (incompressibility) 2nd-PK stress Sigma^b = p_mix J C^{-1}, with p_mix
  // the element hydrostatic pressure (mixed u/p; held fixed under the xi
  // complex-step, its column is added analytically).
  const T bulk = T(p.p_mix);
  Mat3<T> Sb = bulk * J * Cinv;                            // bulk

  Mat3<T> Edot = Mat3<T>::Zero();                           // viscous
  for (int m = 0; m < 6; m++) Edot += k.dE[m] * xidot[m];
  Mat3<T> Sv = p.visc_mode == 1 ? Mat3<T>(p.gamma * Edot.trace() * I)
                                : Mat3<T>(p.gamma * Edot);
  if (p.visc_mode == 2) {  // Sigma^v_ab = gam_ab E_dot_ab in the fiber basis
    const Mat3<T> Qb = fiber_basis(alpha).template cast<T>();
    const Mat3<T> Ep = Qb.transpose() * Edot * Qb;
    const Mat3<T> Sp = (Ep.array() * p.gam.template cast<T>().array()).matrix();
    Sv = Qb * Sp * Qb.transpose();
  }

  // active fiber stress magnitude:
  // sigma_1D = (tau_c + mu * e_c_dot) / (1 + e_fib)^pow, 1 + e_fib = sqrt(I4),
  // pow = 1 (Eq 30) unless the Eq 59 diagnostic is selected.
  // Written via the fiber tension tau_c + mu*e_c_dot (= T_fib) so the stiff
  // series stiffness k_s stays out of the mechanical residual.
  const T Tfib = act.tscale * act.tauc + act.mu * act.ecdot;
  const T sigma_act = act.series == 1
                          ? T(Tfib * (1.0 + 2.0 * act.ec) / I4)
                          : T(Tfib / std::pow(I4, 0.5 * p.stretch_pow));
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
/// piecewise-linear length-dependence curve (n_0 = interp(e_c, abscissas,
/// ordinates); Caruel et al. 2013, \cite genet23 Ref. 13). It represents the
/// actin-myosin overlap (length-tension) relationship: n_0 rises from 0,
/// plateaus at 1 over the physiological contractile strain e_c in [0.20, 0.47],
/// and falls back to 0. The default breakpoints are the MEDISIM PhysioBlocks
/// reference values (see ChamberCylinder::n0_ec_/n0_val_), overridable per block.
/// NOTE: the plateau sits at e_c ~ 0.2, not ~1 as Caruel Fig. 7(a) suggests -
/// that figure is drawn over the wide isotonic papillary-muscle strain range.
/// X, Y are the (double) breakpoints; e_c may be complex (complex-step).
template <typename T>
inline T frank_starling(T e_c, const std::vector<double>& X,
                        const std::vector<double>& Y) {
  const int NP = static_cast<int>(X.size());
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
  double k_s, k_0, mu, alpha, sigma0;
  double fv_reg = 1e-4;  ///< force-velocity regularizer eps (|e_c dot| smoothing)
  /// alpha in the tau_c decay only (k_c keeps `alpha`). Equal to `alpha` for
  /// Eq 33; (alpha + 1)/2 reproduces Chapelle 2012 Eq 26 as printed, whose
  /// tau_c/sqrt(k_c) update omits alpha in its decay term.
  double alpha_tau = 0.0;
  /// Series law (see ActiveInput::series): 0 = additive engineering strain,
  /// e_kin = e_fib; 1 = multiplicative Green-Lagrange, e_kin = e_1D = e_F.E.e_F.
  int series = 0;
  /// Frank-Starling placement: false = n0(e_c) on the sources (genet23 Eq 33);
  /// true = d(e_c) on the output tension, sources without n0 (Sainte-Marie 2006).
  bool starling_out = false;
};

/// Residuals of the BCS internal-variable ODEs (Genet Eqs. 32-33 / Chapelle
/// Eq. 9) at a quadrature point. a = [e_c, tau_c, k_c], ad = its rate; the
/// contractile-element strain rate is e_c_dot = ad[0]. e_kin is the fiber
/// strain of the series law (e_fib or e_1D, see BcsParams::series). n0_x/n0_y
/// are the Frank-Starling force-length breakpoints passed to frank_starling().
template <typename T>
void bcs_residual(const T a[], const T ad[], T e_kin,
                  double nu_abs, double nu_plus, const BcsParams& bp,
                  const std::vector<double>& n0_x,
                  const std::vector<double>& n0_y, T res[]) {
  const T ecdot = ad[0];
  // Smoothed |e_c dot| (eps = bp.fv_reg regularizes the force-velocity kink so
  // Newton and its complex-step tangent behave; negligible vs |nu| physically).
  const T abs_ecd = std::sqrt(ecdot * ecdot + bp.fv_reg);
  const T n0 = frank_starling<T>(a[0], n0_x, n0_y);
  const T n0_src = bp.starling_out ? T(1.0) : n0;       // on the sources
  const T tau_eff = bp.starling_out ? T(n0 * a[1]) : a[1];  // in the tension
  const T decay = nu_abs + bp.alpha * abs_ecd;
  const T decay_tau = nu_abs + bp.alpha_tau * abs_ecd;
  // series-spring force balance, mu * e_c_dot = (series force) - tau_eff:
  //   additive (Eq 32):   series force = k_s (e_fib - e_c)
  //   multiplicative GL:  series force = E_s (e_1D - e_c)(1+2e_1D)/(1+2e_c)^3
  const T fs = bp.series == 1
                   ? T(bp.k_s * (e_kin - a[0]) * (1.0 + 2.0 * e_kin) /
                       ((1.0 + 2.0 * a[0]) * (1.0 + 2.0 * a[0]) * (1.0 + 2.0 * a[0])))
                   : T(bp.k_s * (e_kin - a[0]));
  res[0] = ecdot - (fs - tau_eff) / bp.mu;
  // tau_c dot = -decay*tau_c + e_c_dot*k_c + n0*sigma0*nu_+
  res[1] = ad[1] - (-decay_tau * a[1] + ecdot * a[2] + n0_src * bp.sigma0 * nu_plus);
  // k_c dot = -decay*k_c + n0*k_0*nu_+
  res[2] = ad[2] - (-decay * a[2] + n0_src * bp.k_0 * nu_plus);
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
  // Rigid-body pins: by default only in the static setting (no inertia); see
  // pin_rigid_. rigid_pins = 0 / 1 forces them off / on.
  const double pins =
      model->get_parameter_value(global_param_ids[ParamId::rigid_pins]);
  pin_rigid_ =
      pins < 0.0
          ? model->get_parameter_value(global_param_ids[ParamId::density]) == 0.0
          : pins > 0.5;

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
  // Genet BCS active block: e_c, tau_c, k_c per quadrature point.
  for (int q = 0; q < n_quad; q++) {
    int_names.push_back("ec_" + std::to_string(q));
    int_names.push_back("tauc_" + std::to_string(q));
    int_names.push_back("kc_" + std::to_string(q));
  }
  int_names.push_back("volume");

  // Full dynamics: velocity companion fields w = d(field)/dt (genet23 Eqs. 8,
  // 18, 45, A5-A6). Appended after the volume DOF.
  for (int a = 0; a < n_node; a++) int_names.push_back("vrho_" + std::to_string(a));
  for (int a = 0; a < n_node; a++) int_names.push_back("vphi_" + std::to_string(a));
  for (int a = 0; a < n_node; a++) int_names.push_back("veta_" + std::to_string(a));
  int_names.push_back("vbeta");
  int_names.push_back("veps");

  // Mixed u/p: one element-wise-constant (P0) hydrostatic pressure per element,
  // appended last so the velocity layout above is unchanged.
  for (int e = 0; e < n_ele; e++)
    int_names.push_back("pmix_" + std::to_string(e));

  // 3*n_node field eqns + beta + eps + active eqns + volume + mass + pressure
  // + velocity companion eqns + element incompressibility constraints.
  const int n_vel = 3 * n_node + 2;
  const int n_pmix = n_ele;
  n_eqn = 3 * n_node + 2 + n_active_var() + 3 + n_vel + n_pmix;
  Block::setup_dofs_(dofhandler, n_eqn, int_names);
  n_var = 4 + 3 * n_node + 2 + n_active_var() + 1 + n_vel + n_pmix;

  // Dense block for dC_dy and dC_dydot; a few linear entries for F/E. (These are
  // reservation estimates only; see SparseSystem::reserve.)
  num_triplets.D = n_eqn * n_var;
  num_triplets.F = 6 * n_node + 20;
  num_triplets.E = 21 * n_node + 8;  // consistent mass + velocity companions
}

double ChamberCylinder::inertia_weight(
    const QuadPoint& pt, const std::vector<double>& parameters) const {
  // pt.w = Gauss weight * dR/dxi * R already carries the cylindrical metric R.
  const double metric = parameters[global_param_ids[ParamId::inertia_metric]];
  return metric > 0.5 ? pt.w : pt.w / pt.R;
}

void ChamberCylinder::update_constant(SparseSystem& system,
                                      std::vector<double>& parameters) {
  // Mass conservation: Qin - Qout - Vdot = 0 (Vdot handled in update_solution).
  // The BCS model puts all of its active rate terms in the nonlinear residual C,
  // so it needs no E/F entries here.
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[1]) = 1.0;   // Qin
  system.F.coeffRef(global_eqn_ids[e_mass()], global_var_ids[3]) = -1.0;  // Qout

  // Pressure equality: Pin - Pout = 0.
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[0]) = 1.0;
  system.F.coeffRef(global_eqn_ids[e_pressure()], global_var_ids[2]) = -1.0;

  // Cavity volume definition: V - pi (R_i + rho_i)^2 (1+eps) L = 0 (linear V).
  system.F.coeffRef(global_eqn_ids[e_volume()], global_var_ids[i_vol()]) = 1.0;

  // Rigid-body pins phi(R_i) = 0 and eta(R_i) = 0 (inner node, index 0), only
  // in the static setting (genet23 Sec 2.4: "These constraints are not
  // necessary in the dynamic setting"). With inertia the rigid rotation and
  // axial translation are free; momentum conservation fixes them.
  if (pin_rigid_) {
    system.F.coeffRef(global_eqn_ids[e_phi(0)], global_var_ids[i_phi(0)]) = 1.0;
    system.F.coeffRef(global_eqn_ids[e_eta(0)], global_var_ids[i_eta(0)]) = 1.0;
  }

  // --- Full dynamics: consistent mass matrix + velocity companion equations ---
  // (genet23 Eqs. 8,18,45; App. A.2 Eqs. A5-A6). The velocity field is the
  // push-forward v = Du[zeta] zeta_dot; introducing companion DOFs w = zeta_dot
  // turns the 2nd-order system into 1st order. The inertia force on field k is
  // sum_j M_kj w_dot_j with the consistent mass M = int rho0 (Du)^T Du dOmega.
  // Du (Eq. A5) is orthogonal between the radial/azimuthal (rho,phi,beta) and
  // axial (eta,eps) groups, so after integrating Theta (2 pi) and Z the only
  // couplings are rho-rho, eta-eta, phi-phi, eps-eps, beta-beta, phi-beta and
  // eta-eps. The reference-configuration (R) part of the mass is assembled here;
  // the state-dependent current-config correction R^2 -> (R+rho)^2 for the
  // beta/phi moment arms (genet23 Eq A5) is added as a nonlinear inertia force in
  // update_solution (see "Current-configuration inertia correction"). Inertia is
  // ~1e-4 of the peak internal/pressure forces but the twist untwist is inertia-
  // sensitive. The O(zeta_dot^2) centrifugal term D2u(zd,zd) is likewise a
  // nonlinear force in update_solution (see "Convective inertia").
  {
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
      const double w = inertia_weight(pt, parameters);
      Madd(e_beta(), i_vbeta(), cL3 * w * R2);   // beta-beta
      Madd(e_eps(), i_veps(), cL3 * w);          // eps-eps
      for (int i = 0; i < 2; i++) {
        const bool free_row = !pin_rigid_ || A[i] != 0;  // pinned phi_0/eta_0 rows
        Madd(e_beta(), i_vphi(A[i]), cL2 * w * R2 * N[i]);  // beta-phi
        Madd(e_eps(), i_veta(A[i]), cL2 * w * N[i]);        // eps-eta
        if (free_row) {
          Madd(e_phi(A[i]), i_vbeta(), cL2 * w * R2 * N[i]);
          Madd(e_eta(A[i]), i_veps(), cL2 * w * N[i]);
        }
        for (int j = 0; j < 2; j++) {
          const double mrr = cL1 * w * N[i] * N[j];       // rho-rho / eta-eta
          Madd(e_rho(A[i]), i_vrho(A[j]), mrr);
          if (free_row) {
            Madd(e_eta(A[i]), i_veta(A[j]), mrr);
            Madd(e_phi(A[i]), i_vphi(A[j]), cL1 * w * R2 * N[i] * N[j]);
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
  // Evaluate the ECG-derived activation nu(t) and set the BCS inputs act = |nu|
  // and act_plus = |nu|_+ (used in the nonlinear residual C assembled by
  // update_solution). nu is the SIGNED contraction rate from the activation
  // function (type `nagumo`); it is deliberately not clamped to [0, 1].
  if (activation_function_ == nullptr) {
    throw std::runtime_error(
        "ChamberCylinder block '" + get_name() +
        "': no activation_function set. Provide an 'activation_function' "
        "(type 'nagumo') on the vessel configuration.");
  }
  const double nu = activation_function_->compute(model->time);
  act = std::abs(nu);
  act_plus = std::max(nu, 0.0);
  // Optional transmural activation delay (Chabiniok 2012 Sec 3.1.4: the wave
  // runs endo -> epi; genet23 uses a uniform nu). act_delay > 0 delays the
  // epicardium by act_delay relative to the endocardium, < 0 the reverse; the
  // same shift applies to deactivation.
  const double delay = parameters[global_param_ids[ParamId::act_delay]];
  const double ri = parameters[global_param_ids[ParamId::Ri]];
  const double re = parameters[global_param_ids[ParamId::Re]];
  act_q.assign(n_quad, act);
  act_plus_q.assign(n_quad, act_plus);
  if (delay != 0.0) {
    for (int q = 0; q < n_quad; q++) {
      const double frac = (quad[q].R - ri) / (re - ri);
      const double dq = delay > 0.0 ? delay * frac : -delay * (1.0 - frac);
      const double nuq = activation_function_->compute(model->time - dq);
      act_q[q] = std::abs(nuq);
      act_plus_q[q] = std::max(nuq, 0.0);
    }
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
  p.gamma = parameters[global_param_ids[ParamId::gamma]];
  p.stretch_pow = parameters[global_param_ids[ParamId::active_stretch_pow]];
  p.visc_mode = (int)std::lround(parameters[global_param_ids[ParamId::visc_mode]]);
  {
    auto P = [&](ParamId id) { return parameters[global_param_ids[id]]; };
    const double g_ff = P(ParamId::gam_ff), g_rr = P(ParamId::gam_rr), g_nn = P(ParamId::gam_nn);
    const double g_fr = P(ParamId::gam_fr), g_fn = P(ParamId::gam_fn), g_rn = P(ParamId::gam_rn);
    p.gam << g_ff, g_fr, g_fn, g_fr, g_rr, g_rn, g_fn, g_rn, g_nn;  // basis (f, r, n)
    p.passive_law = (int)std::lround(P(ParamId::passive_law));
    p.fung_C = P(ParamId::fung_C);
    p.fung_blin = P(ParamId::fung_blin);
    const double b_ff = P(ParamId::fung_bff), b_rr = P(ParamId::fung_brr), b_nn = P(ParamId::fung_bnn);
    const double b_fr = P(ParamId::fung_bfr), b_fn = P(ParamId::fung_bfn), b_rn = P(ParamId::fung_brn);
    p.fung_b << b_ff, b_fr, b_fn, b_fr, b_rr, b_rn, b_fn, b_rn, b_nn;
    p.ho_a = P(ParamId::ho_a);
    p.ho_b = P(ParamId::ho_b);
    p.ho_af = P(ParamId::ho_af);
    p.ho_bf = P(ParamId::ho_bf);
    p.ho_as = P(ParamId::ho_as);
    p.ho_bs = P(ParamId::ho_bs);
    p.ho_afs = P(ParamId::ho_afs);
    p.ho_bfs = P(ParamId::ho_bfs);
    p.ho_sheet = (int)std::lround(P(ParamId::ho_sheet));
  }
  const double sigma_max = parameters[global_param_ids[ParamId::sigma_max]];
  const double Rip = parameters[global_param_ids[ParamId::Ri]];
  const double Lp = parameters[global_param_ids[ParamId::length]];

  BcsParams bp;
  bp.k_s = parameters[global_param_ids[ParamId::k_s]];
  bp.k_0 = parameters[global_param_ids[ParamId::k_0]];
  bp.mu = parameters[global_param_ids[ParamId::mu]];
  bp.alpha = parameters[global_param_ids[ParamId::bcs_alpha]];
  bp.sigma0 = sigma_max;
  bp.fv_reg = parameters[global_param_ids[ParamId::fv_reg]];
  const double tauc_alpha = parameters[global_param_ids[ParamId::tauc_alpha]];
  bp.alpha_tau = tauc_alpha < 0.0 ? bp.alpha : tauc_alpha;
  bp.series = (int)std::lround(parameters[global_param_ids[ParamId::series_law]]);
  bp.starling_out = parameters[global_param_ids[ParamId::starling_out]] > 0.5;

  // Gather block DOFs and their rates.
  auto Y = [&](int lv) { return y[global_var_ids[lv]]; };
  auto Yd = [&](int lv) { return dy[global_var_ids[lv]]; };
  const double beta = Y(i_beta()), eps = Y(i_eps());
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

    // Active-stress input (Genet BCS): this quadrature point's contractile-
    // element bond stress tau_c and strain rate e_c_dot; the fiber stress is
    // computed from the fiber stretch internally.
    ActiveInput ai;
    ai.tauc = Y(i_tauc(q));
    ai.ecdot = Yd(i_ec(q));
    ai.mu = bp.mu;
    ai.series = bp.series;
    ai.ec = Y(i_ec(q));
    // Sainte-Marie placement of the Frank-Starling factor: d(e_c) scales tau_c
    // in the tension; its slope enters the e_c column below.
    double dn0_dec = 0.0;
    if (bp.starling_out) {
      const cplx n0c = frank_starling<cplx>(cplx(ai.ec, 1e-30), n0_ec_, n0_val_);
      ai.tscale = n0c.real();
      dn0_dec = n0c.imag() / 1e-30;
    }

    // Mixed u/p: this element's hydrostatic pressure enters the bulk stress.
    p.p_mix = Y(i_pmix(pt.elem));

    // Base (double) forces plus the analytic active/viscous rate derivatives.
    Kinematics<double> kin = compute_kinematics<double>(xi, pt.R);
    Matrix3d M = fiber_tensor<double>(pt.alpha);
    Matrix3d S = compute_stress<double>(kin, xidot, ai, pt.alpha, p);
    const double I4q = (kin.C * M).trace();  // fiber invariant e_F . C . e_F
    double g[6], dgdtau[6], dgdp[6], Kd[6][6], K[6][6];
    // Mixed u/p: d(g[k])/d(p_mix) = (J C^{-1}) : dE[k], since Sigma^b = p_mix J C^{-1}.
    const Matrix3d JCinv = (kin.J * kin.C.inverse()).eval();
    for (int i = 0; i < 6; i++) {
      g[i] = contract<double>(S, kin.dE[i]);
      dgdtau[i] = contract<double>(M, kin.dE[i]);
      dgdp[i] = contract<double>(JCinv, kin.dE[i]);
      for (int j = 0; j < 6; j++)
        Kd[i][j] = p.visc_mode == 1
                       ? p.gamma * kin.dE[i].trace() * kin.dE[j].trace()
                       : p.gamma * contract<double>(kin.dE[i], kin.dE[j]);
    }
    if (p.visc_mode == 2) {  // Kd_ij = sum_ab gam_ab dE'_i,ab dE'_j,ab (fiber basis)
      const Matrix3d Qb = fiber_basis(pt.alpha);
      Matrix3d dEp[6];
      for (int i = 0; i < 6; i++) dEp[i] = Qb.transpose() * kin.dE[i] * Qb;
      for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
          Kd[i][j] = (p.gam.array() * dEp[i].array() * dEp[j].array()).sum();
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

    // Active-DOF coupling of the mechanical forces: BCS couples to this point's
    // tau_c via d sigma_1D/d tau_c = 1/(1 + e_fib)^pow, 1 + e_fib = sqrt(I4).
    const int active_col = i_tauc(q);
    // sigma_1D = fac * T_fib, T_fib = tscale tau_c + mu e_c_dot, with fac =
    // 1/(1+e_fib)^pow (additive series law) or (1+2e_c)/I4 (multiplicative).
    const double fac = bp.series == 1 ? (1.0 + 2.0 * ai.ec) / I4q
                                      : 1.0 / std::pow(I4q, 0.5 * p.stretch_pow);
    const double active_coef = ai.tscale * fac;  // d sigma_1D / d tau_c
    const double Tfib_q = ai.tscale * ai.tauc + ai.mu * ai.ecdot;
    const double dsig_dec = (bp.series == 1 ? 2.0 * Tfib_q / I4q : 0.0) +
                            dn0_dec * ai.tauc * fac;  // d sigma_1D / d e_c

    // Row (equation) and column (variable) contributions of each local force.
    // xi index: 0=rho(N), 1=rho'(dN), 2=phi'(dN), 3=beta, 4=eps, 5=eta'(dN).
    int req[6][2], cvr[6][2], rn[6], cn[6];
    double rw[6][2], cw[6][2];
    req[0][0] = e_rho(A0); rw[0][0] = N0; req[0][1] = e_rho(A1); rw[0][1] = N1; rn[0] = 2;
    cvr[0][0] = i_rho(A0); cw[0][0] = N0; cvr[0][1] = i_rho(A1); cw[0][1] = N1; cn[0] = 2;
    req[1][0] = e_rho(A0); rw[1][0] = dN0; req[1][1] = e_rho(A1); rw[1][1] = dN1; rn[1] = 2;
    cvr[1][0] = i_rho(A0); cw[1][0] = dN0; cvr[1][1] = i_rho(A1); cw[1][1] = dN1; cn[1] = 2;
    {  // phi': rows skip the inner node when it is pinned (static setting)
      int c = 0;
      if (!pin_rigid_ || A0 != 0) { req[2][c] = e_phi(A0); rw[2][c] = dN0; c++; }
      req[2][c] = e_phi(A1); rw[2][c] = dN1; c++; rn[2] = c;
    }
    cvr[2][0] = i_phi(A0); cw[2][0] = dN0; cvr[2][1] = i_phi(A1); cw[2][1] = dN1; cn[2] = 2;
    req[3][0] = e_beta(); rw[3][0] = 1.0; rn[3] = 1;
    cvr[3][0] = i_beta(); cw[3][0] = 1.0; cn[3] = 1;
    req[4][0] = e_eps(); rw[4][0] = 1.0; rn[4] = 1;
    cvr[4][0] = i_eps(); cw[4][0] = 1.0; cn[4] = 1;
    {  // eta': rows skip the inner node when it is pinned (static setting)
      int c = 0;
      if (!pin_rigid_ || A0 != 0) { req[5][c] = e_eta(A0); rw[5][c] = dN0; c++; }
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
        Kat(re, i_pmix(pt.elem)) += rwt * dgdp[k];
        // BCS: sigma_1D also depends on e_c_dot via the mu*e_c_dot term.
        Kdat(re, i_ec(q)) += rwt * (bp.mu * fac) * dgdtau[k];
        if (dsig_dec != 0.0) Kat(re, i_ec(q)) += rwt * dsig_dec * dgdtau[k];
        for (int m = 0; m < 6; m++) {
          for (int ci = 0; ci < cn[m]; ci++) {
            const int cv = cvr[m][ci];
            Kat(re, cv) += rwt * K[k][m] * cw[m][ci];
            Kdat(re, cv) += rwt * Kd[k][m] * cw[m][ci];
          }
        }
      }
    }

    // --- Convective (centrifugal) inertia D2u(zeta_dot,zeta_dot) (genet23 Eqs 8,10,18) ---
    // The O(zeta_dot^2) part of the acceleration a = D2u(zd,zd) + Du(zdd), omitted
    // from the linear consistent mass. With the angular velocity w_om = beta_dot*Z +
    // phi_dot, contracting D2u (Eq 10) with the virtual velocity Du (Eq 9) and
    // integrating over Theta (2 pi) and Z (moments L, L^2/2, L^3/3) gives the forces
    // conjugate to the rho / phi / beta virtual fields below. Rates are the velocity
    // companion DOFs (w = zeta_dot). With a_c = 2 rho_dot w_om e_theta - (R+rho)
    // w_om^2 e_r and the virtual velocity rho^ e_r + (beta^ Z + phi^)(R+rho) e_theta,
    // every term carries one factor of the CURRENT radius R+rho.
    {
      const double rho0 = parameters[global_param_ids[ParamId::density]];
      const double Wr = N0 * Y(i_vrho(A0)) + N1 * Y(i_vrho(A1));  // rho_dot
      const double Wp = N0 * Y(i_vphi(A0)) + N1 * Y(i_vphi(A1));  // phi_dot
      const double Wb = Y(i_vbeta());                             // beta_dot
      const double c0 = 2.0 * M_PI * rho0 * inertia_weight(pt, parameters);
      const double base = c0 * (pt.R + rho);  // d(base)/d(rho) = c0
      const double L1 = Lp, L2 = Lp * Lp / 2.0, L3 = Lp * Lp * Lp / 3.0;
      const double Iw2 = Wb * Wb * L3 + 2.0 * Wb * Wp * L2 + Wp * Wp * L1;  // int w_om^2
      const double Qp = Wb * L2 + Wp * L1;   // int w_om
      const double Qb = Wb * L3 + Wp * L2;   // int w_om Z
      const int An[2] = {A0, A1};
      const double Nn[2] = {N0, N1};
      // beta (global twist): + 2 rho_dot (int w_om Z)
      Cloc[e_beta()] += base * 2.0 * Wr * Qb;
      Kat(e_beta(), i_vbeta()) += base * 2.0 * Wr * L3;
      for (int j = 0; j < 2; j++) {
        Kat(e_beta(), i_vrho(An[j])) += base * 2.0 * Nn[j] * Qb;
        Kat(e_beta(), i_vphi(An[j])) += base * 2.0 * Wr * L2 * Nn[j];
        Kat(e_beta(), i_rho(An[j])) += c0 * Nn[j] * 2.0 * Wr * Qb;
      }
      for (int i = 0; i < 2; i++) {
        const int a = An[i];
        const double Na = Nn[i];
        // rho: - (int w_om^2)
        Cloc[e_rho(a)] += -base * Iw2 * Na;
        Kat(e_rho(a), i_vbeta()) += -base * Na * (2.0 * Wb * L3 + 2.0 * Wp * L2);
        for (int j = 0; j < 2; j++) {
          Kat(e_rho(a), i_vphi(An[j])) +=
              -base * Na * (2.0 * Wb * L2 + 2.0 * Wp * L1) * Nn[j];
          Kat(e_rho(a), i_rho(An[j])) += -c0 * Nn[j] * Iw2 * Na;
        }
        // phi (skip the inner node when pinned): + 2 rho_dot (int w_om)
        if (!pin_rigid_ || a != 0) {
          Cloc[e_phi(a)] += base * 2.0 * Wr * Qp * Na;
          Kat(e_phi(a), i_vbeta()) += base * 2.0 * Wr * L2 * Na;
          for (int j = 0; j < 2; j++) {
            Kat(e_phi(a), i_vrho(An[j])) += base * 2.0 * Nn[j] * Qp * Na;
            Kat(e_phi(a), i_vphi(An[j])) += base * 2.0 * Wr * L1 * Nn[j] * Na;
            Kat(e_phi(a), i_rho(An[j])) += c0 * Nn[j] * 2.0 * Wr * Qp * Na;
          }
        }
      }
    }

    // --- Current-configuration inertia correction (genet23 Eq A5, Eq 52) ---
    // The consistent mass M = rho0 int (Du[zeta])^T Du[zeta] dOmega uses the
    // CURRENT kinematics gradient Du[zeta], whose twist/azimuthal moment arm is
    // the current radius (R+rho) (Eq A5: entries (R+rho) sin/cos(beta Z+phi)),
    // NOT the reference R. Since Du^T Du is rotation-invariant (sin^2+cos^2=1),
    // only the beta/phi moment arms change: R^2 -> (R+rho)^2. The reference-R^2
    // mass is assembled once in E (update_constant); here we add the state-
    // dependent (R+rho)^2 - R^2 correction to the RESIDUAL (it vanishes at rho=0,
    // a pure addition to the reference dynamics). Rates are the velocity
    // companions (w_dot = zeta_ddot), so the mass enters dC_dydot (Kdat).
    //
    // Jacobian: we deliberately OMIT the d(mass)/d(rho) coupling from dC_dy. That
    // term is a small, non-symmetric rho <-> momentum block (the mass depends on
    // the radial position); keeping it tips the stiff cases (e.g.
    // chamber_cylinder_mixed) past the Newton iteration limit. Dropping it is an
    // inexact-Newton choice: the RESIDUAL is still the exact current-config mass,
    // so the converged solution is unchanged (only the local convergence rate is,
    // and it stays robustly convergent). The total beta/phi mass is (R+rho)^2 > 0.
    {
      const double rho0m = parameters[global_param_ids[ParamId::density]];
      const double cL1 = rho0m * 2.0 * M_PI * Lp;                 // int_Z 1   = L
      const double cL2 = rho0m * 2.0 * M_PI * Lp * Lp / 2.0;      // int_Z Z   = L^2/2
      const double cL3 = rho0m * 2.0 * M_PI * Lp * Lp * Lp / 3.0; // int_Z Z^2 = L^3/3
      const double rcur = pt.R + rho;                 // current radius R + rho
      const double dr2 = rcur * rcur - pt.R * pt.R;   // (R+rho)^2 - R^2
      const double wI = inertia_weight(pt, parameters);
      const int Am[2] = {A0, A1};
      const double Nm[2] = {N0, N1};
      const double aBeta = Yd(i_vbeta());                            // beta_ddot
      const double aPhi = N0 * Yd(i_vphi(A0)) + N1 * Yd(i_vphi(A1)); // phi-field ddot
      // beta momentum row: couples to beta (cL3) and the phi field (cL2).
      const double fBeta = cL3 * aBeta + cL2 * aPhi;
      Cloc[e_beta()] += wI * dr2 * fBeta;
      Kdat(e_beta(), i_vbeta()) += wI * dr2 * cL3;
      for (int j = 0; j < 2; j++)
        Kdat(e_beta(), i_vphi(Am[j])) += wI * dr2 * cL2 * Nm[j];
      // phi momentum rows (skip the inner node when pinned): beta (cL2), phi (cL1).
      const double fPhi = cL2 * aBeta + cL1 * aPhi;
      for (int i = 0; i < 2; i++) {
        const int a = Am[i];
        if (pin_rigid_ && a == 0) continue;
        Cloc[e_phi(a)] += wI * dr2 * Nm[i] * fPhi;
        Kdat(e_phi(a), i_vbeta()) += wI * dr2 * Nm[i] * cL2;
        for (int j = 0; j < 2; j++)
          Kdat(e_phi(a), i_vphi(Am[j])) += wI * dr2 * Nm[i] * cL1 * Nm[j];
      }
    }

    // --- DIAGNOSTIC (genet23 App. A7/A9 as printed): Z-independent velocity memory ---
    // A7 writes the inertia as 4 pi int rho0 zeta^T (M dzeta/dt^2 - N^T v^n/dt) dR
    // with N = int_0^L Du dZ, i.e. the stored velocity v^n(R) is pulled out of the
    // Z-integral as if it did not depend on Z. The true velocity Du zeta_dot does
    // (twist and axial stretch grow linearly in Z). Each midpoint step then keeps
    // only the Z-mean velocity and discards the momentum
    //   dM zeta_dot,  dM = int Du^T Du dZ - (1/L) N^T N   (>= 0),
    // which acts as a dashpot 2 dM / dt (dt = the paper's step, vmem_dt). At small
    // twist only dM_beta,beta = dM_eps,eps = (R+rho)^2 L^3/12 resp. L^3/12 survive.
    // vmem_dt = 0 (default) disables it. Z-integrals use 8-point Gauss-Legendre;
    // the state dependence of dM is lagged in the tangent (diagnostic only).
    {
      const double vdt = parameters[global_param_ids[ParamId::vmem_dt]];
      if (vdt > 0.0) {
        const double rho0v = parameters[global_param_ids[ParamId::density]];
        const double cvm = (2.0 / vdt) * 2.0 * M_PI * rho0v * inertia_weight(pt, parameters);
        const double phi_q = N0 * Y(i_phi(A0)) + N1 * Y(i_phi(A1));
        const double r = pt.R + rho;
        static const double gz[8] = {-0.9602898564975363, -0.7966664774136267,
                                     -0.5255324099163290, -0.1834346424956498,
                                     0.1834346424956498,  0.5255324099163290,
                                     0.7966664774136267,  0.9602898564975363};
        static const double gwz[8] = {0.1012285362903763, 0.2223810344533745,
                                      0.3137066458778873, 0.3626837833783620,
                                      0.3626837833783620, 0.3137066458778873,
                                      0.2223810344533745, 0.1012285362903763};
        // vmem_mask selects which velocity components lose their Z-deviation:
        // bit 0 the in-plane (e_r, e_theta: twist) part, bit 1 the axial part.
        const int vmask = (int)std::lround(parameters[global_param_ids[ParamId::vmem_mask]]);
        double Mi[5][5] = {}, Ni[3][5] = {};  // local order: rho, beta, phi, eps, eta
        for (int g = 0; g < 8; g++) {
          const double Z = 0.5 * Lp * (1.0 + gz[g]), wz = 0.5 * Lp * gwz[g];
          const double ps = beta * Z + phi_q, c = std::cos(ps), s = std::sin(ps);
          const double Du[3][5] = {{c, -s * Z * r, -s * r, 0.0, 0.0},
                                   {s, c * Z * r, c * r, 0.0, 0.0},
                                   {0.0, 0.0, 0.0, Z, 1.0}};
          for (int k = 0; k < 3; k++) {
            if (!(vmask & (k < 2 ? 1 : 2))) continue;
            for (int m = 0; m < 5; m++) {
              Ni[k][m] += wz * Du[k][m];
              for (int n = 0; n < 5; n++) Mi[m][n] += wz * Du[k][m] * Du[k][n];
            }
          }
        }
        double dM[5][5];
        for (int m = 0; m < 5; m++)
          for (int n = 0; n < 5; n++) {
            double nn = 0.0;
            for (int k = 0; k < 3; k++) nn += Ni[k][m] * Ni[k][n];
            dM[m][n] = Mi[m][n] - nn / Lp;
          }
        // Local generalized DOFs -> global rows/velocity columns (value fields use N).
        auto rows = [&](int m, int* e, double* wt) -> int {
          if (m == 1) { e[0] = e_beta(); wt[0] = 1.0; return 1; }
          if (m == 3) { e[0] = e_eps(); wt[0] = 1.0; return 1; }
          int n = 0;
          const int An[2] = {A0, A1};
          const double Nn[2] = {N0, N1};
          for (int i = 0; i < 2; i++) {
            if (pin_rigid_ && An[i] == 0 && m != 0) continue;  // pinned phi_0 / eta_0
            e[n] = m == 0 ? e_rho(An[i]) : m == 2 ? e_phi(An[i]) : e_eta(An[i]);
            wt[n++] = Nn[i];
          }
          return n;
        };
        auto cols = [&](int m, int* v, double* wt) -> int {
          if (m == 1) { v[0] = i_vbeta(); wt[0] = 1.0; return 1; }
          if (m == 3) { v[0] = i_veps(); wt[0] = 1.0; return 1; }
          v[0] = m == 0 ? i_vrho(A0) : m == 2 ? i_vphi(A0) : i_veta(A0);
          v[1] = m == 0 ? i_vrho(A1) : m == 2 ? i_vphi(A1) : i_veta(A1);
          wt[0] = N0; wt[1] = N1;
          return 2;
        };
        double zd[5];
        for (int n = 0; n < 5; n++) {
          int v[2]; double wt[2]; const int nc = cols(n, v, wt);
          zd[n] = 0.0;
          for (int j = 0; j < nc; j++) zd[n] += wt[j] * Y(v[j]);
        }
        for (int m = 0; m < 5; m++) {
          int e[2]; double we[2]; const int nr = rows(m, e, we);
          double fm = 0.0;
          for (int n = 0; n < 5; n++) fm += dM[m][n] * zd[n];
          for (int i = 0; i < nr; i++) {
            Cloc[e[i]] += cvm * we[i] * fm;
            for (int n = 0; n < 5; n++) {
              if (dM[m][n] == 0.0) continue;
              int v[2]; double wv[2]; const int nc = cols(n, v, wv);
              for (int j = 0; j < nc; j++) Kat(e[i], v[j]) += cvm * we[i] * dM[m][n] * wv[j];
            }
          }
        }
      }
    }

    // --- Optional end-cap support (Chabiniok 2012 Eq 13, Sainte-Marie 2006 base) ---
    // Tangential visco-elastic traction t_theta = -(k u_theta + c v_theta) on the
    // wall's end annuli Z_e in {0 (base), L (apex)}, with u_theta = r theta,
    // theta = beta Z_e + phi(R), r = R + rho. Its virtual power is
    // -2 pi int r^2 (k theta + c theta_dot)(beta^ Z_e + phi^) R dR per end.
    // genet23 has no such support (Sec 2.4); support_ends = 0 (default) is off.
    {
      const int ends = (int)std::lround(parameters[global_param_ids[ParamId::support_ends]]);
      const double ks = parameters[global_param_ids[ParamId::support_k]];
      const double cs = parameters[global_param_ids[ParamId::support_c]];
      if (ends != 0 && (ks != 0.0 || cs != 0.0)) {
        const int An[2] = {A0, A1};
        const double Nn[2] = {N0, N1};
        const double r = pt.R + rho;
        const double phq = N0 * Y(i_phi(A0)) + N1 * Y(i_phi(A1));
        const double phq_d = N0 * Yd(i_phi(A0)) + N1 * Yd(i_phi(A1));
        const double base = 2.0 * M_PI * pt.w * r * r;
        for (int e = 0; e < 2; e++) {
          if (!(ends & (1 << e))) continue;  // bit 0 = base (Z=0), bit 1 = apex (Z=L)
          const double Ze = e == 0 ? 0.0 : Lp;
          const double th = beta * Ze + phq, thd = beta_d * Ze + phq_d;
          const double f = ks * th + cs * thd;
          // rows: beta (weight Ze) and phi nodes (weight N_a)
          int rows[3]; double rw_[3]; int nr = 0;
          if (Ze != 0.0) { rows[nr] = e_beta(); rw_[nr++] = Ze; }
          for (int i = 0; i < 2; i++) {
            if (pin_rigid_ && An[i] == 0) continue;
            rows[nr] = e_phi(An[i]); rw_[nr++] = Nn[i];
          }
          for (int ir = 0; ir < nr; ir++) {
            const double w_r = base * rw_[ir];
            Cloc[rows[ir]] += w_r * f;
            Kat(rows[ir], i_beta()) += w_r * ks * Ze;
            Kdat(rows[ir], i_beta()) += w_r * cs * Ze;
            for (int j = 0; j < 2; j++) {
              Kat(rows[ir], i_phi(An[j])) += w_r * ks * Nn[j];
              Kdat(rows[ir], i_phi(An[j])) += w_r * cs * Nn[j];
              Kat(rows[ir], i_rho(An[j])) += 2.0 * M_PI * pt.w * 2.0 * r * Nn[j] * rw_[ir] * f;
            }
          }
        }
      }
    }

    // --- Mixed u/p: weak incompressibility constraint for this element (P0) ---
    // Residual C[e_pmix] = sum_{q in e} pref (J-1) = 0. J depends only on
    // rho, rho', eps; map the analytic dJ/dxi through the same columns as the
    // force assembly. The pressure-pressure block is zero (saddle-point form).
    {
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
    {
      const double sfib = std::sqrt(I4q);
      const double e_fib = sfib - 1.0;
      const double e_1D = 0.5 * (I4q - 1.0);  // Green-Lagrange fiber strain
      const double e_kin = bp.series == 1 ? e_1D : e_fib;
      const int erow[3] = {e_ec(q), e_tauc(q), e_kc(q)};
      const int avar[3] = {i_ec(q), i_tauc(q), i_kc(q)};
      double a[3] = {Y(i_ec(q)), Y(i_tauc(q)), Y(i_kc(q))};
      double ad[3] = {Yd(i_ec(q)), Yd(i_tauc(q)), Yd(i_kc(q))};
      double Ra[3];
      bcs_residual<double>(a, ad, e_kin, act_q[q], act_plus_q[q], bp, n0_ec_, n0_val_, Ra);
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
          bcs_residual<cplx>(ac, adc, cplx(e_kin), act_q[q], act_plus_q[q], bp, n0_ec_,
                             n0_val_, Rc);
          for (int ii = 0; ii < 3; ii++)
            Kat(erow[ii], avar[jj]) += std::imag(Rc[ii]) / h;
          ac[jj] = a[jj];
          adc[jj] += cplx(0.0, h);
          bcs_residual<cplx>(ac, adc, cplx(e_kin), act_q[q], act_plus_q[q], bp, n0_ec_,
                             n0_val_, Rc);
          for (int ii = 0; ii < 3; ii++)
            Kdat(erow[ii], avar[jj]) += std::imag(Rc[ii]) / h;
          adc[jj] = ad[jj];
        }
      }
      // Only the e_c ODE couples to the wall kinematics, through the fiber
      // strain: dR_ec/de_fib = -k_s/mu, de_fib/dxi_m = dgdtau[m]/sqrt(I4). The
      // tau_c/k_c ODEs couple to the mechanics only via e_c (handled above).
      const double ecq = Y(i_ec(q));
      const double gl3 = (1.0 + 2.0 * ecq) * (1.0 + 2.0 * ecq) * (1.0 + 2.0 * ecq);
      // d res_ec / d e_kin, times d e_kin / d xi_m = dgdtau[m] / sqrt(I4) (e_fib)
      // or dgdtau[m] (e_1D = e_F.E.e_F).
      const double dRec_defib =
          bp.series == 1 ? -bp.k_s * ((1.0 + 2.0 * e_1D) + 2.0 * (e_1D - ecq)) / gl3 / bp.mu
                         : -bp.k_s / bp.mu;
      const double dekin_scale = bp.series == 1 ? 1.0 : 1.0 / sfib;
      for (int m = 0; m < 6; m++) {
        // (default path keeps the original operation order bit-for-bit)
        const double coef = bp.series == 1 ? dRec_defib * dgdtau[m] * dekin_scale
                                           : dRec_defib * dgdtau[m] / sfib;
        for (int ci = 0; ci < cn[m]; ci++)
          Kat(e_ec(q), cvr[m][ci]) += coef * cw[m][ci];
      }
    }
  }

  // --- Optional epicardial spring-dashpot support (off in genet23) ---
  // Traction t = -(k u + c u_dot) per component on the lateral epicardium
  // R = Re, Z in [0, L] (reference area 2 pi Re dZ): normal u_r = rho_e
  // (Pfaller 2019 Eq 5), circumferential u_th = r_e theta with theta = beta Z +
  // phi_e and r_e = Re + rho_e, axial u_z = eps Z + eta_e (Moireau 2011 Eq 9).
  {
    auto P = [&](ParamId id) { return parameters[global_param_ids[id]]; };
    const double kn = P(ParamId::epi_kn), cn = P(ParamId::epi_cn);
    const double kt = P(ParamId::epi_kt), ct = P(ParamId::epi_ct);
    const double kz = P(ParamId::epi_kz), cz = P(ParamId::epi_cz);
    const int ae = n_node - 1;
    const double A = 2.0 * M_PI * node_R[ae];
    const double L1 = Lp, L2 = Lp * Lp / 2.0, L3 = Lp * Lp * Lp / 3.0;
    if (kn != 0.0 || cn != 0.0) {
      const int r = e_rho(ae), c = i_rho(ae);
      Cloc[r] += A * L1 * (kn * Y(c) + cn * Yd(c));
      Kat(r, c) += A * L1 * kn;
      Kdat(r, c) += A * L1 * cn;
    }
    // Z-moments of a field u = x Z + w against the virtual fields Z and 1.
    auto lateral = [&](int ex, int ew, int ix, int iw, double k, double c, double s,
                       double ds_drho) {
      const double x = Y(ix), w = Y(iw), xd = Yd(ix), wd = Yd(iw);
      const double fx = k * (x * L3 + w * L2) + c * (xd * L3 + wd * L2);
      const double fw = k * (x * L2 + w * L1) + c * (xd * L2 + wd * L1);
      Cloc[ex] += A * s * fx;
      Cloc[ew] += A * s * fw;
      Kat(ex, ix) += A * s * k * L3; Kat(ex, iw) += A * s * k * L2;
      Kat(ew, ix) += A * s * k * L2; Kat(ew, iw) += A * s * k * L1;
      Kdat(ex, ix) += A * s * c * L3; Kdat(ex, iw) += A * s * c * L2;
      Kdat(ew, ix) += A * s * c * L2; Kdat(ew, iw) += A * s * c * L1;
      if (ds_drho != 0.0) {
        Kat(ex, i_rho(ae)) += A * ds_drho * fx;
        Kat(ew, i_rho(ae)) += A * ds_drho * fw;
      }
    };
    if (kt != 0.0 || ct != 0.0) {
      const double r_e = node_R[ae] + Y(i_rho(ae));
      lateral(e_beta(), e_phi(ae), i_beta(), i_phi(ae), kt, ct, r_e * r_e, 2.0 * r_e);
    }
    if (kz != 0.0 || cz != 0.0)
      lateral(e_eps(), e_eta(ae), i_eps(), i_eta(ae), kz, cz, 1.0, 0.0);
  }

  // --- Probe (off in genet23): generalized twist dashpot D beta_dot, power D beta_dot^2 ---
  {
    const double Dt = parameters[global_param_ids[ParamId::twist_c]];
    if (Dt != 0.0) {
      Cloc[e_beta()] += Dt * beta_d;
      Kdat(e_beta(), i_beta()) += Dt;
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

  // --- Scatter local contributions into the sparse system (assign) ---
  for (int i = 0; i < n_eqn; i++) {
    system.C(global_eqn_ids[i]) = Cloc[i];
    for (int j = 0; j < n_var; j++) {
      system.dC_dy.coeffRef(global_eqn_ids[i], global_var_ids[j]) = Kat(i, j);
      system.dC_dydot.coeffRef(global_eqn_ids[i], global_var_ids[j]) = Kdat(i, j);
    }
  }
}

void ChamberCylinder::set_frank_starling_curve(std::vector<double> ec,
                                               std::vector<double> val) {
  if (ec.size() != val.size() || ec.size() < 2) {
    throw std::runtime_error(
        "ChamberCylinder::set_frank_starling_curve: e_c and n0 arrays must have "
        "equal length >= 2 (got " +
        std::to_string(ec.size()) + " and " + std::to_string(val.size()) + ").");
  }
  n0_ec_ = std::move(ec);
  n0_val_ = std::move(val);
}
