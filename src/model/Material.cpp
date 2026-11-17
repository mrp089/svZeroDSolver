// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "Material.h"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace {

using cplx = std::complex<double>;
template <typename T>
using Mat3 = Eigen::Matrix<T, 3, 3>;
using Vec3 = Eigen::Matrix<double, 3, 1>;

/// Real part of a scalar, for comparisons on complex-step arguments.
inline double re(double x) { return x; }
inline double re(const cplx& x) { return x.real(); }

/// Upper clamp compatible with complex-step: x below m, the constant m above.
template <typename T>
inline T clamp_max(T x, double m) {
  return re(x) < m ? x : T(m);
}

/// Stress of the energy
/// W = W1 (I1b - 3) + W2 (I2b - 3) + C0 exp(C1 (I1b - 3)^2)
///   + C2 exp(C3 (I4b - 1)^2)
/// with the coefficients c = {W1, W2, C0, C1, C2, C3}.
template <typename T>
Mat3<T> stress_t(const Mat3<T>& C, const T& J, const Vec3& f,
                 const double c[6]) {
  const Eigen::Matrix<T, 3, 1> ef = f.cast<T>();
  const Mat3<T> M = ef * ef.transpose();
  const Mat3<T> I = Mat3<T>::Identity();
  const Mat3<T> Cinv = C.inverse();
  const T I1 = C.trace();
  const T I2 = 0.5 * (I1 * I1 - (C.array() * C.array()).sum());
  const T I4 = (C * M).trace();
  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const T Jm43 = Jm23 * Jm23;
  const T I1b = Jm23 * I1;
  const T I4b = Jm23 * I4;

  // Derivatives of the isochoric invariants with respect to C
  const Mat3<T> dI1b = Jm23 * (I - (I1 / 3.0) * Cinv);
  const Mat3<T> dI2b = Jm43 * (I1 * I - C - (2.0 / 3.0) * I2 * Cinv);
  const Mat3<T> dI4b = Jm23 * (M - (I4 / 3.0) * Cinv);

  // dW/d(isochoric invariant); the exponents are capped so that non-physical
  // trial states stay finite.
  const T a1 = c[3] * (I1b - 3.0) * (I1b - 3.0);
  const T a4 = c[5] * (I4b - 1.0) * (I4b - 1.0);
  const T w1 =
      c[0] + 2.0 * c[2] * c[3] * (I1b - 3.0) * std::exp(clamp_max(a1, 300.0));
  const T w2 = T(c[1]);
  const T w4 = 2.0 * c[4] * c[5] * (I4b - 1.0) * std::exp(clamp_max(a4, 300.0));
  return 2.0 * (w1 * dI1b + w2 * dI2b + w4 * dI4b);
}

/// Double contraction A:B (without complex conjugation, for complex-step).
template <typename T>
inline T contract(const Mat3<T>& A, const Mat3<T>& B) {
  return (A.array() * B.array()).sum();
}

/// dW/dI of a Holzapfel-Ogden fiber or sheet energy
/// W = a/(2b) chi (exp(b E^2) - 1) with E = I - 1 and the smoothed Heaviside
/// chi = 1/(1 + exp(-k E)).
template <typename T>
inline T ho_family(T E, double a, double b, double k) {
  if (a == 0.0) return T(0.0);
  const T o = 1.0 / (std::exp(clamp_max(k * E, 300.0)) + 1.0);  // 1 - chi
  const T rexp = std::exp(clamp_max(b * E * E, 300.0));
  return a * ((1.0 - o) * E * rexp + 0.5 * k * (o - o * o) / b * (rexp - 1.0));
}

/// Holzapfel-Ogden constants a, b, a4f, b4f, a4s, b4s, afs, bfs, k.
struct HolzapfelOgdenConstants {
  double a, b, af, bf, as, bs, afs, bfs, k;
};

/// Stress of the svMultiPhysics HolzapfelOgden energy with the isochoric
/// invariants I1b, I4b = J^-2/3 e_F . C . e_F, I4sb = J^-2/3 e_S . C . e_S
/// and I8b = J^-2/3 e_F . C . e_S.
template <typename T>
Mat3<T> holzapfel_ogden_stress(const Mat3<T>& C, const T& J, const Vec3& f,
                               const Vec3& s,
                               const HolzapfelOgdenConstants& h) {
  const Eigen::Matrix<T, 3, 1> ef = f.cast<T>(), es = s.cast<T>();
  const Mat3<T> I = Mat3<T>::Identity();
  const Mat3<T> Cinv = C.inverse();
  const Mat3<T> Hf = ef * ef.transpose();
  const Mat3<T> Hs = es * es.transpose();
  const Mat3<T> Hfs = 0.5 * (ef * es.transpose() + es * ef.transpose());
  const T Jm23 = std::pow(J, -2.0 / 3.0);
  const T I1 = C.trace();
  const T I4 = contract<T>(C, Hf);
  const T I6 = contract<T>(C, Hs);
  const T I8 = contract<T>(C, Hfs);
  const Mat3<T> dI1b = Jm23 * (I - (I1 / 3.0) * Cinv);
  const Mat3<T> dI4b = Jm23 * (Hf - (I4 / 3.0) * Cinv);
  const Mat3<T> dI6b = Jm23 * (Hs - (I6 / 3.0) * Cinv);
  const Mat3<T> dI8b = Jm23 * (Hfs - (I8 / 3.0) * Cinv);
  const T I8b = Jm23 * I8;
  const T wa = 0.5 * h.a * std::exp(clamp_max(h.b * (Jm23 * I1 - 3.0), 300.0));
  const T wf = ho_family<T>(Jm23 * I4 - 1.0, h.af, h.bf, h.k);
  const T ws = ho_family<T>(Jm23 * I6 - 1.0, h.as, h.bs, h.k);
  const T wfs = h.afs * I8b * std::exp(clamp_max(h.bfs * I8b * I8b, 300.0));
  return 2.0 * (wa * dI1b + wf * dI4b + ws * dI6b + wfs * dI8b);
}

}  // namespace

Material::Material(
    const std::vector<std::pair<std::string, InputParameter>>& props)
    : ParameterHolder(props) {}

std::unique_ptr<Material> Material::create(const std::string& type_str) {
  if (type_str == "mooney_rivlin") {
    return std::make_unique<MooneyRivlinMaterial>();
  }
  if (type_str == "exponential") {
    return std::make_unique<ExponentialMaterial>();
  }
  if (type_str == "mooney_rivlin_exponential") {
    return std::make_unique<MooneyRivlinExponentialMaterial>();
  }
  if (type_str == "holzapfel_ogden") {
    return std::make_unique<HolzapfelOgdenMaterial>();
  }
  throw std::runtime_error(
      "Unknown material type '" + type_str +
      "'. Must be one of: mooney_rivlin, exponential, "
      "mooney_rivlin_exponential, holzapfel_ogden");
}

double MooneyRivlinExponentialMaterial::get(const std::string& name) const {
  const auto it = params_.find(name);
  return it == params_.end() ? 0.0 : it->second;
}

Material::Mat3 MooneyRivlinExponentialMaterial::stress(const Mat3& C, double J,
                                                       const Vec3& f,
                                                       const Vec3& s) const {
  (void)s;
  const double c[6] = {get("W1"), get("W2"), get("C0"),
                       get("C1"), get("C2"), get("C3")};
  return stress_t<double>(C, J, f, c);
}

Material::Mat3c MooneyRivlinExponentialMaterial::stress(const Mat3c& C,
                                                        cplx J, const Vec3& f,
                                                        const Vec3& s) const {
  (void)s;
  const double c[6] = {get("W1"), get("W2"), get("C0"),
                       get("C1"), get("C2"), get("C3")};
  return stress_t<cplx>(C, J, f, c);
}

// --- Holzapfel-Ogden -------------------------------------------------------

void HolzapfelOgdenMaterial::finalize() {
  if ((params_.at("a4f") != 0.0 && params_.at("b4f") <= 0.0) ||
      (params_.at("a4s") != 0.0 && params_.at("b4s") <= 0.0)) {
    throw std::runtime_error(
        "holzapfel_ogden material: b4f and b4s must be positive for nonzero "
        "a4f and a4s.");
  }
}

namespace {
HolzapfelOgdenConstants ho_constants(const std::map<std::string, double>& p) {
  return {p.at("a"),   p.at("b"),   p.at("a4f"), p.at("b4f"), p.at("a4s"),
          p.at("b4s"), p.at("afs"), p.at("bfs"), p.at("k")};
}
}  // namespace

Material::Mat3 HolzapfelOgdenMaterial::stress(const Mat3& C, double J,
                                              const Vec3& f,
                                              const Vec3& s) const {
  return holzapfel_ogden_stress<double>(C, J, f, s, ho_constants(params_));
}

Material::Mat3c HolzapfelOgdenMaterial::stress(const Mat3c& C, cplx J,
                                               const Vec3& f,
                                               const Vec3& s) const {
  return holzapfel_ogden_stress<cplx>(C, J, f, s, ho_constants(params_));
}
