// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ChamberSpherePQ.h"

#include <Eigen/Dense>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ChamberSphere parameter order. Must match ChamberSphere::ParamId.
const std::vector<std::string> kParamNames = {
    "rho",       "thick0",    "radius0", "W1",    "W2",    "eta",
    "sigma_max", "alpha_max", "alpha_min", "tsys", "tdias", "steepness"};
enum P { RHO, THICK0, RADIUS0, W1, W2, ETA, SIGMA, AMAX, AMIN, TSYS, TDIAS, STEEP };

using Vec = std::vector<double>;

// Cumulative trapezoidal integral with I[0] = 0.
Vec cumint(const Vec& f, const Vec& t) {
  Vec out(f.size(), 0.0);
  for (size_t i = 1; i < f.size(); i++) {
    out[i] = out[i - 1] + 0.5 * (f[i] + f[i - 1]) * (t[i] - t[i - 1]);
  }
  return out;
}

// Second-order central derivative on a (possibly non-uniform) grid; first-order
// one-sided at the boundaries (matches numpy.gradient).
Vec ddt(const Vec& x, const Vec& t) {
  const size_t n = x.size();
  Vec g(n, 0.0);
  if (n < 2) return g;
  g[0] = (x[1] - x[0]) / (t[1] - t[0]);
  g[n - 1] = (x[n - 1] - x[n - 2]) / (t[n - 1] - t[n - 2]);
  for (size_t i = 1; i + 1 < n; i++) {
    const double hs = t[i] - t[i - 1];
    const double hd = t[i + 1] - t[i];
    g[i] = -hd / (hs * (hs + hd)) * x[i - 1] +
           (hd - hs) / (hs * hd) * x[i] +
           hs / (hd * (hs + hd)) * x[i + 1];
  }
  return g;
}

// Per-time-point activation value and its parameter derivatives, mirroring
// ChamberSphere::get_elastance_values.
struct Activation {
  double act, act_plus;
  double dact[12], dact_plus[12];  // only AMAX..STEEP nonzero
};
Activation activation(double tc, const double* a) {
  const double steep = a[STEEP];
  const double as = (tc - a[TSYS]) / steep;
  const double ad = (tc - a[TDIAS]) / steep;
  const double ths = std::tanh(as), thd = std::tanh(ad);
  const double Sp = 0.5 * (1.0 + ths), Sm = 0.5 * (1.0 - thd);
  const double f = Sp * Sm;
  const double act_t = a[AMAX] * f + a[AMIN] * (1.0 - f);
  const double sech2s = 1.0 - ths * ths, sech2d = 1.0 - thd * thd;
  // df/d(timing params)
  const double df_dtsys = (-0.5 * sech2s / steep) * Sm;
  const double df_dtdias = Sp * (0.5 * sech2d / steep);
  const double df_dsteep =
      (0.5 * sech2s * (-(tc - a[TSYS]) / (steep * steep))) * Sm +
      Sp * (0.5 * sech2d * ((tc - a[TDIAS]) / (steep * steep)));
  const double damin_factor = a[AMAX] - a[AMIN];
  double dact_t[12] = {0};
  dact_t[AMAX] = f;
  dact_t[AMIN] = 1.0 - f;
  dact_t[TSYS] = damin_factor * df_dtsys;
  dact_t[TDIAS] = damin_factor * df_dtdias;
  dact_t[STEEP] = damin_factor * df_dsteep;

  Activation r{};
  r.act = std::abs(act_t);
  r.act_plus = std::fmax(act_t, 0.0);
  const double sgn = (act_t > 0.0) - (act_t < 0.0);
  const double heav = act_t > 0.0 ? 1.0 : 0.0;
  for (int k = 0; k < 12; k++) {
    r.dact[k] = sgn * dact_t[k];
    r.dact_plus[k] = heav * dact_t[k];
  }
  return r;
}

// Reconstructed internal trajectory and the active-stress residual G.
struct Recon {
  Vec R, velo, dvelo, stress, tau, dtau, G;
  std::vector<Activation> act;
};

Recon reconstruct(const double* a, const Vec& t, const Vec& P, const Vec& F,
                  const Vec& Iflow, double period) {
  const size_t n = t.size();
  const double r0 = a[RADIUS0];
  Recon rc;
  rc.R.resize(n);
  rc.velo.resize(n);
  rc.stress.resize(n);
  rc.tau.resize(n);
  rc.act.resize(n);
  for (size_t i = 0; i < n; i++) {
    // Reference-state V0 = 4/3 pi r0^3 ; absolute cavity volume from net flow.
    const double Vcav = (4.0 / 3.0) * M_PI * r0 * r0 * r0 + Iflow[i];
    const double R = std::cbrt(3.0 * Vcav / (4.0 * M_PI));
    const double velo = F[i] / (4.0 * M_PI * R * R);
    rc.R[i] = R;
    rc.velo[i] = velo;
  }
  rc.dvelo = ddt(rc.velo, t);
  for (size_t i = 0; i < n; i++) {
    const double R = rc.R[i], velo = rc.velo[i], dvelo = rc.dvelo[i];
    const double CG = (R / r0) * (R / r0);
    const double dCG = 2.0 * R * velo / (r0 * r0);
    const double stress =
        (P[i] * CG - a[RHO] * a[THICK0] * dvelo) * r0 * r0 / (a[THICK0] * R);
    const double passive = 4.0 * (1.0 - std::pow(CG, -3)) * (a[W1] + CG * a[W2]) +
                           a[ETA] * dCG * (1.0 + 2.0 * std::pow(CG, -6));
    rc.stress[i] = stress;
    rc.tau[i] = stress - passive;
    const double tc = period > 0.0 ? std::fmod(t[i], period) : t[i];
    rc.act[i] = activation(tc, a);
  }
  rc.dtau = ddt(rc.tau, t);
  rc.G.resize(n);
  for (size_t i = 0; i < n; i++) {
    rc.G[i] = rc.dtau[i] + rc.act[i].act * rc.tau[i] - a[SIGMA] * rc.act[i].act_plus;
  }
  return rc;
}

// d(tau)/d(param) array for a single non-activation parameter (RHO..ETA,
// RADIUS0). Returns zeros for SIGMA and activation params.
Vec dtau_dparam(int p, const double* a, const Vec& t, const Vec& P,
                const Vec& F, const Recon& rc) {
  const size_t n = t.size();
  const double r0 = a[RADIUS0];
  Vec d(n, 0.0);
  if (p == SIGMA || p >= AMAX) return d;  // tau independent of these

  // For radius0, tau depends on r0 both explicitly and through R, velo, dvelo.
  Vec dvelo_dr0, ddvelo_dr0;
  if (p == RADIUS0) {
    dvelo_dr0.resize(n);
    for (size_t i = 0; i < n; i++) {
      const double R = rc.R[i];
      dvelo_dr0[i] = -2.0 * rc.velo[i] * r0 * r0 / (R * R * R);
    }
    ddvelo_dr0 = ddt(dvelo_dr0, t);
  }

  for (size_t i = 0; i < n; i++) {
    const double R = rc.R[i], velo = rc.velo[i], dvelo = rc.dvelo[i];
    switch (p) {
      case RHO:
        d[i] = -dvelo * r0 * r0 / R;
        break;
      case THICK0:
        d[i] = -P[i] * R / (a[THICK0] * a[THICK0]);
        break;
      case W1:
        d[i] = -4.0 + 4.0 * std::pow(r0, 6) / std::pow(R, 6);
        break;
      case W2:
        d[i] = 4.0 * (-std::pow(R, 6) + std::pow(r0, 6)) /
               (std::pow(R, 4) * r0 * r0);
        break;
      case ETA:
        d[i] = -2.0 * velo * (std::pow(R, 12) + 2.0 * std::pow(r0, 12)) /
               (std::pow(R, 11) * r0 * r0);
        break;
      case RADIUS0: {
        const double dR_dr0 = r0 * r0 / (R * R);
        const double dtau_dr0_expl =
            8.0 * R * R * a[W2] / std::pow(r0, 3) +
            4.0 * R * a[ETA] * velo / std::pow(r0, 3) -
            2.0 * dvelo * r0 * a[RHO] / R +
            16.0 * a[W2] * std::pow(r0, 3) / std::pow(R, 4) +
            24.0 * a[W1] * std::pow(r0, 5) / std::pow(R, 6) -
            40.0 * a[ETA] * std::pow(r0, 9) * velo / std::pow(R, 11);
        const double dtau_dR =
            P[i] / a[THICK0] - 8.0 * R * a[W2] / (r0 * r0) -
            2.0 * a[ETA] * velo / (r0 * r0) +
            dvelo * r0 * r0 * a[RHO] / (R * R) -
            16.0 * a[W2] * std::pow(r0, 4) / std::pow(R, 5) -
            24.0 * a[W1] * std::pow(r0, 6) / std::pow(R, 7) +
            44.0 * a[ETA] * std::pow(r0, 10) * velo / std::pow(R, 12);
        const double dtau_dvelo =
            -2.0 * a[ETA] * (std::pow(R, 12) + 2.0 * std::pow(r0, 12)) /
            (std::pow(R, 11) * r0 * r0);
        const double dtau_ddvelo = -r0 * r0 * a[RHO] / R;
        d[i] = dtau_dr0_expl + dtau_dR * dR_dr0 +
               dtau_dvelo * dvelo_dr0[i] + dtau_ddvelo * ddvelo_dr0[i];
        break;
      }
    }
  }
  return d;
}

}  // namespace

bool is_chamber_pq_mode(const nlohmann::json& config) {
  // A ChamberSphere is always calibrated from port pressure/flow only.
  if (!config.contains("vessels")) return false;
  for (const auto& v : config["vessels"]) {
    if (v.value("zero_d_element_type", std::string()) == "ChamberSphere")
      return true;
  }
  return false;
}

nlohmann::json calibrate_chamber_pq(const nlohmann::json& config) {
  auto out = nlohmann::json(config);

  // Locate the ChamberSphere vessel.
  const nlohmann::json* vessel = nullptr;
  for (const auto& v : config["vessels"]) {
    if (v.value("zero_d_element_type", std::string()) == "ChamberSphere") {
      vessel = &v;
      break;
    }
  }
  if (vessel == nullptr) {
    throw std::runtime_error(
        "[svzerodcalibrator] P/Q ChamberSphere calibration: no ChamberSphere "
        "vessel found.");
  }
  const std::string name = (*vessel)["vessel_name"].get<std::string>();

  if (!vessel->contains("boundary_conditions") ||
      !(*vessel)["boundary_conditions"].contains("inlet") ||
      !(*vessel)["boundary_conditions"].contains("outlet")) {
    throw std::runtime_error(
        "[svzerodcalibrator] P/Q ChamberSphere calibration currently requires "
        "the ChamberSphere vessel to have inlet and outlet boundary "
        "conditions.");
  }
  const std::string inbc = (*vessel)["boundary_conditions"]["inlet"];
  const std::string outbc = (*vessel)["boundary_conditions"]["outlet"];
  const std::string innode = inbc + ":" + name;
  const std::string outnode = name + ":" + outbc;

  // Read observations: port pressure (P = outlet = inlet pressure) and net flow.
  const auto& y = config["y"];
  auto get = [&](const std::string& key) {
    if (!y.contains(key)) {
      throw std::runtime_error(
          "[svzerodcalibrator] P/Q ChamberSphere calibration: missing "
          "observation '" + key + "'.");
    }
    return y[key].get<Vec>();
  };
  const Vec Pout = get("pressure:" + outnode);
  const Vec Qin = get("flow:" + innode);
  const Vec Qout = get("flow:" + outnode);
  if (!config.contains("time")) {
    throw std::runtime_error(
        "[svzerodcalibrator] P/Q ChamberSphere calibration requires a 'time' "
        "vector.");
  }
  const Vec t = config["time"].get<Vec>();
  const size_t n = t.size();
  if (Pout.size() != n || Qin.size() != n || Qout.size() != n) {
    throw std::runtime_error(
        "[svzerodcalibrator] P/Q ChamberSphere calibration: observation and "
        "time lengths differ.");
  }
  Vec F(n);
  for (size_t i = 0; i < n; i++) F[i] = Qin[i] - Qout[i];
  const Vec Iflow = cumint(F, t);  // parameter-independent

  // Calibration parameters.
  const auto& cp = config["calibration_parameters"];
  const double tol_grad = cp.value("tolerance_gradient", 1e-6);
  const double tol_inc = cp.value("tolerance_increment", 1e-10);
  const int max_iter = cp.value("maximum_iterations", 100);
  double lambda = cp.value("initial_damping_factor", 1.0);
  const double period = cp.value("cardiac_cycle_period", -1.0);

  // Initial parameter vector and the active (calibrated) subset.
  double alpha[12];
  const auto& vals = (*vessel)["zero_d_element_values"];
  for (int k = 0; k < 12; k++) alpha[k] = vals.value(kParamNames[k], 0.0);

  std::set<std::string> calib;
  if (vessel->contains("calibrate")) {
    for (const auto& s : (*vessel)["calibrate"])
      calib.insert(s.get<std::string>());
  }
  std::vector<int> active;
  for (int k = 0; k < 12; k++)
    if (calib.count(kParamNames[k])) active.push_back(k);
  if (active.empty()) {
    throw std::runtime_error(
        "[svzerodcalibrator] No parameters selected for calibration. Add a "
        "'calibrate' list to the ChamberSphere vessel.");
  }
  const int na = static_cast<int>(active.size());

  // Levenberg-Marquardt loop on the active-stress residual G.
  Eigen::VectorXd vec_old;
  for (int it = 0; it < max_iter; it++) {
    Recon rc = reconstruct(alpha, t, Pout, F, Iflow, period);

    // Assemble Jacobian J (n x na) of G with respect to active parameters.
    Eigen::MatrixXd J(n, na);
    for (int c = 0; c < na; c++) {
      const int p = active[c];
      const Vec dtau = dtau_dparam(p, alpha, t, Pout, F, rc);
      Vec dG(n, 0.0);
      if (p == SIGMA) {
        for (size_t i = 0; i < n; i++) dG[i] = -rc.act[i].act_plus;
      } else if (p >= AMAX) {
        for (size_t i = 0; i < n; i++)
          dG[i] = rc.act[i].dact[p] * rc.tau[i] -
                  alpha[SIGMA] * rc.act[i].dact_plus[p];
      } else {
        const Vec dditau = ddt(dtau, t);
        for (size_t i = 0; i < n; i++)
          dG[i] = dditau[i] + rc.act[i].act * dtau[i];
      }
      for (size_t i = 0; i < n; i++) J(i, c) = dG[i];
    }
    Eigen::Map<const Eigen::VectorXd> G(rc.G.data(), n);

    // Normal equations with Levenberg-Marquardt damping (mirrors
    // LevenbergMarquardtOptimizer::update_delta).
    Eigen::VectorXd vec = J.transpose() * G;
    if (it > 0 && vec_old.norm() > 0.0) lambda *= vec.norm() / vec_old.norm();
    vec_old = vec;
    Eigen::MatrixXd JtJ = J.transpose() * J;
    Eigen::MatrixXd mat = JtJ + lambda * Eigen::MatrixXd(JtJ.diagonal().asDiagonal());
    Eigen::VectorXd delta = mat.ldlt().solve(vec);

    for (int c = 0; c < na; c++) alpha[active[c]] -= delta[c];

    if (vec.norm() < tol_grad && delta.norm() < tol_inc) break;
  }

  // Write calibrated parameters back.
  for (auto& v : out["vessels"]) {
    if (v.value("zero_d_element_type", std::string()) == "ChamberSphere" &&
        v["vessel_name"] == name) {
      for (int k = 0; k < 12; k++) v["zero_d_element_values"][kParamNames[k]] = alpha[k];
    }
  }
  out.erase("y");
  out.erase("dy");
  out.erase("time");
  out.erase("calibration_parameters");
  return out;
}
