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

// Cumulative trapezoidal integral with I[0] = 0. This is the only numerical
// operation on the time series: the absolute cavity volume is the time-integral
// of the net flow (flow alone gives volume change). All time-derivatives are
// taken from the supplied ``ydot`` analytically, so no finite differences are
// used.
Vec cumint(const Vec& f, const Vec& t) {
  Vec out(f.size(), 0.0);
  for (size_t i = 1; i < f.size(); i++) {
    out[i] = out[i - 1] + 0.5 * (f[i] + f[i - 1]) * (t[i] - t[i - 1]);
  }
  return out;
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
  const double df_dtsys = (-0.5 * sech2s / steep) * Sm;
  const double df_dtdias = Sp * (0.5 * sech2d / steep);
  const double df_dsteep =
      (0.5 * sech2s * (-(tc - a[TSYS]) / (steep * steep))) * Sm +
      Sp * (0.5 * sech2d * ((tc - a[TDIAS]) / (steep * steep)));
  const double da = a[AMAX] - a[AMIN];
  double dact_t[12] = {0};
  dact_t[AMAX] = f;
  dact_t[AMIN] = 1.0 - f;
  dact_t[TSYS] = da * df_dtsys;
  dact_t[TDIAS] = da * df_dtdias;
  dact_t[STEEP] = da * df_dsteep;

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

// Reconstructed internal trajectory and the active-stress residual G. The wall
// inertia term (rho) is neglected, so stress = P*R/thick0 and every quantity is
// algebraic in (P, dP, R, velo, dvelo) with velo/dvelo obtained analytically
// from the net flow and its supplied derivative.
struct Recon {
  Vec R, velo, dvelo, tau, dtau, G;
  std::vector<Activation> act;
};

Recon reconstruct(const double* a, const Vec& t, const Vec& Pp, const Vec& dPp,
                  const Vec& F, const Vec& dF, const Vec& Iflow, double period) {
  const size_t n = t.size();
  const double r0 = a[RADIUS0], thick0 = a[THICK0], W1v = a[W1], W2v = a[W2],
               eta = a[ETA], sigma = a[SIGMA];
  Recon rc;
  rc.R.resize(n); rc.velo.resize(n); rc.dvelo.resize(n);
  rc.tau.resize(n); rc.dtau.resize(n); rc.G.resize(n); rc.act.resize(n);
  for (size_t i = 0; i < n; i++) {
    const double Vcav = (4.0 / 3.0) * M_PI * r0 * r0 * r0 + Iflow[i];
    const double R = std::cbrt(3.0 * Vcav / (4.0 * M_PI));
    const double velo = F[i] / (4.0 * M_PI * R * R);
    const double dvelo = dF[i] / (4.0 * M_PI * R * R) - 2.0 * velo * velo / R;
    const double Pi = Pp[i], dP = dPp[i];

    const double tau =
        Pi * R / thick0 -
        2 * R * eta * velo * (1 + 2 * pow(r0, 12) / pow(R, 12)) / pow(r0, 2) -
        (4 - 4 * pow(r0, 6) / pow(R, 6)) * (pow(R, 2) * W2v / pow(r0, 2) + W1v);
    const double dtau =
        -2 * R * velo *
            (4 * W2v * (1 - pow(r0, 6) / pow(R, 6)) +
             12 * pow(r0, 8) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) / pow(R, 8)) /
            pow(r0, 2) -
        eta * ((1 + 2 * pow(r0, 12) / pow(R, 12)) *
                   (2 * R * dvelo + 2 * pow(velo, 2)) / pow(r0, 2) -
               48 * pow(r0, 10) * pow(velo, 2) / pow(R, 12)) +
        (Pi * velo + R * dP) / thick0;

    rc.R[i] = R; rc.velo[i] = velo; rc.dvelo[i] = dvelo;
    rc.tau[i] = tau; rc.dtau[i] = dtau;
    const double tc = period > 0.0 ? std::fmod(t[i], period) : t[i];
    rc.act[i] = activation(tc, a);
    rc.G[i] = dtau + rc.act[i].act * tau - sigma * rc.act[i].act_plus;
  }
  return rc;
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

  // Read port observations and their supplied derivatives (ydot). P = outlet
  // (= inlet) pressure; net flow F = Qin - Qout.
  if (!config.contains("y") || !config.contains("dy") ||
      !config.contains("time")) {
    throw std::runtime_error(
        "[svzerodcalibrator] P/Q ChamberSphere calibration requires 'y', 'dy' "
        "and 'time'.");
  }
  const auto& y = config["y"];
  const auto& dy = config["dy"];
  auto get = [](const nlohmann::json& src, const std::string& key) {
    if (!src.contains(key)) {
      throw std::runtime_error(
          "[svzerodcalibrator] P/Q ChamberSphere calibration: missing "
          "observation '" + key + "'.");
    }
    return src[key].get<Vec>();
  };
  const Vec Pp = get(y, "pressure:" + outnode);
  const Vec Qin = get(y, "flow:" + innode);
  const Vec Qout = get(y, "flow:" + outnode);
  const Vec dPp = get(dy, "pressure:" + outnode);
  const Vec dQin = get(dy, "flow:" + innode);
  const Vec dQout = get(dy, "flow:" + outnode);
  const Vec t = config["time"].get<Vec>();
  const size_t n = t.size();
  if (Pp.size() != n || Qin.size() != n || Qout.size() != n ||
      dPp.size() != n || dQin.size() != n || dQout.size() != n) {
    throw std::runtime_error(
        "[svzerodcalibrator] P/Q ChamberSphere calibration: observation and "
        "time lengths differ.");
  }
  Vec F(n), dF(n);
  for (size_t i = 0; i < n; i++) {
    F[i] = Qin[i] - Qout[i];
    dF[i] = dQin[i] - dQout[i];
  }
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
  if (calib.count("rho")) {
    throw std::runtime_error(
        "[svzerodcalibrator] 'rho' cannot be calibrated from P/Q: the wall "
        "inertia term it scales is neglected, so it has no effect on the "
        "observed pressure/flow. Remove it from 'calibrate'.");
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

  // Levenberg-Marquardt loop on the active-stress residual G. The Jacobian is
  // fully analytic.
  Eigen::VectorXd vec_old;
  for (int it = 0; it < max_iter; it++) {
    Recon rc = reconstruct(alpha, t, Pp, dPp, F, dF, Iflow, period);
    const double r0 = alpha[RADIUS0], thick0 = alpha[THICK0], W1v = alpha[W1],
                 W2v = alpha[W2], eta = alpha[ETA], sigma = alpha[SIGMA];

    Eigen::MatrixXd J(n, na);
    for (int c = 0; c < na; c++) {
      const int p = active[c];
      for (size_t i = 0; i < n; i++) {
        const double R = rc.R[i], velo = rc.velo[i], dvelo = rc.dvelo[i];
        const double Pi = Pp[i], dP = dPp[i];
        const double act = rc.act[i].act;
        double tau_dp = 0.0, dtau_dp = 0.0, dG = 0.0;
        switch (p) {
          case THICK0:
            tau_dp = -Pi * R / pow(thick0, 2);
            dtau_dp = -(Pi * velo + R * dP) / pow(thick0, 2);
            dG = dtau_dp + act * tau_dp;
            break;
          case W1:
            tau_dp = -4 + 4 * pow(r0, 6) / pow(R, 6);
            dtau_dp = -24 * pow(r0, 6) * velo / pow(R, 7);
            dG = dtau_dp + act * tau_dp;
            break;
          case W2:
            tau_dp = -pow(R, 2) * (4 - 4 * pow(r0, 6) / pow(R, 6)) / pow(r0, 2);
            dtau_dp =
                -2 * R * velo * (4 + 8 * pow(r0, 6) / pow(R, 6)) / pow(r0, 2);
            dG = dtau_dp + act * tau_dp;
            break;
          case ETA:
            tau_dp =
                -2 * R * velo * (1 + 2 * pow(r0, 12) / pow(R, 12)) / pow(r0, 2);
            dtau_dp = -(1 + 2 * pow(r0, 12) / pow(R, 12)) *
                          (2 * R * dvelo + 2 * pow(velo, 2)) / pow(r0, 2) +
                      48 * pow(r0, 10) * pow(velo, 2) / pow(R, 12);
            dG = dtau_dp + act * tau_dp;
            break;
          case RADIUS0: {
            // Chain through R, velo, dvelo (all functions of r0).
            const double dR = r0 * r0 / (R * R);
            const double dvelo_dr0 = -2 * velo * r0 * r0 / (R * R * R);
            const double ddvelo_dr0 =
                -dF[i] * r0 * r0 / (2 * M_PI * pow(R, 5)) +
                10 * velo * velo * r0 * r0 / pow(R, 4);
            const double tau_dr0e =
                2 * pow(R, 2) * W2v * (4 - 4 * pow(r0, 6) / pow(R, 6)) /
                    pow(r0, 3) +
                4 * R * eta * velo * (1 + 2 * pow(r0, 12) / pow(R, 12)) /
                    pow(r0, 3) +
                24 * pow(r0, 5) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                    pow(R, 6) -
                48 * eta * pow(r0, 9) * velo / pow(R, 11);
            const double tau_dR =
                Pi / thick0 -
                2 * R * W2v * (4 - 4 * pow(r0, 6) / pow(R, 6)) / pow(r0, 2) -
                2 * eta * velo * (1 + 2 * pow(r0, 12) / pow(R, 12)) /
                    pow(r0, 2) -
                24 * pow(r0, 6) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                    pow(R, 7) +
                48 * eta * pow(r0, 10) * velo / pow(R, 12);
            const double tau_dvelo =
                -2 * R * eta * (1 + 2 * pow(r0, 12) / pow(R, 12)) / pow(r0, 2);
            const double dtau_dr0e =
                -2 * R * velo *
                    (-48 * W2v * pow(r0, 5) / pow(R, 6) +
                     96 * pow(r0, 7) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                         pow(R, 8)) /
                    pow(r0, 2) +
                4 * R * velo *
                    (4 * W2v * (1 - pow(r0, 6) / pow(R, 6)) +
                     12 * pow(r0, 8) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                         pow(R, 8)) /
                    pow(r0, 3) -
                eta * (-2 * (1 + 2 * pow(r0, 12) / pow(R, 12)) *
                           (2 * R * dvelo + 2 * pow(velo, 2)) / pow(r0, 3) -
                       480 * pow(r0, 9) * pow(velo, 2) / pow(R, 12) +
                       24 * pow(r0, 9) * (2 * R * dvelo + 2 * pow(velo, 2)) /
                           pow(R, 12));
            const double dtau_dR =
                -2 * R * velo *
                    (48 * W2v * pow(r0, 6) / pow(R, 7) -
                     96 * pow(r0, 8) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                         pow(R, 9)) /
                    pow(r0, 2) +
                dP / thick0 -
                eta * (2 * dvelo * (1 + 2 * pow(r0, 12) / pow(R, 12)) /
                           pow(r0, 2) +
                       576 * pow(r0, 10) * pow(velo, 2) / pow(R, 13) -
                       24 * pow(r0, 10) * (2 * R * dvelo + 2 * pow(velo, 2)) /
                           pow(R, 13)) -
                2 * velo *
                    (4 * W2v * (1 - pow(r0, 6) / pow(R, 6)) +
                     12 * pow(r0, 8) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                         pow(R, 8)) /
                    pow(r0, 2);
            const double dtau_dvelo =
                Pi / thick0 -
                2 * R *
                    (4 * W2v * (1 - pow(r0, 6) / pow(R, 6)) +
                     12 * pow(r0, 8) * (pow(R, 2) * W2v / pow(r0, 2) + W1v) /
                         pow(R, 8)) /
                    pow(r0, 2) -
                eta * (4 * velo * (1 + 2 * pow(r0, 12) / pow(R, 12)) /
                           pow(r0, 2) -
                       96 * pow(r0, 10) * velo / pow(R, 12));
            const double dtau_ddvelo =
                -2 * R * eta * (1 + 2 * pow(r0, 12) / pow(R, 12)) / pow(r0, 2);
            tau_dp = tau_dr0e + tau_dR * dR + tau_dvelo * dvelo_dr0;
            dtau_dp = dtau_dr0e + dtau_dR * dR + dtau_dvelo * dvelo_dr0 +
                      dtau_ddvelo * ddvelo_dr0;
            dG = dtau_dp + act * tau_dp;
            break;
          }
          case SIGMA:
            dG = -rc.act[i].act_plus;
            break;
          default:  // AMAX, AMIN, TSYS, TDIAS, STEEP
            dG = rc.act[i].dact[p] * rc.tau[i] -
                 sigma * rc.act[i].dact_plus[p];
            break;
        }
        J(i, c) = dG;
      }
    }
    Eigen::Map<const Eigen::VectorXd> G(rc.G.data(), n);

    // Normal equations with Levenberg-Marquardt damping (mirrors
    // LevenbergMarquardtOptimizer::update_delta).
    Eigen::VectorXd vec = J.transpose() * G;
    if (it > 0 && vec_old.norm() > 0.0) lambda *= vec.norm() / vec_old.norm();
    vec_old = vec;
    Eigen::MatrixXd JtJ = J.transpose() * J;
    Eigen::MatrixXd mat =
        JtJ + lambda * Eigen::MatrixXd(JtJ.diagonal().asDiagonal());
    Eigen::VectorXd delta = mat.ldlt().solve(vec);

    for (int c = 0; c < na; c++) alpha[active[c]] -= delta[c];

    if (vec.norm() < tol_grad && delta.norm() < tol_inc) break;
  }

  // Write calibrated parameters back.
  for (auto& v : out["vessels"]) {
    if (v.value("zero_d_element_type", std::string()) == "ChamberSphere" &&
        v["vessel_name"] == name) {
      for (int k = 0; k < 12; k++)
        v["zero_d_element_values"][kParamNames[k]] = alpha[k];
    }
  }
  out.erase("y");
  out.erase("dy");
  out.erase("time");
  out.erase("calibration_parameters");
  return out;
}
