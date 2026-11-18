// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ChamberSphereThick.h
 * @brief model::ChamberSphereThick source file
 */
#ifndef SVZERODSOLVER_MODEL_ChamberSphereThick_HPP_
#define SVZERODSOLVER_MODEL_ChamberSphereThick_HPP_

#include <math.h>

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "ActivationFunction.h"
#include "ActiveStress.h"
#include "Block.h"
#include "Material.h"
#include "SparseSystem.h"

/**
 * @brief Thick-walled spherical heart chamber model
 *
 * Reduced-order ventricular model on a thick-walled sphere of reference
 * internal radius \f$R_i\f$ and external radius \f$R_e\f$. It is the
 * thick-walled counterpart of \ref ChamberSphere (\cite caruel13), which
 * neglects the radial stress and assumes the wall state to be uniform across
 * the thickness, and the spherical counterpart of \ref ChamberCylinder
 * (\cite genet23). Like both, the block models the ventricular wall only and
 * exposes a single cavity pressure at its inlet/outlet; valves and the
 * circulation are separate 0D blocks.
 *
 * Compared to \ref ChamberSphere, this block carries the radial stress and
 * resolves the stretch, the stress and the state of the active stress law
 * **through the wall thickness**, so it stays valid for wall-thickness ratios
 * at which the shell assumption of \ref ChamberSphere breaks down.
 *
 * ### Kinematics
 *
 * The motion is spherically symmetric, \f$r = r(R)\f$, \f$\theta = \Theta\f$,
 * \f$\varphi = \Phi\f$, with the deformation gradient
 * \f$F = \mathrm{diag}(r', r/R, r/R)\f$ in the orthonormal
 * \f$(\mathbf{e}_R,\mathbf{e}_\Theta,\mathbf{e}_\Phi)\f$ basis and the volume
 * ratio \f$J = r'(r/R)^2\f$. Incompressibility \f$J\equiv1\f$ integrates
 * **exactly** to
 * \f[
 * r(R)^3 = R^3 + r_i^3 - R_i^3,
 * \f]
 * so the whole wall is determined by the single unknown
 * \f$\rho = r_i - R_i\f$, the endocardial radial displacement, at any wall
 * thickness. This is not a kinematic approximation: it is the exact
 * incompressible spherically symmetric motion. Because every admissible
 * variation keeps \f$J = 1\f$ pointwise (\f$\delta J\equiv0\f$), the
 * hydrostatic pressure does no virtual work and never enters the equations,
 * so no mixed u/p formulation and no pressure unknown are needed (unlike
 * \ref ChamberCylinder).
 *
 * With \f$r = (R^3 + r_i^3 - R_i^3)^{1/3}\f$ the circumferential stretch is
 * \f$\lambda(R) = r/R\f$, the radial stretch is \f$\lambda^{-2}\f$ and
 * \f[
 * C = \mathrm{diag}(\lambda^{-4}, \lambda^{2}, \lambda^{2}),\quad
 * \frac{\partial r}{\partial\rho} = \frac{r_i^2}{r^2},\quad
 * \lambda_{,\rho} = \frac{r_i^2}{R\,r^2}.
 * \f]
 *
 * ### Myofibers
 *
 * Spherical symmetry requires the myofibers to be distributed isotropically
 * in the wall plane (\cite caruel13 Remark 1); a single helix angle would
 * break it, which is why \cite genet23 moved to a cylinder to resolve the
 * fiber orientation and ventricular twist. Since \f$C\f$ has equal in-plane
 * eigenvalues, \f$I_4 = \mathbf{e}_F\cdot C\cdot\mathbf{e}_F = \lambda^2\f$
 * for **any** in-plane fiber direction, so the passive \ref Material is
 * evaluated with \f$\mathbf{e}_F = \mathbf{e}_\Theta\f$ without loss of
 * generality, and the active fiber stress homogenizes to the in-plane average
 * \f$\langle\mathbf{e}_F\otimes\mathbf{e}_F\rangle =
 * \tfrac12(\mathbf{e}_\Theta\otimes\mathbf{e}_\Theta +
 * \mathbf{e}_\Phi\otimes\mathbf{e}_\Phi)\f$. This is the same convention as
 * \ref ChamberSphere, whose scalar spherical stress is the sum of the two
 * orthoradial components.
 *
 * ### Constitutive law
 *
 * The second Piola-Kirchhoff stress is
 * \f[
 * \Sigma = \Sigma^d + \gamma\,\dot E
 * + \frac{\tau}{2}\left(\mathbf{e}_\Theta\otimes\mathbf{e}_\Theta
 * + \mathbf{e}_\Phi\otimes\mathbf{e}_\Phi\right),
 * \f]
 * with the deviatoric stress \f$\Sigma^d = 2\,\partial W/\partial C\f$ of the
 * wall \ref Material (`material` object), a viscous term \f$\gamma\dot E\f$
 * and the active fiber stress \f$\tau\f$ of the \ref ActiveStress law
 * (`active_stress` object), evaluated **at every quadrature point** at the
 * local fiber stretch \f$\lambda(R)\f$ and driven by the
 * `activation_function`. Each quadrature point therefore carries its own copy
 * of the law's internal state, which is what makes the transmural
 * heterogeneity of the contraction visible.
 *
 * ### Governing equations
 *
 * The dynamic principle of virtual work on the one-parameter family above,
 * with \f$\mathbf{u} = (r-R)\mathbf{e}_R\f$,
 * \f$\mathrm{d}\Omega = 4\pi R^2\mathrm{d}R\f$ and the cavity volume
 * \f$V = \tfrac43\pi r_i^3\f$, reads
 * \f[
 * 4\pi\int_{R_i}^{R_e}\left[\varrho_0\left(
 * \left(\frac{\partial r}{\partial\rho}\right)^2\dot v
 * + \frac{\partial r}{\partial\rho}
 * \frac{\partial^2 r}{\partial\rho^2}v^2\right)
 * + \Sigma:\frac{\partial E}{\partial\lambda}\,\lambda_{,\rho}
 * \right]R^2\,\mathrm{d}R = 4\pi P_v r_i^2,
 * \f]
 * where \f$v = \dot\rho\f$ is a velocity companion unknown and the second
 * inertia term is the convective part
 * \f$\tfrac12\,\mathrm{d}M/\mathrm{d}\rho\;v^2\f$ of the configuration
 * dependent mass \f$M(\rho) = \varrho_0\int(\partial r/\partial\rho)^2
 * \mathrm{d}\Omega\f$. The integral is evaluated with composite three-point
 * Gauss quadrature over `num_elements` sub-intervals of \f$[R_i,R_e]\f$. The
 * wall equation is complemented by the active stress law at each quadrature
 * point and by
 * \f[
 * V - \tfrac43\pi r_i^3 = 0,\quad
 * Q_\text{in} - Q_\text{out} - \dot V = 0,\quad
 * P_\text{in} - P_\text{out} = P_v,\quad
 * v - \dot\rho = 0.
 * \f]
 *
 * ### Relation to ChamberSphere
 *
 * As \f$d_0 = R_e - R_i \to 0\f$ the wall equation reduces term by term to
 * \ref ChamberSphere: the integral becomes \f$4\pi R_0^2 d_0\,(\cdot)\f$ and
 * \f$\lambda_{,\rho}\to1/R_0\f$, so dividing by \f$4\pi R_0^2\f$ gives
 * \f$\varrho_0 d_0\dot v + (d_0/R_0)\,\Sigma:\partial E/\partial\lambda
 * - P_v\lambda^2 = 0\f$, which is \cite caruel13 Eq. (19) with the spherical
 * stress \f$S = \lambda^{-1}\Sigma:\partial E/\partial\lambda\f$ of
 * \ref ChamberSphere. The viscosity \f$\gamma\f$ of this block is the
 * \f$\eta\f$ of \ref ChamberSphere: the term
 * \f$\eta\dot C(1+2C^{-6})\f$ of \ref ChamberSphere is exactly
 * \f$\lambda^{-1}(\gamma\dot E):\partial E/\partial\lambda\f$. A
 * \ref ChamberSphere model of radius \f$r_0\f$ and thickness \f$d_0\f$
 * therefore corresponds to \f$R_i = r_0 - d_0/2\f$, \f$R_e = r_0 + d_0/2\f$
 * here, and the two agree to \f$O((d_0/r_0)^2)\f$.
 *
 * ### Tangent
 *
 * The block assembles exact derivatives. With \f$A = \partial E/\partial
 * \lambda = \mathrm{diag}(-2\lambda^{-5},\lambda,\lambda)\f$,
 * \f$B = \partial^2E/\partial\lambda^2 =
 * \mathrm{diag}(10\lambda^{-6},1,1)\f$, the wall integrand
 * \f$g = \Sigma:A\,\lambda_{,\rho}\f$ has
 * - \f$\partial g/\partial\tau = \lambda\,\lambda_{,\rho}\f$,
 * - \f$\partial g/\partial v = \gamma\,\lambda_{,\rho}^2\,(A:A)\f$,
 * - \f$\partial g/\partial\rho = \lambda_{,\rho}\left[\left(
 *   \mathrm{D}\Sigma^d\lambda_{,\rho}
 *   + \gamma v(B\lambda_{,\rho}^2 + A\lambda_{,\rho\rho})\right):A
 *   + \Sigma:B\,\lambda_{,\rho}\right]
 *   + \Sigma:A\;\lambda_{,\rho\rho}\f$,
 *
 * where \f$\mathrm{D}\Sigma^d\f$ is the \ref Material derivative along
 * \f$\delta C = 2A\f$, \f$\delta J = 0\f$. The active stress law depends on
 * the wall through \f$\partial\lambda/\partial\rho = \lambda_{,\rho}\f$. With
 * \f$u = r_i^3 - R_i^3\f$ the radial derivatives are
 * \f$r_{,\rho} = r_i^2/r^2\f$,
 * \f$r_{,\rho\rho} = 2r_i/r^2 - 2r_i^4/r^5\f$ and
 * \f$r_{,\rho\rho\rho} = 10r_i^6/r^8 - 12r_i^3/r^5 + 2/r^2\f$, and
 * \f$\lambda_{,\rho\cdots} = r_{,\rho\cdots}/R\f$.
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block:
 *
 * * `Ri` - Reference internal (endocardial) radius \f$R_i\f$
 * * `Re` - Reference external (epicardial) radius \f$R_e\f$
 * * `gamma` - Material viscosity \f$\gamma\f$ (the `eta` of
 *   \ref ChamberSphere)
 * * `num_elements` - Number of sub-intervals of the three-point Gauss
 *   quadrature through the wall (optional, default 10)
 * * `density` - Reference mass density \f$\varrho_0 \ge 0\f$ (optional,
 *   default 1000 kg/m^3); 0 gives a quasi-static wall
 *
 * Like \ref ChamberSphere, the block is declared in `chambers` with these
 * parameters in `values`, a `material` object (see \ref Material), an
 * `active_stress` object (see \ref ActiveStress) and an
 * `activation_function` object (see \ref ActivationFunction).
 *
 * ### Usage in json configuration file
 *
 *     "chambers": [
 *        {
 *            "type": "ChamberSphereThick",
 *            "name": "ventricle",
 *            "values": {
 *                "Ri" : 0.025,
 *                "Re" : 0.035,
 *                "gamma" : 10.0,
 *                "num_elements" : 10,
 *                "density" : 1e3
 *            },
 *            "activation_function": {},
 *            "material": {},
 *            "active_stress": {}
 *        }
 *     ]
 *
 * ### Internal variables
 *
 * Names of internal variables in this block's output:
 *
 * * `radius` - Endocardial radial displacement \f$\rho = r_i - R_i\f$
 * * `<state>_<q>` - States of the active stress law at quadrature point
 *   \f$q\f$ (state 0 is the active stress \f$\tau\f$), ordered from the
 *   endocardium (\f$q = 0\f$) to the epicardium
 * * `volume` - Cavity volume \f$V\f$
 * * `velo` - Endocardial radial velocity \f$v = \dot\rho\f$
 *
 * `volume` is defined algebraically, \f$V = \tfrac43\pi r_i^3\f$, and is a
 * pure output: mass conservation uses \f$\dot V = 4\pi r_i^2 v\f$ from the
 * velocity companion, so the volume unknown feeds back nowhere. Its default
 * initial value 0 is inconsistent with \f$\rho = 0\f$, and generalized-\f$
 * \alpha\f$ then rings it down over the first few steps (a factor
 * \f$1-1/\alpha_f\f$ per step) before it locks onto the geometric volume.
 * Give it the consistent `initial_condition`
 * \f$\texttt{"volume:<name>"} = \tfrac43\pi R_i^3\f$ to avoid that startup
 * transient.
 *
 */
class ChamberSphereThick : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    Ri = 0,
    Re = 1,
    gamma = 2,
    num_elements = 3,
    density = 4,
  };

  /**
   * @brief Construct a new ChamberSphereThick object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ChamberSphereThick(int id, Model* model)
      : Block(id, model, BlockType::chamber_sphere_thick, BlockClass::chamber,
              {{"Ri", InputParameter()},
               {"Re", InputParameter()},
               {"gamma", InputParameter()},
               {"num_elements", InputParameter(true, false, true, 10.0)},
               {"density", InputParameter(true, false, true, 1000.0)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem& system,
                   std::vector<double>& parameters) override;

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivate of the solution
   */
  void update_solution(
      SparseSystem& system, std::vector<double>& parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) override;

  /**
   * @brief Number of triplets of element
   */
  TripletsContributions num_triplets;

  /**
   * @brief Get number of triplets of element
   *
   * @return TripletsContributions Number of triplets of element
   */
  TripletsContributions get_num_triplets() override { return num_triplets; }

  /**
   * @brief Set the activation function driving the active stress law.
   *
   * @param af Unique pointer to the activation function
   */
  void set_activation_function(
      std::unique_ptr<ActivationFunction> af) override {
    activation_function_ = std::move(af);
  }

  /**
   * @brief Set the wall material.
   *
   * @param m Unique pointer to the material
   */
  void set_material(std::unique_ptr<Material> m) override {
    material_ = std::move(m);
  }

  /**
   * @brief Set the active stress law.
   *
   * @param as Unique pointer to the active stress law
   */
  void set_active_stress(std::unique_ptr<ActiveStress> as) override {
    active_stress_ = std::move(as);
  }

 private:
  // --- Quadrature of the wall thickness ---
  int n_quad = 0;  ///< Number of quadrature points through the wall
  int n_law = 0;   ///< Number of active stress states per quadrature point
  int n_eqn = 0;   ///< Number of block equations

  /// Check the geometry, density and that material, active stress law and
  /// activation function are set and compatible.
  void check_inputs() const;

  /**
   * @brief A single Gauss quadrature point in the thickness
   */
  struct QuadPoint {
    double R;  ///< Reference radius of the quadrature point
    double w;  ///< Integration weight (Gauss weight * dR/dxi * R^2)
  };
  std::vector<QuadPoint> quad;  ///< Quadrature points across the thickness

  // --- Local variable / equation index helpers ---
  // Layout:
  //   vars: [Pin, Qin, Pout, Qout, radius, (active stress states)_q..,
  //          volume, velo]
  //   eqns: [wall, (active stress equations)_q.., volume, mass, pressure,
  //          velo]
  int i_rho() const { return 4; }
  /// Start of the active stress variables
  int i_active0() const { return 5; }
  /// State i of the active stress law at quadrature point q (state 0: tau)
  int i_law(int q, int i) const { return i_active0() + n_law * q + i; }
  int n_active_var() const { return n_law * n_quad; }
  int i_vol() const { return i_active0() + n_active_var(); }
  int i_velo() const { return i_vol() + 1; }

  int e_wall() const { return 0; }
  /// Start of the active stress equations
  int e_active0() const { return 1; }
  /// Equation i of the active stress law at quadrature point q
  int e_law(int q, int i) const { return e_active0() + n_law * q + i; }
  int e_volume() const { return e_active0() + n_active_var(); }
  int e_mass() const { return e_volume() + 1; }
  int e_pressure() const { return e_volume() + 2; }
  int e_velo() const { return e_volume() + 3; }

  /// Wall material
  std::unique_ptr<Material> material_;

  /// Active stress law
  std::unique_ptr<ActiveStress> active_stress_;

  /// Activation function driving the active stress law
  std::unique_ptr<ActivationFunction> activation_function_;
};

#endif  // SVZERODSOLVER_MODEL_ChamberSphereThick_HPP_
