// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 * Model file for the SM + real singlet with the tree-level potential
 */

#pragma once

#include <string> // for string
#include <vector> // for vector

#include <BSMPT/models/ClassPotentialOrigin.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

/**
 * @brief The Class_Potential_RxSM class
 * Implementation of the RxSM according to https://arxiv.org/pdf/1512.05355 Eq.
 * (11)
 *
 * * \f$ -L_S = (msq / 2) H^\dagger H + (\lambda / 4) (H^\dagger H)^2 +
 * (\lambda_{HS} / 2) (H^\dagger H) S^2 + (mSsq / 2) S^2 + (\lambda_S / 24) S^4
 * \f$
 * * \f$ -L_Y = \overline{U}_L V \text{diag}(y_d,y_s,y_b) D_R H^+  +
 * \overline{D}_L \text{diag}(y_d,y_s,y_b) D_R H^0 + \overline{U}_L
 * \text{diag}(y_u, y_c, y_t) U_R (H^0)^\ast - \overline{D}_L V^\dagger
 * \text{diag}(y_u, y_c, y_t) U_R H^- + \overline{E}_L \text{diag}(y_e, y_\mu,
 * y_\tau) E_R H^0 +\overline{\nu}_L \text{diag}(y_e, y_\mu, y_\tau) E_R H^+ +
 * c.c.  \f$
 * * \f$ L_G = (D_\mu H)^\dagger (D^\mu H) \f$
 *
 * with
 *
 * * \f$ H = \begin{pmatrix} H^+ \\  H^0 \end{pmatrix} = 1/\sqrt{2}
 * \begin{pmatrix} higgsfield[0] + I*higgsfield[1] \\  higgsfield[3] + I
 * *higgsfield[2] \end{pmatrix}\,, \f$
 * * \f$ S = higgsfield[4] \f$ with higgsfield = [rho1,eta1,psi1,zeta1,zetaS]
 * * \f$ D_\mu = -I C\_{}g/2 * W_\mu^a \sigma^a -I C\_{}gs/2 B_\mu =
 * -\frac{I}{2} \begin{pmatrix} C\_{}gs B + C\_{}g W3  & C\_{}g (W1 -I W2) \\
 * C\_{}g (W1 +I W2) & C\_{}gs B - C\_{}g W3 \end{pmatrix} \f$ with the gauge
 * base = [W1,W2,W3,B]
 *
 * * \f$ \overline{U}_L = \begin{pmatrix} \overline{u}_L &  \overline{c}_L &
 * \overline{t}_L \end{pmatrix}\,, \f$
 * * \f$ U_R = \begin{pmatrix} u_R \\ c_R \\ t_R \end{pmatrix}\,, \f$
 * * \f$ \overline{D}_L = \begin{pmatrix} \overline{d}_L & \overline{s}_L &
 * \overline{b}_L \end{pmatrix}\,, \f$
 * * \f$ D_R = \begin{pmatrix} d_R \\ s_R \\ b_R \end{pmatrix}\,, \f$
 * * \f$ \overline{E}_L = \begin{pmatrix} \overline{e}_L & \overline{\mu}_L &
 * \overline{\tau}_L \end{pmatrix}\,, \f$
 * * \f$ E_R = \begin{pmatrix} e_R \\ \mu_R \\ \tau_R \end{pmatrix}\,, \f$
 * * \f$ \overline{\nu}_L = \begin{pmatrix} \overline{\nu}_{e,L} &
 * \overline{\nu}_{\mu,L} & \overline{\nu}_{\tau,L} \end{pmatrix} \f$
 *
 * The basis of the quarks given by \f$ [u_R, c_R, t_R, d_R, s_R, b_R,
 * \overline{u}_L, \overline{c}_L, \overline{t}_L, \overline{d}_L,
 * \overline{s}_L, \overline{b}_L ] \f$ and for the leptons \f$ [\overline{e}_L,
 * e_R, \overline{\mu}_L , \mu_R, \overline{\tau}_L, \tau_R,
 * \overline{\nu}_{e,L}, \overline{\nu}_{\mu,L}, \overline{\nu}_{\tau,L} ] \f$
 *
 */
class Class_Potential_RxSM : public Class_Potential_Origin
{
public:
  Class_Potential_RxSM(const ISMConstants &smConstants);
  virtual ~Class_Potential_RxSM();

  // Choice of parameters of Lagrangian from https://arxiv.org/pdf/1512.05355
  // Eq. (11)
  double lambdaS, lambdaHS, vS;

  // Not an input parameter; lambda is fixed via the requirement of having
  // one of the Higgs bosons as the SM one with 125.09 GeV
  double lambda;

  // Not an input parameter; set to the SM value
  double vH;

  // Not input parameters; set through the tadpole equations
  double msq, mSsq;

  double alpha;

  bool UnbrokenSingletPhase;

  double dmsq, dlambda, dmSsq, dlambdaS, dlambdaHS, dT1, dT2, dT3, dT4, dT5;

  std::size_t pos_Gp, pos_Gm, pos_G0, pos_h, pos_H;
  std::size_t pos_h_SM, pos_h_H;

  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  /*
   * RxSM interaction basis:
   * 0     1     2     3      4
   * rho1, eta1, psi1, zeta1, zetaS
   */
  const std::size_t pos_rho1 = 0, pos_eta1 = 1, pos_psi1 = 2, pos_zeta1 = 3,
                    pos_zetaS = 4;

  /**
   * Helper function to determine mass indices of rotation matrix
   * @param HiggsMasses : vector with squared Higgs masses allocated
   *                      in AdjustRotationMatrix
   * @param HiggsRot : rotation matrix from interaction to mass basis
   *                   as calculated in AdjustRotationMatrix
   */
  void FindMassBasisIndices(const std::vector<double> &HiggsMasses,
                            const MatrixXd &HiggsRot);

  void AdjustRotationMatrix() override;
  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
