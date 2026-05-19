// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using Approx = Catch::Approx;

#include "C2HDM.h"
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/NumericalDerivatives.h>
namespace
{

const std::vector<double> example_point_C2HDM{/* lambda_1 = */ 3.29771,
                                              /* lambda_2 = */ 0.274365,
                                              /* lambda_3 = */ 4.71019,
                                              /* lambda_4 = */ -2.23056,
                                              /* Re(lambda_5) = */ -2.43487,
                                              /* Im(lambda_5) = */ 0.124948,
                                              /* Re(m_{12}^2) = */ 2706.86,
                                              /* tan(beta) = */ 4.64487,
                                              /* Yukawa Type = */ 1};

} // namespace

/**
 * @test Check if the automatic Debye corrections match the SM one in the SM
 * case. This should be y_t^2/4.
 */
TEST_CASE("Test Calculate Debye", "[origin]")
{

  using namespace BSMPT;
  ISMConstants SM;
  SM.C_MassTop = 172;
  SM.C_vev0    = 246;

  // This is not a legal point but just a dummy point to only have the top
  // coupling in the Debye Contributions
  const std::vector<double> example_point_CXSM{/* vh = */ SM.C_vev0,
                                               /* vs = */ 0,
                                               /* va = */ 0,
                                               /* ms = */ 41.67,
                                               /* lambda = */ 0,
                                               /* delta2 = */ 0,
                                               /* b2 = */ 0,
                                               /* d2 = */ 0,
                                               /* Reb1 = */ 0,
                                               /* Imb1 = */ 0,
                                               /* Rea1 = */ 0,
                                               /* Ima1 = */ 0};

  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SM);
  modelPointer->initModel(example_point_CXSM);
  modelPointer->CalculateDebye(true);
  auto debye = modelPointer->get_DebyeHiggs();

  double topCoupling = SM.C_MassTop * std::sqrt(2) / SM.C_vev0;
  double expected    = std::pow(topCoupling, 2) / 4.0;

  auto calculated = debye.at(3).at(3);

  REQUIRE(calculated != 0);

  REQUIRE(calculated == Approx(expected).margin(1e-4));
}

TEST_CASE("Check f_{abcd}", "[origin]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);
  modelPointer->resetScale(200);

  {
    double expected = -0.30575009e-3;
    double result   = modelPointer->fbaseFour(10, 20, 30, 40);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }
  {
    double expected = -3.518370054e-7;
    double result   = modelPointer->fbaseFour(1e4, 200, 300, 5);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -0.23248145e-3;
    double result   = modelPointer->fbaseFour(20, 20, 30, 40);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -0.284264097e-3;
    double result   = modelPointer->fbaseFour(20, 20, 20, 40);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -0.4166666667e-3;
    double result   = modelPointer->fbaseFour(20, 20, 20, 20);
    REQUIRE(result == Approx(expected).margin(1e-4));
  }

  {
    double expected = -3.696784963e-9;
    double result   = modelPointer->fbaseFour(
        0, std::pow(100, 2), std::pow(200, 2), std::pow(50, 2));
    REQUIRE(result == Approx(expected).margin(1e-4));
  }
}

TEST_CASE("Check potential int -> Order cast", "[origin]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  modelPointer->initModel(example_point_C2HDM);

  const std::vector<double> vev =
      modelPointer->MinimizeOrderVEV({0., 20.4, 30.7, 40.2}); // random
  const std::vector<double> T_list    = {0, 10, 50, 100, 300};
  const std::vector<double> diff_list = {0};

  for (const auto &T : T_list)
    for (const auto &diff : diff_list)
    {
      REQUIRE(modelPointer->VEff(vev, T, diff, 0) ==
              Approx(modelPointer->VEff(vev, T, diff, Order::TreeLevel))
                  .margin(1e-8));
      REQUIRE(modelPointer->VEff(vev, T, diff, 1) ==
              Approx(modelPointer->VEff(vev, T, diff, Order::OneLoop))
                  .margin(1e-8));
      // In case we add 2-loop we get a warning from here
      REQUIRE(modelPointer->VEff(vev, T, diff, 2) ==
              Approx(modelPointer->VEff(vev, T, diff, Order::OneLoop))
                  .margin(1e-8));
    }
}

TEST_CASE("Check of analytical T-derivative of Veff (CxSM example point)",
          "[origin]")
{
  const std::vector<double> example_point_CXSM{/* v = */ 245.34120667410863,
                                               /* vs = */ 0,
                                               /* va = */ 0,
                                               /* msq = */ -15650,
                                               /* lambda = */ 0.52,
                                               /* delta2 = */ 0.55,
                                               /* b2 = */ -8859,
                                               /* d2 = */ 0.5,
                                               /* Reb1 = */ 0,
                                               /* Imb1 = */ 0,
                                               /* Rea1 = */ 0,
                                               /* Ima1 = */ 0};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SMConstants);
  modelPointer->initModel(example_point_CXSM);

  std::vector<double> T_values{20., 50., 100., 200., 500., 1000.};
  std::vector<std::vector<double>> phi_values{{0., 0., 0.},
                                              {100., 0., 0.},
                                              {0., 100., 0.},
                                              {0., 0., 100.},
                                              {100., 200., 0.},
                                              {50., 100., 200.}};

  const double eps = 1e-3;
  std::function<double(std::vector<double>, double)> dTV_num =
      [=](auto const &vev, auto const &T)
  {
    return NablaNumerical(
               {T},
               [&](std::vector<double> Tv)
               {
                 // Potential wrapper
                 return modelPointer->VEff(modelPointer->MinimizeOrderVEV(vev),
                                           Tv.at(0));
               },
               eps)
        .at(0);
  };

  for (double T_val : T_values)
  {
    for (auto &phi_val : phi_values)
    {
      double x_num = dTV_num(phi_val, T_val);
      double x_an  = modelPointer->VEff(
          modelPointer->MinimizeOrderVEV(phi_val), T_val, -1);
      REQUIRE(x_num == Approx(x_an).epsilon(5e-3));
    }
  }
}

TEST_CASE("Check of analytical T-derivative of Veff (R2HDM example point)",
          "[origin]")
{
  const std::vector<double> example_point_R2HDM{
      /* lambda_1 = */ 6.9309437685026,
      /* lambda_2 = */ 0.26305141403285998,
      /* lambda_3 = */ 1.2865950045595,
      /* lambda_4 = */ 4.7721306931875001,
      /* lambda_5 = */ 4.7275722046239004,
      /* m_{12}^2 = */ 18933.440789693999,
      /* tan(beta) = */ 16.577896825227999,
      /* Yukawa Type = */ 1};

  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  modelPointer->initModel(example_point_R2HDM);

  std::vector<double> T_values{20., 50., 100., 200., 500., 1000.};
  std::vector<std::vector<double>> phi_values{{0., 0., 0., 0.},
                                              {100., 0., 0., 0.},
                                              {0., 100., 0., 0.},
                                              {0., 0., 100., 0.},
                                              {0., 50., 200., 0.},
                                              {50., 100., 200., 10.}};

  const double eps = 1e-3;
  std::function<double(std::vector<double>, double)> dTV_num =
      [=](auto const &vev, auto const &T)
  {
    return NablaNumerical(
               {T},
               [&](std::vector<double> Tv)
               {
                 // Potential wrapper
                 return modelPointer->VEff(modelPointer->MinimizeOrderVEV(vev),
                                           Tv.at(0));
               },
               eps)
        .at(0);
  };

  for (double T_val : T_values)
  {
    for (auto &phi_val : phi_values)
    {
      double x_num = dTV_num(phi_val, T_val);
      double x_an  = modelPointer->VEff(
          modelPointer->MinimizeOrderVEV(phi_val), T_val, -1);
      REQUIRE(x_num == Approx(x_an).epsilon(5e-3));
    }
  }
}
