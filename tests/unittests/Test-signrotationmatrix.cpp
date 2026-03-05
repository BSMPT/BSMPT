// SPDX-FileCopyrightText: 2026 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <BSMPT/minimizer/Minimizer.h>
#include <BSMPT/models/ClassPotentialC2HDM.h>       // for Class_Potential_C2HDM
#include <BSMPT/models/ClassPotentialCPintheDark.h> // for Class_Potential_CPintheDark
#include <BSMPT/models/ClassPotentialCxSM.h>        // for Class_CxSM
#include <BSMPT/models/ClassPotentialN2HDM.h>       // for Class_Potential_N2HDM
#include <BSMPT/models/ClassPotentialOrigin.h> // for Class_Potential_Origin
#include <BSMPT/models/ClassPotentialR2HDM.h>  // for Class_Potential_R2HDM
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/models/modeltests/ModelTestfunctions.h>
#include <BSMPT/utility/Logger.h>
#include <fstream>

#include <iostream>

using Approx = Catch::Approx;

TEST_CASE("Checking sign convention of rotation matrix for C2HDM",
          "[signrotationmatrix]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::C2HDM, SMConstants);
  auto nPar = modelPointer->get_nPar();

  std::ifstream testData(
      TEST_DATA_PATH
      "/signrotationmatrix_data/c2hdm_unittest_signrotationmatrix.tsv");

  auto Check = [](auto result, auto expected)
  {
    if (std::abs(expected) > 1e-4)
    {
      REQUIRE(result == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  std::string line;
  std::getline(testData, line); // header, ignore

  while (std::getline(testData, line))
  {
    // Read in test data and parse properly: first the input parameters, then
    // the rotation matrix elements
    std::stringstream ss(line);
    std::vector<double> parameter_point;
    while (!ss.eof() && parameter_point.size() < nPar)
    {
      double x;
      ss >> x;
      parameter_point.push_back(x);
    }

    std::vector<double> rotation_matrix_angles;
    while (!ss.eof())
    {
      double x;
      ss >> x;
      rotation_matrix_angles.push_back(x);
    }

    modelPointer->initModel(parameter_point);

    // Since we need to access model-specific variables, cast modelPointer to
    // the actual child class
    auto modelSpecificPointer =
        std::static_pointer_cast<BSMPT::Models::Class_Potential_C2HDM>(
            modelPointer);

    // Define some abbreviations
    auto CosBeta = modelSpecificPointer->C_CosBeta;
    auto SinBeta = modelSpecificPointer->C_SinBeta;

    std::vector<double> alphas;
    alphas.push_back(modelSpecificPointer->alpha1);
    alphas.push_back(modelSpecificPointer->alpha2);
    alphas.push_back(modelSpecificPointer->alpha3);

    auto pos_h1 = modelSpecificPointer->pos_h1;
    auto pos_h2 = modelSpecificPointer->pos_h2;
    auto pos_h3 = modelSpecificPointer->pos_h3;

    auto pos_zeta1 = modelSpecificPointer->pos_zeta1;
    auto pos_zeta2 = modelSpecificPointer->pos_zeta2;
    auto pos_psi1  = modelSpecificPointer->pos_psi1;
    auto pos_psi2  = modelSpecificPointer->pos_psi2;

    const std::size_t n         = 3;
    const std::size_t numAngles = n * (n - 1) / 2;
    // Ensure that the test data is in the correct format, i.e. first all
    // rotation matrix elements, then n*(n - 1)/2 angles
    REQUIRE(rotation_matrix_angles.size() == n * n + numAngles);

    // Construct the neutral 3x3 mixing matrix
    std::vector<std::vector<double>> HiggsRotNeutralFixed(
        n, std::vector<double>(n));
    HiggsRotNeutralFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[0][2] =
        -modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                                pos_psi1) *
            SinBeta +
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_psi2) *
            CosBeta;

    HiggsRotNeutralFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[1][2] =
        -modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                                pos_psi1) *
            SinBeta +
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_psi2) *
            CosBeta;

    HiggsRotNeutralFixed[2][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[2][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[2][2] =
        -modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                                pos_psi1) *
            SinBeta +
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_psi2) *
            CosBeta;

    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        auto Rij_EnsuredConvention = HiggsRotNeutralFixed[i][j];
        auto expected              = rotation_matrix_angles[i * n + j];
        Check(Rij_EnsuredConvention, expected);
      }
    }

    // Check also the angles
    for (std::size_t i = 0; i < numAngles; ++i)
    {
      auto expected = rotation_matrix_angles[n * n + i];
      Check(alphas[i], expected);
    }

    // Do the checks also for the remaining elements (Goldstone, charged Higgs)
    auto pos_Gp = modelSpecificPointer->pos_Gp;
    auto pos_Gm = modelSpecificPointer->pos_Gm;
    auto pos_Hp = modelSpecificPointer->pos_Hp;
    auto pos_Hm = modelSpecificPointer->pos_Hm;
    auto pos_G0 = modelSpecificPointer->pos_G0;

    auto pos_rho1 = modelSpecificPointer->pos_rho1;
    auto pos_rho2 = modelSpecificPointer->pos_rho2;
    auto pos_eta1 = modelSpecificPointer->pos_eta1;
    auto pos_eta2 = modelSpecificPointer->pos_eta2;

    std::vector<std::vector<double>> HiggsRotFixed(2, std::vector<double>(2));
    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gp,
                                                               pos_rho1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gp,
                                                               pos_rho2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hp,
                                                               pos_rho1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hp,
                                                               pos_rho2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);

    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gm,
                                                               pos_eta1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gm,
                                                               pos_eta2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hm,
                                                               pos_eta1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hm,
                                                               pos_eta2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);

    auto G0Fixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
                       pos_G0, pos_psi1) *
                       CosBeta +
                   modelPointer->get_HiggsRotationMatrixEnsuredConvention(
                       pos_G0, pos_psi2) *
                       SinBeta;

    Check(G0Fixed, 1.);
  }
}

TEST_CASE("Checking sign convention of rotation matrix for CP in the Dark",
          "[signrotationmatrix]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CPINTHEDARK, SMConstants);
  auto nPar = modelPointer->get_nPar();

  std::ifstream testData(
      TEST_DATA_PATH
      "/signrotationmatrix_data/cpinthedark_unittest_signrotationmatrix.tsv");

  auto Check = [](auto result, auto expected)
  {
    if (std::abs(expected) > 1e-4)
    {
      REQUIRE(result == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  std::string line;
  std::getline(testData, line); // header, ignore

  while (std::getline(testData, line))
  {
    // Read in test data and parse properly: first the input parameters, then
    // the rotation matrix elements
    std::stringstream ss(line);
    std::vector<double> parameter_point;
    while (!ss.eof() && parameter_point.size() < nPar)
    {
      double x;
      ss >> x;
      parameter_point.push_back(x);
    }

    std::vector<double> rotation_matrix_angles;
    while (!ss.eof())
    {
      double x;
      ss >> x;
      rotation_matrix_angles.push_back(x);
    }

    modelPointer->initModel(parameter_point);

    // Since we need to access model-specific variables, cast modelPointer to
    // the actual child class
    auto modelSpecificPointer =
        std::static_pointer_cast<BSMPT::Models::Class_Potential_CPintheDark>(
            modelPointer);

    std::vector<double> alphas;
    alphas.push_back(modelSpecificPointer->alpha1);
    alphas.push_back(modelSpecificPointer->alpha2);
    alphas.push_back(modelSpecificPointer->alpha3);

    // Define some abbreviations
    auto pos_h1 = modelSpecificPointer->pos_h1;
    auto pos_h2 = modelSpecificPointer->pos_h2;
    auto pos_h3 = modelSpecificPointer->pos_h3;

    auto pos_zeta2 = modelSpecificPointer->pos_zeta2;
    auto pos_psi2  = modelSpecificPointer->pos_psi2;
    auto pos_rhoS  = modelSpecificPointer->pos_rhoS;

    const std::size_t n         = 3;
    const std::size_t numAngles = n * (n - 1) / 2;
    // Ensure that the test data is in the correct format, i.e. first all
    // rotation matrix elements, then n*(n - 1)/2 angles
    REQUIRE(rotation_matrix_angles.size() == n * n + numAngles);

    // Construct the neutral 3x3 mixing matrix
    std::vector<std::vector<double>> HiggsRotNeutralFixed(
        n, std::vector<double>(n));
    HiggsRotNeutralFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_psi2);
    HiggsRotNeutralFixed[0][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_rhoS);

    HiggsRotNeutralFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_psi2);
    HiggsRotNeutralFixed[1][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_rhoS);

    HiggsRotNeutralFixed[2][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[2][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_psi2);
    HiggsRotNeutralFixed[2][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_rhoS);

    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        auto Rij_EnsuredConvention = HiggsRotNeutralFixed[i][j];
        auto expected              = rotation_matrix_angles[i * n + j];
        Check(Rij_EnsuredConvention, expected);
      }
    }

    // Do not test the angles for CP in the Dark, as a1, a2, a3 from the data
    // input is potentially for a
    //   different mass ordering, as they are input values for ScannerS

    // Do the checks also for the remaining elements (Goldstone, charged Higgs)
    auto pos_Gp  = modelSpecificPointer->pos_Gp;
    auto pos_Gm  = modelSpecificPointer->pos_Gm;
    auto pos_Hp  = modelSpecificPointer->pos_Hp;
    auto pos_Hm  = modelSpecificPointer->pos_Hm;
    auto pos_HSM = modelSpecificPointer->pos_HSM;
    auto pos_G0  = modelSpecificPointer->pos_G0;

    auto pos_rho1  = modelSpecificPointer->pos_rho1;
    auto pos_eta1  = modelSpecificPointer->pos_eta1;
    auto pos_rho2  = modelSpecificPointer->pos_rho2;
    auto pos_eta2  = modelSpecificPointer->pos_eta2;
    auto pos_zeta1 = modelSpecificPointer->pos_zeta1;
    auto pos_psi1  = modelSpecificPointer->pos_psi1;

    auto GpFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_Gp, pos_rho1);
    auto GmFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_Gm, pos_eta1);
    auto HpFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_Hp, pos_rho2);
    auto HmFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_Hm, pos_eta2);
    auto HSMFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_HSM, pos_zeta1);
    auto G0Fixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_G0, pos_psi1);

    Check(GpFixed, 1.);
    Check(GmFixed, 1.);
    Check(HpFixed, 1.);
    Check(HmFixed, 1.);
    Check(HSMFixed, 1.);
    Check(G0Fixed, 1.);
  }
}

TEST_CASE("Checking sign convention of rotation matrix for CxSM",
          "[signrotationmatrix]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::CXSM, SMConstants);
  auto nPar = modelPointer->get_nPar();

  std::ifstream testData(
      TEST_DATA_PATH
      "/signrotationmatrix_data/cxsm_unittest_signrotationmatrix.tsv");

  auto Check = [](auto result, auto expected)
  {
    if (std::abs(expected) > 1e-4)
    {
      REQUIRE(result == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  std::string line;
  std::getline(testData, line); // header, ignore

  while (std::getline(testData, line))
  {
    // Read in test data and parse properly: first the input parameters, then
    // the rotation matrix elements
    std::stringstream ss(line);
    std::vector<double> parameter_point;
    while (!ss.eof() && parameter_point.size() < nPar)
    {
      double x;
      ss >> x;
      parameter_point.push_back(x);
    }

    std::vector<double> rotation_matrix_angles;
    while (!ss.eof())
    {
      double x;
      ss >> x;
      rotation_matrix_angles.push_back(x);
    }

    modelPointer->initModel(parameter_point);

    // Since we need to access model-specific variables, cast modelPointer to
    // the actual child class
    auto modelSpecificPointer =
        std::static_pointer_cast<BSMPT::Models::Class_CxSM>(modelPointer);

    std::vector<double> alphas;
    alphas.push_back(modelSpecificPointer->alpha1);
    alphas.push_back(modelSpecificPointer->alpha2);
    alphas.push_back(modelSpecificPointer->alpha3);

    // Define some abbreviations
    auto pos_h1 = modelSpecificPointer->pos_h1;
    auto pos_h2 = modelSpecificPointer->pos_h2;
    auto pos_h3 = modelSpecificPointer->pos_h3;

    auto pos_zeta1 = modelSpecificPointer->pos_zeta1;
    auto pos_zeta2 = modelSpecificPointer->pos_zeta2;
    auto pos_zeta3 = modelSpecificPointer->pos_zeta3;

    const std::size_t n         = 3;
    const std::size_t numAngles = n * (n - 1) / 2;
    // Ensure that the test data is in the correct format, i.e. first all
    // rotation matrix elements, then n*(n - 1)/2 angles
    REQUIRE(rotation_matrix_angles.size() == n * n + numAngles);

    // Construct the neutral 3x3 mixing matrix
    std::vector<std::vector<double>> HiggsRotNeutralFixed(
        n, std::vector<double>(n));
    HiggsRotNeutralFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[0][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta3);

    HiggsRotNeutralFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[1][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta3);

    HiggsRotNeutralFixed[2][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[2][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[2][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta3);

    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        auto Rij_EnsuredConvention = HiggsRotNeutralFixed[i][j];
        auto expected              = rotation_matrix_angles[i * n + j];
        Check(Rij_EnsuredConvention, expected);
      }
    }

    // Do not test the angles for CxSM, as a1, a2, a3 from the data input is
    // potentially for a
    //   different mass ordering, as they are input values for ScannerS

    // Do the checks also for the remaining elements (Goldstone, charged Higgs)
    auto pos_Gp = modelSpecificPointer->pos_Gp;
    auto pos_Gm = modelSpecificPointer->pos_Gm;
    auto pos_G0 = modelSpecificPointer->pos_G0;

    auto pos_i_Gp = modelSpecificPointer->pos_i_Gp;
    auto pos_i_Gm = modelSpecificPointer->pos_i_Gm;
    auto pos_i_G0 = modelSpecificPointer->pos_i_G0;

    auto GpFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_Gp, pos_i_Gp);
    auto GmFixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_Gm, pos_i_Gm);
    auto G0Fixed = modelPointer->get_HiggsRotationMatrixEnsuredConvention(
        pos_G0, pos_i_G0);

    Check(GpFixed, 1.);
    Check(GmFixed, 1.);
    Check(G0Fixed, 1.);
  }
}

TEST_CASE("Checking sign convention of rotation matrix for N2HDM",
          "[signrotationmatrix]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::N2HDM, SMConstants);
  auto nPar = modelPointer->get_nPar();

  std::ifstream testData(
      TEST_DATA_PATH
      "/signrotationmatrix_data/n2hdm_unittest_signrotationmatrix.tsv");

  auto Check = [](auto result, auto expected)
  {
    if (std::abs(expected) > 1e-4)
    {
      REQUIRE(result == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  std::string line;
  std::getline(testData, line); // header, ignore

  while (std::getline(testData, line))
  {
    // Read in test data and parse properly: first the input parameters, then
    // the rotation matrix elements
    std::stringstream ss(line);
    std::vector<double> parameter_point;
    while (!ss.eof() && parameter_point.size() < nPar)
    {
      double x;
      ss >> x;
      parameter_point.push_back(x);
    }

    std::vector<double> rotation_matrix_angles;
    while (!ss.eof())
    {
      double x;
      ss >> x;
      rotation_matrix_angles.push_back(x);
    }

    modelPointer->initModel(parameter_point);

    // Since we need to access model-specific variables, cast modelPointer to
    // the actual child class
    auto modelSpecificPointer =
        std::static_pointer_cast<BSMPT::Models::Class_Potential_N2HDM>(
            modelPointer);

    // Define some abbreviations
    auto CosBeta = modelSpecificPointer->C_CosBeta;
    auto SinBeta = modelSpecificPointer->C_SinBeta;

    std::vector<double> alphas;
    alphas.push_back(modelSpecificPointer->alpha1);
    alphas.push_back(modelSpecificPointer->alpha2);
    alphas.push_back(modelSpecificPointer->alpha3);

    auto pos_h1 = modelSpecificPointer->pos_h1;
    auto pos_h2 = modelSpecificPointer->pos_h2;
    auto pos_h3 = modelSpecificPointer->pos_h3;

    auto pos_zeta1 = modelSpecificPointer->pos_zeta1;
    auto pos_zeta2 = modelSpecificPointer->pos_zeta2;
    auto pos_rhoS  = modelSpecificPointer->pos_rhoS;

    const std::size_t n         = 3;
    const std::size_t numAngles = n * (n - 1) / 2;
    // Ensure that the test data is in the correct format, i.e. first all
    // rotation matrix elements, then n*(n - 1)/2 angles
    REQUIRE(rotation_matrix_angles.size() == n * n + numAngles);

    // Construct the neutral 3x3 mixing matrix
    std::vector<std::vector<double>> HiggsRotNeutralFixed(
        n, std::vector<double>(n));
    HiggsRotNeutralFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[0][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h1,
                                                               pos_rhoS);

    HiggsRotNeutralFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[1][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h2,
                                                               pos_rhoS);

    HiggsRotNeutralFixed[2][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[2][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_zeta2);
    HiggsRotNeutralFixed[2][2] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h3,
                                                               pos_rhoS);

    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        auto Rij_EnsuredConvention = HiggsRotNeutralFixed[i][j];
        auto expected              = rotation_matrix_angles[i * n + j];
        Check(Rij_EnsuredConvention, expected);
      }
    }

    // Check also the angles
    for (std::size_t i = 0; i < numAngles; ++i)
    {
      auto expected = rotation_matrix_angles[n * n + i];
      Check(alphas[i], expected);
    }

    // Do the checks also for the remaining elements (Goldstone, charged Higgs)
    auto pos_G0 = modelSpecificPointer->pos_G0;
    auto pos_A  = modelSpecificPointer->pos_A;
    auto pos_Gp = modelSpecificPointer->pos_Gp;
    auto pos_Gm = modelSpecificPointer->pos_Gm;
    auto pos_Hp = modelSpecificPointer->pos_Hp;
    auto pos_Hm = modelSpecificPointer->pos_Hm;

    auto pos_psi1 = modelSpecificPointer->pos_psi1;
    auto pos_psi2 = modelSpecificPointer->pos_psi2;
    auto pos_rho1 = modelSpecificPointer->pos_rho1;
    auto pos_eta1 = modelSpecificPointer->pos_eta1;
    auto pos_rho2 = modelSpecificPointer->pos_rho2;
    auto pos_eta2 = modelSpecificPointer->pos_eta2;

    std::vector<std::vector<double>> HiggsRotFixed(2, std::vector<double>(2));
    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_G0,
                                                               pos_psi1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_G0,
                                                               pos_psi2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_A, pos_psi1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_A, pos_psi2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);

    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gp,
                                                               pos_rho1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gp,
                                                               pos_rho2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hp,
                                                               pos_rho1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hp,
                                                               pos_rho2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);

    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gm,
                                                               pos_eta1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gm,
                                                               pos_eta2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hm,
                                                               pos_eta1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hm,
                                                               pos_eta2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);
  }
}

TEST_CASE("Checking sign convention of rotation matrix for R2HDM",
          "[signrotationmatrix]")
{
  using namespace BSMPT;
  const auto SMConstants = GetSMConstants();
  std::shared_ptr<BSMPT::Class_Potential_Origin> modelPointer =
      ModelID::FChoose(ModelID::ModelIDs::R2HDM, SMConstants);
  auto nPar = modelPointer->get_nPar();

  std::ifstream testData(
      TEST_DATA_PATH
      "/signrotationmatrix_data/r2hdm_unittest_signrotationmatrix.tsv");

  auto Check = [](auto result, auto expected)
  {
    if (std::abs(expected) > 1e-4)
    {
      REQUIRE(result == Approx(expected).epsilon(1e-4));
    }
    else
    {
      REQUIRE(std::abs(result) < 1e-4);
    }
  };

  std::string line;
  std::getline(testData, line); // header, ignore

  while (std::getline(testData, line))
  {
    // Read in test data and parse properly: first the input parameters, then
    // the rotation matrix elements
    std::stringstream ss(line);
    std::vector<double> parameter_point;
    while (!ss.eof() && parameter_point.size() < nPar)
    {
      double x;
      ss >> x;
      parameter_point.push_back(x);
    }

    // Only the mixing angle is given in the test data for the R2HDM
    double mixing_angle;
    ss >> mixing_angle;

    std::vector<double> rotation_matrix_angles;
    rotation_matrix_angles.push_back(std::cos(mixing_angle));
    rotation_matrix_angles.push_back(std::sin(mixing_angle));
    rotation_matrix_angles.push_back(-std::sin(mixing_angle));
    rotation_matrix_angles.push_back(std::cos(mixing_angle));
    rotation_matrix_angles.push_back(mixing_angle);

    modelPointer->initModel(parameter_point);

    // Since we need to access model-specific variables, cast modelPointer to
    // the actual child class
    auto modelSpecificPointer =
        std::static_pointer_cast<BSMPT::Models::Class_Potential_R2HDM>(
            modelPointer);

    // Define some abbreviations
    auto CosBeta = modelSpecificPointer->C_CosBeta;
    auto SinBeta = modelSpecificPointer->C_SinBeta;

    std::vector<double> alphas;
    alphas.push_back(modelSpecificPointer->alpha);

    auto pos_H = modelSpecificPointer->pos_H;
    auto pos_h = modelSpecificPointer->pos_h;

    auto pos_zeta1 = modelSpecificPointer->pos_zeta1;
    auto pos_zeta2 = modelSpecificPointer->pos_zeta2;

    const std::size_t n         = 2;
    const std::size_t numAngles = n * (n - 1) / 2;
    // Ensure that the test data is in the correct format, i.e. first all
    // rotation matrix elements, then n*(n - 1)/2 angles
    REQUIRE(rotation_matrix_angles.size() == n * n + numAngles);

    // Construct the neutral 2x2 mixing matrix
    std::vector<std::vector<double>> HiggsRotNeutralFixed(
        n, std::vector<double>(n));
    HiggsRotNeutralFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_H,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_H,
                                                               pos_zeta2);

    HiggsRotNeutralFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h,
                                                               pos_zeta1);
    HiggsRotNeutralFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_h,
                                                               pos_zeta2);

    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        auto Rij_EnsuredConvention = HiggsRotNeutralFixed[i][j];
        auto expected              = rotation_matrix_angles[i * n + j];
        Check(Rij_EnsuredConvention, expected);
      }
    }

    // Check also the angles
    for (std::size_t i = 0; i < numAngles; ++i)
    {
      auto expected = rotation_matrix_angles[n * n + i];
      Check(alphas[i], expected);
    }

    // Do the checks also for the remaining elements (Goldstone, charged Higgs)
    auto pos_Gp = modelSpecificPointer->pos_Gp;
    auto pos_Gm = modelSpecificPointer->pos_Gm;
    auto pos_Hp = modelSpecificPointer->pos_Hp;
    auto pos_Hm = modelSpecificPointer->pos_Hm;
    auto pos_G0 = modelSpecificPointer->pos_G0;
    auto pos_A  = modelSpecificPointer->pos_A;

    auto pos_rho1 = modelSpecificPointer->pos_rho1;
    auto pos_eta1 = modelSpecificPointer->pos_eta1;
    auto pos_rho2 = modelSpecificPointer->pos_rho2;
    auto pos_eta2 = modelSpecificPointer->pos_eta2;
    auto pos_psi1 = modelSpecificPointer->pos_psi1;
    auto pos_psi2 = modelSpecificPointer->pos_psi2;

    std::vector<std::vector<double>> HiggsRotFixed(2, std::vector<double>(2));
    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gp,
                                                               pos_rho1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gp,
                                                               pos_rho2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hp,
                                                               pos_rho1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hp,
                                                               pos_rho2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);

    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gm,
                                                               pos_eta1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Gm,
                                                               pos_eta2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hm,
                                                               pos_eta1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_Hm,
                                                               pos_eta2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);

    HiggsRotFixed[0][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_G0,
                                                               pos_psi1);
    HiggsRotFixed[0][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_G0,
                                                               pos_psi2);
    HiggsRotFixed[1][0] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_A, pos_psi1);
    HiggsRotFixed[1][1] =
        modelPointer->get_HiggsRotationMatrixEnsuredConvention(pos_A, pos_psi2);

    Check(HiggsRotFixed[0][0], CosBeta);
    Check(HiggsRotFixed[0][1], SinBeta);
    Check(HiggsRotFixed[1][0], -SinBeta);
    Check(HiggsRotFixed[1][1], CosBeta);
  }
}
