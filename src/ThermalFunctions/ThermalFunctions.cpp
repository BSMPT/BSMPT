// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/ThermalFunctions/thermalcoefficientcalculator.h>
#include <BSMPT/models/SMparam.h>
#include <complex>
#include <map>

#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

namespace BSMPT
{
namespace ThermalFunctions
{

namespace
{
ThermalCoefficientCalculator FermionInterpolatedLowCoefficientCalculator(
    [](int l) -> double
    {
      if (l < 2)
      {
        return 0;
      }
      else
      {
        return gsl_sf_doublefact(2 * l - 3) * gsl_sf_zeta(2 * l - 1) /
               (gsl_sf_doublefact(2 * l) * (l + 1)) * (pow(2, 2 * l - 1) - 1);
      }
    },
    5);

ThermalCoefficientCalculator BosonInterpolatedLowCoefficientCalculator(
    [](int l) -> double
    {
      if (l < 2)
      {
        return 0;
      }
      else
      {
        return gsl_sf_doublefact(2 * l - 3) * gsl_sf_zeta(2 * l - 1) /
               (gsl_sf_doublefact(2 * l) * (l + 1));
      }
    },
    5);

ThermalCoefficientCalculator JInterpolatedHighCoefficientCalculator(
    [](int l) -> double
    {
      return 1 / (std::pow(2, l) * gsl_sf_fact(l)) * gsl_sf_gamma(2.5 + l) /
             gsl_sf_gamma(2.5 - l);
    },
    5);
} // namespace

double JfermionInterpolatedLow(const double &x, const int &n, int diff)
{
  if (x == 0 and diff == 0)
  {
    return -7 * pow(M_PI, 4) / 360.0;
  }
  else if (x == 0 and diff == 1)
  {
    return pow(M_PI, 2) / 24;
  }
  using std::log;
  using std::pow;
  double res = 0;
  double cf  = 1.5 + 2 * log(4 * M_PI) - 2 * C_euler_gamma - 2 * log(4);
  if (diff == 0)
  {
    res = -7 * pow(M_PI, 4) / 360.0;
    res += pow(M_PI, 2) / 24 * x;
    res += 1 / 32.0 * pow(x, 2) * (log(x) - cf);
    double sum = 0;

    for (int l = 2; l <= n; l++)
    {
      double Kl =
          FermionInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += pow(-x / (4 * pow(M_PI, 2)), l) * Kl;
    }
    res += -pow(M_PI, 2) * x * sum;
  }
  else if (diff == 1)
  {
    res = pow(M_PI, 2) / 24.0;
    res += x * (-6 * cf + 3) / 96;
    res += x * log(x) / 16;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      double Kl =
          FermionInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += -Kl * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
    }
    res += sum;
  }

  return res;
}

double JbosonInterpolatedLow(const double &x, const int &n, int diff)
{
  if (x == 0 and diff == 0)
  {
    return -pow(M_PI, 4) / 45.0;
  }
  else if (x == 0 and diff == 1)
  {
    return pow(M_PI, 2) / 12.0;
  }
  using std::log;
  using std::pow;
  using std::sqrt;
  double cb  = 1.5 + 2 * std::log(4 * M_PI) - 2 * C_euler_gamma;
  double res = 0;
  if (diff == 0)
  {
    res = -pow(M_PI, 4) / 45.0;
    res += pow(M_PI, 2) * x / 12.0;
    res += -M_PI * pow(x, 1.5) / 6;
    res += -pow(x, 2) * (log(x) - cb) / 32.0;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      double Kl =
          BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += pow(-x / (4 * pow(M_PI, 2)), l) * Kl;
    }
    res += pow(M_PI, 2) * x * sum;
  }
  else if (diff == 1)
  {
    res = pow(M_PI, 2) / 12.0;
    res += x * (6 * cb - 3) / 96.0;
    res += -x * log(x) / 16.0;
    res += -M_PI * sqrt(x) / 4.0;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      double Kl =
          BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += Kl * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
    }
    res += sum;
  }
  return res;
}

double JInterpolatedHigh(const double &x, const int &n, int diff)
{
  using std::exp;
  using std::pow;
  using std::sqrt;

  double res = 0;
  if (diff == 0)
  {
    double sum = 0;
    for (int l = 0; l <= n; l++)
    {
      double Kl =
          JInterpolatedHighCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += Kl * pow(x, -l / 2.0);
    }
    res = -exp(-sqrt(x)) * sqrt(M_PI / 2 * pow(x, 1.5)) * sum;
  }
  else if (diff == 1)
  {
    double sum = 0;
    for (int l = 0; l <= n; l++)
    {
      double Kl =
          JInterpolatedHighCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += Kl * pow(x, (1 - l) / 2) * (2 * l + 2 * sqrt(x) - 3);
    }
    res = exp(-sqrt(x)) * sqrt(2 * M_PI) / (8 * pow(x, 3.0 / 4.0)) * sum;
  }
  return res;
}

double JfermionInterpolated(const double &x, int diff)
{
  double res = 0;
  if (x >= C_FermionTheta)
  {
    res = -JInterpolatedHigh(x, 3, diff);
  }
  else
  {
    res = -JfermionInterpolatedLow(x, 4, diff);
    if (diff == 0) res += -C_FermionShift;
  }
  return res;
}

double JbosonInterpolated(const double &x, int diff)
{
  double res = 0;
  if (x >= C_BosonTheta)
  {
    res = JInterpolatedHigh(x, 3, diff);
    if (diff == 0) res -= C_BosonShift;
  }
  else if (x >= 0)
  {
    res = JbosonInterpolatedLow(x, 3, diff);
  }
  else if (x < 0)
  {
    res = JbosonInterpolatedNegative(x, diff);
  }
  return res;
}

double JbosonInterpolatedNegative(const double &x, int diff)
{
  if (x >= 0) return 0;
  double PotVal = 0;

  if (diff == 0)
  {
    PotVal = JbosonNegativeSpline(-x) + 3.533375127e-06;
  }
  else if (diff == 1)
  {
    PotVal = -JbosonNegativeSpline.deriv(1, -x);
  }

  return PotVal;
}

} // namespace ThermalFunctions
} // namespace BSMPT
