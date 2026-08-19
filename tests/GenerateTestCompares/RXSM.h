// SPDX-FileCopyrightText: 2026 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include <BSMPT/minimizer/Minimizer.h>
#include <map>
#include <vector>
class Compare_RXSM
{
public:
  using Matrix3D = std::vector<std::vector<std::vector<double>>>;
  using Matrix2D = std::vector<std::vector<double>>;
  Compare_RXSM();
  Matrix3D CheckTripleCT;
  Matrix3D CheckTripleCW;
  Matrix3D CheckTripleTree;
  std::map<int, BSMPT::Minimizer::EWPTReturnType> EWPTPerSetting;
};
