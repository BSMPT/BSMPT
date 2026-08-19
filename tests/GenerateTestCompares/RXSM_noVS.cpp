// SPDX-FileCopyrightText: 2026 Lisa Biermann, Margarete Mühlleitner, Rui
// Santos, João Viana
//
// SPDX-License-Identifier: GPL-3.0-or-later
#include "RXSM_noVS.h"
Compare_RXSM_noVS::Compare_RXSM_noVS()
{
  std::size_t NHiggs = 5;
  CheckTripleTree =
      Matrix3D{NHiggs, Matrix2D{NHiggs, std::vector<double>(NHiggs, 0)}};
  CheckTripleCW =
      Matrix3D{NHiggs, Matrix2D{NHiggs, std::vector<double>(NHiggs, 0)}};
  CheckTripleCT =
      Matrix3D{NHiggs, Matrix2D{NHiggs, std::vector<double>(NHiggs, 0)}};
  EWPTPerSetting[4].Tc = 159.073;
  EWPTPerSetting[4].vc = 22.4306;
  EWPTPerSetting[4].EWMinimum.push_back(-22.4306);
  EWPTPerSetting[4].EWMinimum.push_back(0);
  EWPTPerSetting[1].Tc = 159.073;
  EWPTPerSetting[1].vc = 22.4305;
  EWPTPerSetting[1].EWMinimum.push_back(22.4305);
  EWPTPerSetting[1].EWMinimum.push_back(0);
  EWPTPerSetting[5].Tc = 159.073;
  EWPTPerSetting[5].vc = 22.4306;
  EWPTPerSetting[5].EWMinimum.push_back(-22.4306);
  EWPTPerSetting[5].EWMinimum.push_back(0);
  EWPTPerSetting[2].Tc = 159.073;
  EWPTPerSetting[2].vc = 22.4303;
  EWPTPerSetting[2].EWMinimum.push_back(22.4303);
  EWPTPerSetting[2].EWMinimum.push_back(0);
  EWPTPerSetting[6].Tc = 159.073;
  EWPTPerSetting[6].vc = 22.4306;
  EWPTPerSetting[6].EWMinimum.push_back(-22.4306);
  EWPTPerSetting[6].EWMinimum.push_back(0);
  EWPTPerSetting[3].Tc = 159.073;
  EWPTPerSetting[3].vc = 22.4303;
  EWPTPerSetting[3].EWMinimum.push_back(22.4303);
  EWPTPerSetting[3].EWMinimum.push_back(0);
  EWPTPerSetting[7].Tc = 159.073;
  EWPTPerSetting[7].vc = 22.4306;
  EWPTPerSetting[7].EWMinimum.push_back(-22.4306);
  EWPTPerSetting[7].EWMinimum.push_back(0);
  CheckTripleTree.at(0).at(0).at(3) = 63.551;
  CheckTripleCT.at(0).at(0).at(3)   = -5.21724;
  CheckTripleCW.at(0).at(0).at(3)   = 5.21724;
  CheckTripleTree.at(0).at(3).at(0) = 63.551;
  CheckTripleCT.at(0).at(3).at(0)   = -5.21724;
  CheckTripleCW.at(0).at(3).at(0)   = 5.21724;
  CheckTripleTree.at(1).at(1).at(3) = 63.551;
  CheckTripleCT.at(1).at(1).at(3)   = -5.21724;
  CheckTripleCW.at(1).at(1).at(3)   = 5.21724;
  CheckTripleTree.at(1).at(3).at(1) = 63.551;
  CheckTripleCT.at(1).at(3).at(1)   = -5.21724;
  CheckTripleCW.at(1).at(3).at(1)   = 5.21724;
  CheckTripleTree.at(2).at(2).at(3) = 63.551;
  CheckTripleCT.at(2).at(2).at(3)   = -5.21724;
  CheckTripleCW.at(2).at(2).at(3)   = 5.21724;
  CheckTripleTree.at(2).at(3).at(2) = 63.551;
  CheckTripleCT.at(2).at(3).at(2)   = -5.21724;
  CheckTripleCW.at(2).at(3).at(2)   = 5.21724;
  CheckTripleTree.at(3).at(0).at(0) = 63.551;
  CheckTripleCT.at(3).at(0).at(0)   = -5.21724;
  CheckTripleCW.at(3).at(0).at(0)   = 5.21724;
  CheckTripleTree.at(3).at(1).at(1) = 63.551;
  CheckTripleCT.at(3).at(1).at(1)   = -5.21724;
  CheckTripleCW.at(3).at(1).at(1)   = 5.21724;
  CheckTripleTree.at(3).at(2).at(2) = 63.551;
  CheckTripleCT.at(3).at(2).at(2)   = -5.21724;
  CheckTripleCW.at(3).at(2).at(2)   = 5.21724;
  CheckTripleTree.at(3).at(3).at(3) = 190.653;
  CheckTripleCT.at(3).at(3).at(3)   = -15.6517;
  CheckTripleCW.at(3).at(3).at(3)   = -0.199475;
  CheckTripleTree.at(3).at(4).at(4) = -7.36761;
  CheckTripleCW.at(3).at(4).at(4)   = 0.0151512;
  CheckTripleTree.at(4).at(3).at(4) = -7.36761;
  CheckTripleCW.at(4).at(3).at(4)   = 0.0151512;
  CheckTripleTree.at(4).at(4).at(3) = -7.36761;
  CheckTripleCW.at(4).at(4).at(3)   = 0.0151512;
}
