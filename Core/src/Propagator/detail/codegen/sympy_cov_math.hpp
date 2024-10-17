// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Note: This file is generated by generate_sympy_cov.py
//       Do not modify it manually.

#pragma once

#include <cmath>

template <typename T>
void transportCovarianceToBoundImpl(const T* C, const T* J_full, T* new_C) {
  const auto x0 = C[0] * J_full[0] + C[5] * J_full[5] + C[10] * J_full[10] +
                  C[15] * J_full[15] + C[20] * J_full[20];
  const auto x1 = C[5] * J_full[0] + C[6] * J_full[5] + C[11] * J_full[10] +
                  C[16] * J_full[15] + C[21] * J_full[20];
  const auto x2 = C[10] * J_full[0] + C[11] * J_full[5] + C[12] * J_full[10] +
                  C[17] * J_full[15] + C[22] * J_full[20];
  const auto x3 = C[15] * J_full[0] + C[16] * J_full[5] + C[17] * J_full[10] +
                  C[18] * J_full[15] + C[23] * J_full[20];
  const auto x4 = C[20] * J_full[0] + C[21] * J_full[5] + C[22] * J_full[10] +
                  C[23] * J_full[15] + C[24] * J_full[20];
  const auto x5 = C[0] * J_full[1] + C[5] * J_full[6] + C[10] * J_full[11] +
                  C[15] * J_full[16] + C[20] * J_full[21];
  const auto x6 = C[5] * J_full[1] + C[6] * J_full[6] + C[11] * J_full[11] +
                  C[16] * J_full[16] + C[21] * J_full[21];
  const auto x7 = C[10] * J_full[1] + C[11] * J_full[6] + C[12] * J_full[11] +
                  C[17] * J_full[16] + C[22] * J_full[21];
  const auto x8 = C[15] * J_full[1] + C[16] * J_full[6] + C[17] * J_full[11] +
                  C[18] * J_full[16] + C[23] * J_full[21];
  const auto x9 = C[20] * J_full[1] + C[21] * J_full[6] + C[22] * J_full[11] +
                  C[23] * J_full[16] + C[24] * J_full[21];
  const auto x10 = C[0] * J_full[2] + C[5] * J_full[7] + C[10] * J_full[12] +
                   C[15] * J_full[17] + C[20] * J_full[22];
  const auto x11 = C[5] * J_full[2] + C[6] * J_full[7] + C[11] * J_full[12] +
                   C[16] * J_full[17] + C[21] * J_full[22];
  const auto x12 = C[10] * J_full[2] + C[11] * J_full[7] + C[12] * J_full[12] +
                   C[17] * J_full[17] + C[22] * J_full[22];
  const auto x13 = C[15] * J_full[2] + C[16] * J_full[7] + C[17] * J_full[12] +
                   C[18] * J_full[17] + C[23] * J_full[22];
  const auto x14 = C[20] * J_full[2] + C[21] * J_full[7] + C[22] * J_full[12] +
                   C[23] * J_full[17] + C[24] * J_full[22];
  const auto x15 = C[0] * J_full[3] + C[5] * J_full[8] + C[10] * J_full[13] +
                   C[15] * J_full[18] + C[20] * J_full[23];
  const auto x16 = C[5] * J_full[3] + C[6] * J_full[8] + C[11] * J_full[13] +
                   C[16] * J_full[18] + C[21] * J_full[23];
  const auto x17 = C[10] * J_full[3] + C[11] * J_full[8] + C[12] * J_full[13] +
                   C[17] * J_full[18] + C[22] * J_full[23];
  const auto x18 = C[15] * J_full[3] + C[16] * J_full[8] + C[17] * J_full[13] +
                   C[18] * J_full[18] + C[23] * J_full[23];
  const auto x19 = C[20] * J_full[3] + C[21] * J_full[8] + C[22] * J_full[13] +
                   C[23] * J_full[18] + C[24] * J_full[23];
  new_C[0] = x0 * J_full[0] + x1 * J_full[5] + x2 * J_full[10] +
             x3 * J_full[15] + x4 * J_full[20];
  new_C[1] = x5 * J_full[0] + x6 * J_full[5] + x7 * J_full[10] +
             x8 * J_full[15] + x9 * J_full[20];
  new_C[2] = x10 * J_full[0] + x11 * J_full[5] + x12 * J_full[10] +
             x13 * J_full[15] + x14 * J_full[20];
  new_C[3] = x15 * J_full[0] + x16 * J_full[5] + x17 * J_full[10] +
             x18 * J_full[15] + x19 * J_full[20];
  new_C[4] = x4;
  new_C[5] = x0 * J_full[1] + x1 * J_full[6] + x2 * J_full[11] +
             x3 * J_full[16] + x4 * J_full[21];
  new_C[6] = x5 * J_full[1] + x6 * J_full[6] + x7 * J_full[11] +
             x8 * J_full[16] + x9 * J_full[21];
  new_C[7] = x10 * J_full[1] + x11 * J_full[6] + x12 * J_full[11] +
             x13 * J_full[16] + x14 * J_full[21];
  new_C[8] = x15 * J_full[1] + x16 * J_full[6] + x17 * J_full[11] +
             x18 * J_full[16] + x19 * J_full[21];
  new_C[9] = x9;
  new_C[10] = x0 * J_full[2] + x1 * J_full[7] + x2 * J_full[12] +
              x3 * J_full[17] + x4 * J_full[22];
  new_C[11] = x5 * J_full[2] + x6 * J_full[7] + x7 * J_full[12] +
              x8 * J_full[17] + x9 * J_full[22];
  new_C[12] = x10 * J_full[2] + x11 * J_full[7] + x12 * J_full[12] +
              x13 * J_full[17] + x14 * J_full[22];
  new_C[13] = x15 * J_full[2] + x16 * J_full[7] + x17 * J_full[12] +
              x18 * J_full[17] + x19 * J_full[22];
  new_C[14] = x14;
  new_C[15] = x0 * J_full[3] + x1 * J_full[8] + x2 * J_full[13] +
              x3 * J_full[18] + x4 * J_full[23];
  new_C[16] = x5 * J_full[3] + x6 * J_full[8] + x7 * J_full[13] +
              x8 * J_full[18] + x9 * J_full[23];
  new_C[17] = x10 * J_full[3] + x11 * J_full[8] + x12 * J_full[13] +
              x13 * J_full[18] + x14 * J_full[23];
  new_C[18] = x15 * J_full[3] + x16 * J_full[8] + x17 * J_full[13] +
              x18 * J_full[18] + x19 * J_full[23];
  new_C[19] = x19;
  new_C[20] = x4;
  new_C[21] = x9;
  new_C[22] = x14;
  new_C[23] = x19;
  new_C[24] = C[24];
}
