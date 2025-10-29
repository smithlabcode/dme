/* Copyright (C) 2025 Andrew D. Smith
 *
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "ScoringMatrix.hpp"

#include "Matrix.hpp"
#include "dme2_common.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <sstream>

float ScoringMatrix::default_correction = 0.0000000001;

std::string
ScoringMatrix::tostring() const {
  std::ostringstream s;
  if (!matrix.empty()) {
    s << "P0       A       C       G       T\n";
    for (std::size_t i = 0; i < width; i++) {
      s.width(2);
      s.fill('0');
      s << i;
      for (std::size_t j = 0; j < alphabet_size; j++) {
        s.width(8);
        s.fill(' ');
        s.setf(std::ios_base::right);
        s.precision(3);
        s << matrix[i][j];
      }
      s << std::endl;
    }
  }
  return s.str();
}

ScoringMatrix::ScoringMatrix(const Matrix &m,
                             const std::vector<float> &base_comp,
                             float correction) :
  matrix{m.matrix}, width{m.width} {
  if (correction == std::numeric_limits<float>::min())
    correction = default_correction;
  for (std::size_t i = 0; i < width; ++i) {
    const float tot = std::accumulate(std::cbegin(m[i]), std::cend(m[i]), 0.0f);
    for (std::size_t j = 0; j < alphabet_size; ++j)
      matrix[i][j] = std::log2(std::max(m.matrix[i][j], correction) / tot) -
                     std::log2(std::max(base_comp[j], correction));
  }
}

ScoringMatrix
ScoringMatrix::revcomp() const {
  ScoringMatrix sm(*this);
  for (auto &col : sm.matrix)
    std::reverse(std::begin(col), std::end(col));
  std::reverse(std::begin(sm.matrix), std::end(sm.matrix));
  return sm;
}
