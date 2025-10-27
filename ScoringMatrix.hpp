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

#ifndef SCORINGMATRIX_HPP
#define SCORINGMATRIX_HPP

#include "dme2_common.hpp"

#include <array>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

struct Matrix;

struct ScoringMatrix {
  ScoringMatrix() = default;
  ScoringMatrix(const Matrix &mat, const std::vector<float> &base_comp,
                float correction = std::numeric_limits<float>::min());

  [[nodiscard]] std::string
  tostring() const;

  [[nodiscard]] Column &
  operator[](const std::size_t i) {
    return matrix[i];
  }

  [[nodiscard]] const Column &
  operator[](const std::size_t i) const {
    return matrix[i];
  }

  ScoringMatrix
  revcomp() const;

  std::vector<Column> matrix;
  std::uint32_t width{};

  static float default_correction;
};

inline std::ostream &
operator<<(std::ostream &s, const ScoringMatrix &sm) {
  return s << sm.tostring();
}

#endif
