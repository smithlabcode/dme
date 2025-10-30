/* MIT License
 *
 * Copyright (c) 2025 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef SCORINGMATRIX_HPP
#define SCORINGMATRIX_HPP

#include "dme2_common.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

struct Matrix;

struct ScoringMatrix {
  static constexpr float default_correction = 0.0000000001;

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
};

inline std::ostream &
operator<<(std::ostream &s, const ScoringMatrix &sm) {
  return s << sm.tostring();
}

#endif
