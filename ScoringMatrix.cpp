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

#include "ScoringMatrix.hpp"
#include "Matrix.hpp"
#include "dme2_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <numeric>
#include <sstream>

// NOLINTBEGIN (*-avoid-magic-numbers)
// NOLINTBEGIN (*-array-index)

std::string
ScoringMatrix::tostring() const {
  std::ostringstream s;
  if (!matrix.empty()) {
    s << "P0       A       C       G       T\n";
    for (std::size_t i = 0; i < width; i++) {
      s.width(2);
      (void)s.fill('0');
      s << i;
      for (std::size_t j = 0; j < alphabet_size; j++) {
        s.width(8);
        (void)s.fill(' ');
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

// NOLINTEND (*-array-index)
// NOLINTEND (*-avoid-magic-numbers)
