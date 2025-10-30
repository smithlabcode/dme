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

#include "Matrix.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iterator>
#include <numeric>

// NOLINTBEGIN (*-avoid-magic-numbers)
// NOLINTBEGIN (*-array-index)

std::string
Matrix::tostring() const {
  // std::ostringstream s;
  std::string s = std::format("P0       A       C       G       T\n");
  // s << "P0       A       C       G       T\n";
  for (std::size_t i = 0; i < width; ++i) {
    // prepend row number, width=2, zero-padded
    s += std::format("{:02}", i + 1);
    for (std::size_t j = 0; j < alphabet_size; ++j) {
      const float val = matrix[i][j];
      const bool has_fraction = (val - std::floor(val)) > 0.00001;
      if (has_fraction)
        s += std::format("{:8.3f}", val);  // width 8, fixed, 3 decimals
      else
        s += std::format("{:8.0f}", val);  // width 8, fixed, 0 decimals
    }
    if (i != width - 1)
      s += '\n';
  }
  //   for (std::size_t i = 0; i < width; i++) {
  //     s.width(2);
  //     (void)s.fill('0');
  //     s << i + 1;
  //     for (std::size_t j = 0; j < x; j++) {
  //       s.width(8);
  //       (void)s.fill(' ');
  //       s.setf(std::ios_base::right);
  //       s.setf(std::ios_base::fixed, std::ios_base::floatfield);
  //       if (matrix[i][j] - std::floor(matrix[i][j]) > 0.00001)
  //         s.precision(3);
  //       else
  //         s.precision(0);
  //       s << matrix[i][j];
  //     }
  //     if (i != width - 1)
  //       s << '\n';
  //   }
  // }
  return s;  // s.str();
}

float
Matrix::info(const std::vector<float> &f) const {
  float information = 0;
  for (std::size_t i = 0; i < width; i++) {
    const float sum =
      std::accumulate(std::cbegin(matrix[i]), std::cend(matrix[i]), 0.0f);
    for (std::size_t j = 0; j < alphabet_size; j++)
      if (matrix[i][j] > 0)
        information += matrix[i][j] / sum *
                       (std::log2(matrix[i][j] / sum) - std::log2(f[j]));
  }
  return information;
}

Matrix
Matrix::revcomp() const {
  Matrix m = *this;
  for (auto &c : m.matrix)
    std::reverse(std::begin(c), std::end(c));
  std::reverse(std::begin(m.matrix), std::end(m.matrix));
  return m;
}

// NOLINTEND (*-array-index)
// NOLINTEND (*-avoid-magic-numbers)
