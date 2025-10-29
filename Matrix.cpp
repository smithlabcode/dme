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

#include "Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <numeric>

std::string
Matrix::tostring() const {
  std::ostringstream s;
  if (!matrix.empty()) {
    s << "P0       A       C       G       T\n";
    for (std::size_t i = 0; i < width; i++) {
      s.width(2);
      s.fill('0');
      s << i + 1;
      for (std::size_t j = 0; j < alphabet_size; j++) {
        s.width(8);
        s.fill(' ');
        s.setf(std::ios_base::right);
        s.setf(std::ios_base::fixed, std::ios_base::floatfield);
        if (matrix[i][j] - std::floor(matrix[i][j]) > 0.00001)
          s.precision(3);
        else
          s.precision(0);
        s << matrix[i][j];
      }
      if (i != width - 1)
        s << '\n';
    }
  }
  return s.str();
}

float
Matrix::info(const std::vector<float> &f) const {
  float information = 0;
  for (std::size_t i = 0; i < width; i++) {
    const float sum =
      std::accumulate(std::cbegin(matrix[i]), std::cend(matrix[i]), 0.0f);
    for (std::size_t j = 0; j < alphabet_size; j++)
      if (matrix[i][j] > 0)
        information +=
          matrix[i][j] / sum * (log2(matrix[i][j] / sum) - log2(f[j]));
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
