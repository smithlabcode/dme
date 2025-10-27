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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "dme2_common.hpp"

#include "smithlab_utils.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

struct Matrix {
  Matrix() = default;
  Matrix(const std::vector<Column> &matrix, const std::uint32_t width) :
    matrix{matrix}, width{width} {};

  std::string
  tostring() const;

  [[nodiscard]] Column &
  operator[](const std::size_t i) {
    return matrix[i];
  }

  [[nodiscard]] const Column &
  operator[](const std::size_t i) const {
    return matrix[i];
  }

  Matrix
  revcomp() const;

  [[nodiscard]] float
  info(const std::vector<float> &) const;

  std::vector<Column> matrix;
  std::uint32_t width{};
};

inline std::ostream &
operator<<(std::ostream &o, const Matrix &m) {
  return o << m.tostring();
}

#endif
