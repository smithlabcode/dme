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

#ifndef CTSET_HPP
#define CTSET_HPP

#include "Matrix.hpp"
#include "dme2_common.hpp"

#include <array>
#include <cstddef>
#include <vector>

struct DMEPath {
  DMEPath(std::vector<size_t> p, float s) : path(p), score(s) {}
  std::vector<size_t> path;
  float score;
};

struct CTSet {
  CTSet(const std::vector<float> &col, const float newgran,
        const std::vector<float> &base_comp, const float correction = 1e-10);
  CTSet(const Column &col, const float newgran,
        const std::vector<float> &base_comp, const float correction = 1e-10);
  CTSet(const float granularity, const std::vector<float> &base_comp,
        const float correction = 1e-10);
  CTSet(const float granularity, const std::vector<float> &base_comp,
        const size_t n_sites);

  size_t
  size() const {
    return types.size();
  }

  [[nodiscard]] auto
  get_matrix() const -> const std::vector<std::vector<float>> & {
    return scoremat;
  }

  [[nodiscard]] auto
  get_bits() const -> const std::vector<float> & {
    return bits;
  }

  [[nodiscard]] auto
  path_to_matrix(const DMEPath &path) const -> Matrix;

  [[nodiscard]] static auto
  path_to_matrix(const DMEPath &path,
                 const std::vector<CTSet> &column_types) -> Matrix;

  [[nodiscard]] static auto
  get_bits(const std::vector<float> &col,
           const std::vector<float> &base_comp) -> float;

  std::vector<std::vector<float>> types;
  std::vector<std::vector<float>> scoremat;
  std::vector<float> bits;

  void
  sort_types(const std::vector<float> &base_comp);
  size_t
  generate_column(const size_t volume, const size_t max_volume,
                  const size_t depth, size_t index, std::vector<float> &v);

  static const size_t n_degen_nucs = 15;
  static std::array<std::array<float, 4>, 15> fixed_matrix;
  // float fixed_matrix[15][4];

  void
  build_scoremat(const std::vector<float> &base_comp, const float correction);
};

#endif
