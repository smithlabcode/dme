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

#include "CTSet.hpp"
#include "dme2_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>

// NOLINTBEGIN (*-c-arrays)

std::array<std::array<float, 4>, 15> CTSet::fixed_matrix = {{
  {1.000, 0.000, 0.000, 0.000},  // A
  {0.000, 1.000, 0.000, 0.000},  // C
  {0.000, 0.000, 1.000, 0.000},  // G
  {0.000, 0.000, 0.000, 1.000},  // T
  {0.500, 0.500, 0.000, 0.000},  // M
  {0.500, 0.000, 0.500, 0.000},  // R
  {0.500, 0.000, 0.000, 0.500},  // W
  {0.000, 0.500, 0.500, 0.000},  // S
  {0.000, 0.500, 0.000, 0.500},  // Y
  {0.000, 0.000, 0.500, 0.500},  // K
  {0.333, 0.333, 0.333, 0.000},  // V
  {0.333, 0.333, 0.000, 0.333},  // H
  {0.333, 0.000, 0.333, 0.333},  // D
  {0.000, 0.333, 0.333, 0.333},  // B
  {0.250, 0.250, 0.250, 0.250},  // N
}};

auto
CTSet::get_bits(const std::vector<float> &col,
                const std::vector<float> &base_comp) -> float {
  float column_bits = 0.0f;
  for (auto i = 0u; i < alphabet_size; ++i)
    column_bits += col[i] > 0.0f
                     ? col[i] * (std::log2(col[i]) - std::log2(base_comp[i]))
                     : 0.0f;
  return column_bits;
}

CTSet::CTSet(const float granularity, const std::vector<float> &base_comp,
             const float correction) {
  if (granularity == 0) {
    types = std::vector<std::vector<float>>(n_degen_nucs);
    for (size_t i = 0; i < n_degen_nucs; i++)
      std::copy(std::cbegin(fixed_matrix[i]), std::cend(fixed_matrix[i]),
                back_inserter(types[i]));
  }
  else {
    size_t volume = static_cast<size_t>(std::ceil(1.0 / granularity));
    size_t n_types = (volume + 1) * (volume + 2) * (volume + 3) / 6;
    types = std::vector<std::vector<float>>(n_types,
                                            std::vector<float>(alphabet_size));
    std::vector<float> v(alphabet_size);
    generate_column(volume, volume, 0, 0, v);
  }
  sort_types(base_comp);
  build_scoremat(base_comp, correction);
}

CTSet::CTSet(const std::vector<float> &original, const float newgran,
             const std::vector<float> &base_comp, const float correction) {
  // first do exact copy of original
  types.push_back(original);

  std::vector<float> temp(original);

  const auto orig_beg = std::cbegin(original);
  const auto orig_end = std::cend(original);
  const auto temp_beg = std::begin(temp);
  const auto temp_end = std::end(temp);

  // now do radius for all letters
  for (size_t i = 0; i < alphabet_size; ++i) {
    std::copy(orig_beg, orig_end, temp_beg);
    temp[i] += newgran;
    const float total = std::accumulate(temp_beg, temp_end, 0.0f);
    std::transform(temp_beg, temp_end, temp_beg,
                   [&](const auto x) { return x / total; });
    types.push_back(temp);
  }
  build_scoremat(base_comp, correction);
}

CTSet::CTSet(const Column &orig_col, const float newgran,
             const std::vector<float> &base_comp, const float correction) {
  std::vector<float> original(alphabet_size);
  std::copy(std::cbegin(orig_col), std::cend(orig_col), std::begin(original));

  // first do exact copy of original
  types.push_back(original);

  std::vector<float> temp(original);

  const auto orig_beg = std::cbegin(original);
  const auto orig_end = std::cend(original);
  const auto temp_beg = std::begin(temp);
  const auto temp_end = std::end(temp);

  // now do radius for all letters
  for (size_t i = 0; i < alphabet_size; ++i) {
    std::copy(orig_beg, orig_end, temp_beg);
    temp[i] += newgran;
    const float total = std::accumulate(temp_beg, temp_end, 0.0f);
    std::transform(temp_beg, temp_end, temp_beg,
                   [&](const auto x) { return x / total; });
    types.push_back(temp);
  }
  build_scoremat(base_comp, correction);
}

size_t
CTSet::generate_column(const size_t volume, const size_t max_volume,
                       const size_t depth, size_t index,
                       std::vector<float> &v) {
  for (size_t i = 0; i <= volume; ++i) {
    v[depth] = i;
    if (depth == 2) {
      v[depth + 1] = volume - i;
      for (size_t j = 0; j < alphabet_size; ++j)
        types[index][j] = v[j] / max_volume;
      index++;
    }
    else
      index = generate_column(volume - i, max_volume, depth + 1, index, v);
  }
  return index;
}

void
CTSet::sort_types(const std::vector<float> &base_comp) {
  std::vector<std::pair<float, size_t>> sorter;
  for (size_t i = 0; i < types.size(); ++i) {
    float column_bits = 0;
    for (size_t j = 0; j < alphabet_size; j++)
      if (types[i][j] > 0)
        column_bits +=
          types[i][j] * (std::log2(types[i][j]) - std::log2(base_comp[j]));
    sorter.emplace_back(column_bits, i);
  }
  std::sort(sorter.begin(), sorter.end(),
            std::greater<std::pair<float, size_t>>());
  std::vector<std::vector<float>> temp_types;
  for (size_t i = 0; i < types.size(); i++)
    temp_types.push_back(types[sorter[i].second]);
  types.clear();
  unique_copy(temp_types.begin(), temp_types.end(), back_inserter(types));
}

/* convert column types to corresponding scoring matrix columns */
void
CTSet::build_scoremat(const std::vector<float> &base_comp,
                      const float correction) {
  scoremat = types;
  bits.clear();
  for (size_t i = 0; i < types.size(); ++i) {
    float col_bits = 0;
    for (size_t j = 0; j < types[i].size(); j++) {
      scoremat[i][j] =
        ((types[i][j] > 0.0) ? std::log2(types[i][j]) : std::log2(correction)) -
        std::log2(base_comp[j]);
      if (types[i][j] > 0.0)
        col_bits +=
          types[i][j] * (std::log2(types[i][j]) - std::log2(base_comp[j]));
    }
    bits.push_back(col_bits);
  }
}

Matrix
CTSet::path_to_matrix(const DMEPath &path) const {
  const std::uint32_t n = std::size(path.path);
  auto best_matrix = std::vector<Column>(n);
  for (size_t i = 0; i < path.path.size(); ++i)
    std::copy(std::cbegin(types[path.path[i]]), std::cend(types[path.path[i]]),
              std::begin(best_matrix[i]));
  return Matrix(best_matrix, n);
}

Matrix
CTSet::path_to_matrix(const DMEPath &path,
                      const std::vector<CTSet> &column_types) {
  const std::uint32_t n = std::size(path.path);
  auto best_matrix = std::vector<Column>(n);
  for (size_t i = 0; i < path.path.size(); ++i)
    std::copy(std::cbegin(column_types[i].types[path.path[i]]),
              std::cend(column_types[i].types[path.path[i]]),
              std::begin(best_matrix[i]));
  return Matrix(best_matrix, n);
}

// NOLINTEND (*-c-arrays)
