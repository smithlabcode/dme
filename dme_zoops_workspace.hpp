/* MIT License
 *
 * Copyright (c) 2024 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef DME_ZOOPS_WORKSPACE_HPP_
#define DME_ZOOPS_WORKSPACE_HPP_

#include "CTSet.hpp"

#include <cstddef>
#include <string>
#include <vector>

struct ScoringMatrix;
struct dme_zoops_lextree;

struct dme_zoops_workspace {
  dme_zoops_workspace(const std::vector<std::string> &foreground,
                      const std::vector<std::string> &background,
                      const std::size_t motif_width_param,
                      const float adjustment);
  ~dme_zoops_workspace();

  [[nodiscard]] auto
  run_dme_zoops(const std::vector<std::vector<float>> &column_types,
                const std::vector<float> &bits,
                const float min_information) -> DMEPath;

  [[nodiscard]] auto
  run_dme_zoops_local(
    const std::vector<std::vector<std::vector<float>>> &refined_column_types,
    const std::vector<std::vector<float>> &refined_bits,
    const float min_information, const std::size_t n_changes) -> DMEPath;

  auto
  deactivate(const ScoringMatrix &sm) -> void;

  auto
  refined_enumeration(const std::size_t depth, const float surplus_information,
                      const std::size_t remaining_changes) -> void;
  auto
  enumeration(const std::size_t depth, const float surplus_information) -> void;

  float
  get_background_score() const;
  float
  get_background_refined_score() const;

  // These are set initially and not changed
  std::vector<std::vector<std::vector<float>>> fgscore;
  std::vector<std::vector<std::vector<dme_zoops_lextree *>>> fgnodes;
  std::vector<dme_zoops_lextree *> bgtrees;

  std::size_t fgsize{};
  std::size_t bgsize{};

  std::size_t motif_width{};
  std::vector<std::size_t> prefix;

  // standard search variables
  std::vector<std::vector<float>> scoremat;
  std::vector<float> coltype_bits;
  std::size_t n_types{};
  float max_score{};

  // local search variables
  std::vector<std::vector<std::vector<float>>> refined_scoremat;
  std::vector<std::size_t> n_refined_types;
  std::vector<std::vector<float>> refined_coltype_bits;
  float refined_max_score{};

  // for both standard and local search (reset right away at each
  // public function call)
  std::vector<std::size_t> best_path;
  float best_score{};
  float lambda{};
};

#endif  // DME_ZOOPS_WORKSPACE_HPP_
