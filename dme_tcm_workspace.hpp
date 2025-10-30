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

#ifndef DME_TCM_WORKSPACE_HPP_
#define DME_TCM_WORKSPACE_HPP_

#include "CTSet.hpp"

#include <cstddef>
#include <string>
#include <vector>

struct ScoringMatrix;
struct dme_tcm_lextree;

struct dme_tcm_workspace {
  dme_tcm_workspace(const std::vector<std::string> &foreground,
                    const std::vector<std::string> &background,
                    const int motif_width_param, const float adjustment);
  ~dme_tcm_workspace();

  // Pass in the column type sets as precomputed scoring matrix column
  // types. They will be "complemented" internally to speed the
  // computation. Similar for the information content.
  DMEPath
  run_dme_tcm(const std::vector<std::vector<float>> &column_types,
              const std::vector<float> &bits, const float min_information);

  // Pass in the column type sets as precomputed scoring matrix column
  // types. They will be "complemented" internally to speed the
  // computation. Similar for the information content.
  DMEPath
  run_dme_tcm_local(
    const std::vector<std::vector<std::vector<float>>> &refined_column_types,
    const std::vector<std::vector<float>> &refined_bits,
    const float min_information, const size_t n_changes);

  void
  deactivate(const ScoringMatrix &sm);

  void
  refined_enumeration(const size_t depth, const size_t prev_frontier,
                      const float surplus_information,
                      const size_t remaining_changes);
  void
  enumeration(const size_t depth, const size_t prev_frontier,
              const float surplus_information);

  // These are set initially and not changed
  std::vector<std::vector<float>> score;
  std::vector<std::vector<dme_tcm_lextree *>> nodes;
  dme_tcm_lextree *tree{nullptr};
  size_t motif_width;  // width of motif
  std::vector<size_t> prefix;

  // standard search variables
  std::vector<std::vector<float>> scoremat;
  std::vector<float> coltype_bits;
  size_t n_types{};

  // local search variables
  std::vector<std::vector<std::vector<float>>> refined_scoremat;
  std::vector<size_t> n_refined_types;
  std::vector<std::vector<float>> refined_coltype_bits;

  // for both standard and local search (reset right away at each
  // public function call)
  std::vector<size_t> best_path;
  float best_score{};
};

#endif  // DME_TCM_WORKSPACE_HPP_
