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

#include "dme_zoops_workspace.hpp"
#include "CTSet.hpp"
#include "ScoringMatrix.hpp"
#include "dme2_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <ranges>  // IWYU pragma: keep

// NOLINTBEGIN (*-pointer-arithmetic)

// convert column types to corresponding scoring matrix columns
template <typename T>
[[nodiscard]] static auto
complement_scoremat(const std::vector<std::vector<T>> &coltypes,
                    std::vector<std::vector<T>> &scoremat) -> T {
  scoremat = coltypes;
  float max_score{};
  for (const auto &s : scoremat)
    // cppcheck-suppress useStlAlgorithm
    max_score = std::max(max_score, std::ranges::max(s));
  const auto cmpl = [&](const auto x) { return max_score - x; };
  for (auto &s : scoremat)
    std::ranges::transform(s, std::begin(s), cmpl);
  return max_score;
}

/// lexicographic tree stuff

struct dme_zoops_lextree {
  auto
  allocate_child_ptrs();
  auto
  insert(const std::string::const_iterator seq, const std::size_t depth,
         const std::size_t motif_width) -> void;
  ~dme_zoops_lextree();
  auto
  remove_matching(const std::vector<std::vector<float>> &score_matrix,
                  const std::size_t motif_width, const std::size_t depth,
                  const float current_score) -> bool;
  auto
  build(const std::string &sequence, const std::size_t motif_width);

  auto
  get_score(
    const std::vector<std::vector<std::vector<float>>> &refined_score_matrix,
    const std::vector<std::size_t> &path, std::size_t depth,
    float current_score, float &max_score) const -> void;

  auto
  get_score(const std::vector<std::vector<float>> &score_matrix,
            const std::vector<std::size_t> &path, std::size_t depth,
            float current_score, float &max_score) const -> void;

  dme_zoops_lextree **child{nullptr};
};

// allocate space for the array of children pointers of a dme_zoops_lextree
auto
dme_zoops_lextree::allocate_child_ptrs() {
  child = new dme_zoops_lextree *[alphabet_size];
  std::fill_n(child, alphabet_size,
              nullptr);  // static_cast<dme_zoops_lextree *>(0));
}

// insert a sequence below a subtree
auto
dme_zoops_lextree::insert(const std::string::const_iterator seq,
                          const std::size_t depth,
                          const std::size_t motif_width) -> void {
  if (motif_width <= depth)
    return;
  const auto index = encode_base[static_cast<char>(*seq)];
  if (child == nullptr)
    allocate_child_ptrs();
  if (child[index] == nullptr)
    child[index] = new dme_zoops_lextree;
  child[index]->insert(seq + 1, depth + 1, motif_width);
}

// build a dme_zoops_lextree from a sets of foreground and background sequences
auto
dme_zoops_lextree::build(const std::string &sequence,
                         const std::size_t motif_width) {
  allocate_child_ptrs();
  const std::size_t n_substrings = std::size(sequence) - motif_width + 1;
  for (std::size_t i = 0; i < n_substrings; ++i) {
    std::size_t j = 0;
    for (; j < motif_width && sequence[i + j] != 'N'; ++j)
      ;
    if (j == motif_width)
      insert(sequence.begin() + i, 0, motif_width);
  }
}

// recursively free space used by a lexicographic tree
dme_zoops_lextree::~dme_zoops_lextree() {
  if (child) {
    for (std::size_t i = 0; i < alphabet_size; ++i)
      if (child[i]) {
        delete child[i];
        child[i] = 0;
      }
    delete[] child;
    child = 0;
  }
}

bool
dme_zoops_lextree::remove_matching(
  const std::vector<std::vector<float>> &score_matrix,
  const std::size_t motif_width, const std::size_t depth,
  const float current_score) {
  std::size_t n_inactivated = 0;
  for (std::size_t i = 0; i < alphabet_size; ++i) {
    if (child[i]) {
      if (current_score - score_matrix[depth][i] > 0) {
        if (depth == motif_width - 1 ||
            child[i]->remove_matching(score_matrix, motif_width, depth + 1,
                                      current_score - score_matrix[depth][i])) {
          n_inactivated++;
          delete child[i];
          child[i] = 0;
        }
      }
    }
    else
      n_inactivated++;
  }
  // return true if no more children exist for the current node
  return (n_inactivated == alphabet_size);
}

/// erasing previous motifs

// remove all nodes in the tree for which all leaves below them correspond to
// strings that match the given matrix
auto
dme_zoops_workspace::deactivate(const ScoringMatrix &sm) -> void {
  std::vector<std::vector<float>> score_matrix(
    sm.width, std::vector<float>(alphabet_size));
  float deactivation_max_score = 0.0;
  for (std::size_t i = 0; i < sm.width; ++i) {
    float max_column_score = 0.0;
    for (std::size_t j = 0; j < alphabet_size; ++j) {
      score_matrix[i][j] = sm.matrix[i][j];
      max_column_score = std::max(score_matrix[i][j], max_column_score);
    }
    for (std::size_t j = 0; j < alphabet_size; j++)
      score_matrix[i][j] = max_column_score - score_matrix[i][j];
    deactivation_max_score += max_column_score;
  }
  for (std::size_t i = 0; i < fgsize; ++i)
    fgnodes.front()[i][0]->remove_matching(score_matrix, sm.width, 0,
                                           deactivation_max_score);
  for (std::size_t i = 0; i < bgsize; ++i)
    bgtrees[i]->remove_matching(score_matrix, sm.width, 0,
                                deactivation_max_score);
}

/// get background score

auto
dme_zoops_lextree::get_score(
  const std::vector<std::vector<std::vector<float>>> &refined_score_matrix,
  const std::vector<std::size_t> &path, std::size_t depth, float current_score,
  float &current_best) const -> void {
  if (child) {
    for (std::size_t i = 0; i < alphabet_size; ++i)
      if (child[i] &&
          current_score - refined_score_matrix[depth][path[depth]][i] >
            current_best)
        child[i]->get_score(refined_score_matrix, path, depth + 1,
                            current_score -
                              refined_score_matrix[depth][path[depth]][i],
                            current_best);
  }
  else
    current_best = std::max(current_best, current_score);
}

float
dme_zoops_workspace::get_background_refined_score() const {
  float total_score = 0;
  for (std::size_t i = 0; i < bgsize; ++i) {
    float current_best = 0;
    bgtrees[i]->get_score(refined_scoremat, prefix, 0, refined_max_score,
                          current_best);
    if (current_best > 0)
      total_score += current_best;
  }
  return total_score;
}

auto
dme_zoops_lextree::get_score(
  const std::vector<std::vector<float>> &score_matrix,
  const std::vector<std::size_t> &path, std::size_t depth, float current_score,
  float &current_best) const -> void {
  if (child) {
    for (std::size_t i = 0; i < alphabet_size; ++i)
      if (child[i] &&
          current_score - score_matrix[path[depth]][i] > current_best)
        child[i]->get_score(score_matrix, path, depth + 1,
                            current_score - score_matrix[path[depth]][i],
                            current_best);
  }
  else
    current_best = std::max(current_best, current_score);
}

float
dme_zoops_workspace::get_background_score() const {
  float total_score = 0;
  for (std::size_t i = 0; i < bgsize; ++i) {
    float current_best = 0;
    bgtrees[i]->get_score(scoremat, prefix, 0, max_score, current_best);
    if (current_best > 0)
      total_score += current_best;
  }
  return total_score;
}

////////////////////////////////////////////////////////////////////////
//////////////////////////                  ////////////////////////////
//////////////////////////    SEARCH STUFF  ////////////////////////////
//////////////////////////                  ////////////////////////////
////////////////////////////////////////////////////////////////////////

auto
dme_zoops_workspace::refined_enumeration(
  const std::size_t depth, const float surplus_information,
  const std::size_t remaining_changes) -> void {
  for (std::size_t i = 0;
       i == 0 || (remaining_changes > 0 && i < n_refined_types[depth - 1]); ++i)
    if (refined_coltype_bits[depth - 1][i] < surplus_information) {
      float upper_bound = 0.0;
      for (std::size_t j = 0; j < fgsize; ++j) {
        float current_best = 0.0;
        std::size_t frontier_size = 0;
        dme_zoops_lextree **lim = &fgnodes[depth - 1][j].front();
        for (std::size_t k = 0; lim[k] != 0; ++k) {
          const dme_zoops_lextree *n = fgnodes[depth - 1][j][k];
          for (std::size_t l = 0; l < alphabet_size; ++l)
            if (n->child[l]) {
              const float score =
                fgscore[depth - 1][j][k] - refined_scoremat[depth - 1][i][l];
              if (score > 0.0) {
                fgscore[depth][j][frontier_size] = score;
                fgnodes[depth][j][frontier_size] = n->child[l];
                current_best = std::max(current_best, score);
                frontier_size++;
              }
            }
        }
        fgnodes[depth][j][frontier_size] = 0;
        upper_bound += current_best;
      }
      if (upper_bound > best_score) {
        prefix[depth - 1] = i;
        if (depth == motif_width) {
          const float background_score = get_background_refined_score();
          if (upper_bound - background_score * lambda > best_score) {
            best_path = prefix;
            best_score = upper_bound - background_score * lambda;
          }
        }
        else
          refined_enumeration(
            depth + 1, surplus_information - refined_coltype_bits[depth - 1][i],
            remaining_changes - (i != 0));
      }
    }
}

DMEPath
dme_zoops_workspace::run_dme_zoops_local(
  const std::vector<std::vector<std::vector<float>>> &refined_types,
  const std::vector<std::vector<float>> &refined_bits,
  const float min_information, const std::size_t n_changes) {
  n_refined_types.resize(refined_types.size());
  for (std::size_t i = 0; i < refined_types.size(); ++i)
    n_refined_types[i] = refined_types[i].size();

  // get the surplus_information content matrix
  float max_col_type_info = 0.0;
  for (std::size_t i = 0; i < motif_width; ++i)
    for (std::size_t j = 0; j < n_refined_types[i]; j++)
      max_col_type_info = std::max(max_col_type_info, refined_bits[i][j]);

  refined_coltype_bits = refined_bits;
  for (std::size_t i = 0; i < motif_width; ++i)
    for (std::size_t j = 0; j < n_refined_types[i]; ++j)
      refined_coltype_bits[i][j] =
        max_col_type_info - refined_coltype_bits[i][j];

  // Set the minimum bits/column
  const float surplus_information =
    (max_col_type_info - min_information) * motif_width;

  // initialize the log scoring matrix and get the max score
  refined_scoremat.resize(motif_width);

  // float max_score = 0.0;
  refined_max_score = 0.0;
  for (std::size_t i = 0; i < motif_width; ++i)
    refined_max_score +=
      complement_scoremat(refined_types[i], refined_scoremat[i]);

  // initialize the max score for each fg and bg sequence
  for (std::size_t i = 0; i < fgsize; ++i)
    fgscore.front()[i][0] = refined_max_score;

  // make sure the variables holding the best current motif and score
  // are initialized empty and zero
  best_path.clear();
  best_score = 0.0;

  refined_enumeration(1, surplus_information, n_changes);

  return DMEPath(best_path, best_score);
}

auto
dme_zoops_workspace::enumeration(const std::size_t depth,
                                 const float surplus_information) -> void {
  for (std::size_t i = 0; i < n_types && coltype_bits[i] < surplus_information;
       ++i) {
    // get the foreground scores
    float upper_bound = 0.0;
    for (std::size_t j = 0; j < fgsize; ++j) {
      float current_best = 0.0;
      std::size_t frontier_size = 0;
      dme_zoops_lextree **lim = &fgnodes[depth - 1][j].front();
      for (std::size_t k = 0; lim[k] != 0; ++k) {
        const dme_zoops_lextree *n = lim[k];
        for (std::size_t l = 0; l < alphabet_size; ++l)
          if (n->child[l]) {
            const float score = fgscore[depth - 1][j][k] - scoremat[i][l];
            if (score > 0.0) {
              fgscore[depth][j][frontier_size] = score;
              fgnodes[depth][j][frontier_size] = n->child[l];
              current_best = std::max(current_best, score);
              frontier_size++;
            }
          }
      }
      fgnodes[depth][j][frontier_size] = 0;
      upper_bound += current_best;
    }
    if (upper_bound >= best_score) {
      prefix[depth - 1] = i;
      if (depth == motif_width) {
        const float background_score = get_background_score();
        if (upper_bound - background_score * lambda > best_score) {
          best_path = prefix;
          best_score = upper_bound - background_score * lambda;
        }
      }
      else
        enumeration(depth + 1, surplus_information - coltype_bits[i]);
    }
  }
}

DMEPath
dme_zoops_workspace::run_dme_zoops(
  const std::vector<std::vector<float>> &column_types,
  const std::vector<float> &bits, const float min_information) {
  // set the global alphabet and alphabet size variables
  n_types = column_types.size();

  // get the info content matrix
  coltype_bits = bits;

  // get the surplus_information content of each column type and the maximum
  // possible surplus_information in a column
  float max_col_type_info = 0.0;
  for (std::size_t i = 0; i < n_types; ++i)
    max_col_type_info = std::max(max_col_type_info, coltype_bits[i]);
  for (std::size_t i = 0; i < n_types; ++i)
    coltype_bits[i] = max_col_type_info - coltype_bits[i];

  // set a bound on how much the surplus_information can go under the maximum
  const float surplus_information =
    (max_col_type_info - min_information) * motif_width;

  // get scoring matrix
  // const float max_score = complement_scoremat(column_types, scoremat);
  max_score = complement_scoremat(column_types, scoremat) * motif_width;

  // set the values in the first cell in the scores arrays
  for (std::size_t i = 0; i < fgsize; ++i)
    fgscore.front()[i][0] = max_score;

  best_path.clear();
  best_score = 0.0;

  enumeration(1, surplus_information);

  return DMEPath(best_path, best_score);
}

dme_zoops_workspace::dme_zoops_workspace(
  const std::vector<std::string> &foreground,
  const std::vector<std::string> &background, const std::size_t width,
  const float adjustment) : motif_width{width} {
  // initialize the motif prefix
  prefix = std::vector<std::size_t>(motif_width, 0);

  fgsize = foreground.size();
  bgsize = background.size();

  lambda = (!background.empty()) ? static_cast<float>(fgsize) / bgsize : 1;

  lambda *= adjustment;

  // allocate the tables of nodes and scores
  fgscore = std::vector<std::vector<std::vector<float>>>(width + 1);
  fgnodes =
    std::vector<std::vector<std::vector<dme_zoops_lextree *>>>(width + 1);
  for (std::size_t i = 0; i <= width; ++i) {
    fgscore[i] = std::vector<std::vector<float>>(fgsize);
    fgnodes[i] = std::vector<std::vector<dme_zoops_lextree *>>(fgsize);
    for (std::size_t j = 0; j < fgsize; ++j) {
      const std::size_t max_frontier_size =
        std::min(std::size(foreground[j]),
                 static_cast<std::size_t>(std::pow(alphabet_size, i))) +
        1;
      fgscore[i][j] = std::vector<float>(max_frontier_size);
      fgnodes[i][j] = std::vector<dme_zoops_lextree *>(max_frontier_size);
    }
  }
  for (std::size_t i = 0; i < fgsize; ++i) {
    fgnodes.front()[i][0] = new dme_zoops_lextree;
    fgnodes.front()[i][0]->build(foreground[i], width);
    fgnodes.front()[i][1] = 0;
  }
  bgtrees = std::vector<dme_zoops_lextree *>(bgsize);
  for (std::size_t i = 0; i < bgsize; ++i) {
    bgtrees[i] = new dme_zoops_lextree;
    bgtrees[i]->build(background[i], width);
  }
}

dme_zoops_workspace::~dme_zoops_workspace() {
  for (std::size_t i = 0; i < fgnodes.front().size(); ++i)
    if (fgnodes.front()[i][0])
      delete fgnodes.front()[i][0];
  for (std::size_t i = 0; i < bgtrees.size(); ++i)
    if (bgtrees[i]) {
      delete bgtrees[i];
      bgtrees[i] = 0;
    }
}

// NOLINTEND (*-pointer-arithmetic)
