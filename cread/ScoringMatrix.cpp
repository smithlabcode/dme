/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith, Pavel Sumazin and Michael Q. Zhang
 *
 * This file is part of CREAD.
 *
 * CREAD is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * CREAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CREAD; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "ScoringMatrix.hpp"

using std::string;
using std::accumulate;
using std::copy;
using std::vector;
using std::max;

string 
ScoringMatrix::tostring() const {
  std::ostringstream s;
  if (matrix != static_cast<float **>(0)) {
    s << "P0       A       C       G       T\n";
    for (size_t i = 0; i < width; i++) {
      s.width(2); 
      s.fill('0');
      s << i;
      for (size_t j = 0; j < alphabet_size; j++) {
	s.width(8); s.fill(' '); s.setf(std::ios_base::right); s.precision(3);
	s << matrix[i][j];
      }
      s << std::endl;
    }
  }
  return s.str();
}

void
ScoringMatrix::swap(ScoringMatrix &sm) {
  std::swap(matrix, sm.matrix);
  std::swap(width, sm.width);
}

ScoringMatrix&
ScoringMatrix::operator=(const ScoringMatrix& M) {
  ScoringMatrix tmp(M);
  swap(tmp);
  return *this;
}

ScoringMatrix::ScoringMatrix(float** m, const size_t w) :
  matrix(0), width(w) {
  if (m != 0) {
    matrix = new float *[width];
    for (size_t i = 0; i < width; i++) {
      matrix[i] = new float[alphabet_size];
      std::copy(m[i], m[i] + alphabet_size, matrix[i]);
    }
  }
}

ScoringMatrix::ScoringMatrix(const Matrix &m, 
			     const vector<float> &base_comp, 
			     float correction) : 
  matrix(new float *[m.get_width()]), 
  width(m.get_width()) {
  
  if (correction == std::numeric_limits<float>::min())
    correction = default_correction();
  for (size_t i = 0; i < width; ++i) {
    matrix[i] = new float[alphabet_size];
    const float count = accumulate(m[i], m[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; ++j)
      matrix[i][j] = cread::log2(max(m[i][j], correction)/count) - 
	cread::log2(max(base_comp[j], correction));
  }
}

ScoringMatrix::ScoringMatrix(const Matrix &m, const float base_comp[], float correction) :
  matrix(new float *[m.get_width()]), width(m.get_width()) {
  if (correction == std::numeric_limits<float>::min())
    correction = default_correction();
  for (size_t i = 0; i < width; i++) {
    matrix[i] = new float[alphabet_size];
    float count = accumulate(m[i], m[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; ++j)
      matrix[i][j] = cread::log2(max(m[i][j], correction)/count) - 
	cread::log2(max(base_comp[j], correction));
  }
}

ScoringMatrix::ScoringMatrix(const ScoringMatrix& mat) :
  matrix(0), width(mat.width) {
  if (mat.matrix != 0) {
    matrix = new float *[width];
    for (size_t i = 0; i < width; i++) {
      matrix[i] = new float[alphabet_size];
      std::copy(mat.matrix[i], mat.matrix[i] + alphabet_size, matrix[i]);
    }
  }
}

ScoringMatrix 
ScoringMatrix::revcomp() const {
  ScoringMatrix sm;
  sm.width = width;
  sm.matrix = new float*[width];
  for (size_t i = 0; i < width; i++) {
    sm.matrix[i] = new float[alphabet_size];
    std::reverse_copy(matrix[width - i - 1], 
		      matrix[width - i - 1] + alphabet_size, sm.matrix[i]);
  }
  return sm;
}

float
ScoringMatrix::functional_depth(float score) const {
  float min_score = 0, max_score = 0;
  for (size_t i = 0; i < width; ++i) {
    float temp_min_score = std::numeric_limits<float>::max();
    float temp_max_score = -std::numeric_limits<float>::max();
    for (size_t j = 0; j < alphabet_size; ++j) {
      temp_min_score = std::min(temp_min_score, matrix[i][j]);
      temp_max_score = std::max(temp_max_score, matrix[i][j]);
    }
    min_score += temp_min_score;
    max_score += temp_max_score;
  }
  return (score - min_score)/(max_score - min_score);
}

float
ScoringMatrix::functional_depth_to_score(const float fd) const {
  float min_score = 0, max_score = 0;
  for (size_t i = 0; i < width; ++i) {
    float temp_min_score = std::numeric_limits<float>::max();
    float temp_max_score = -std::numeric_limits<float>::max();
    for (size_t j = 0; j < alphabet_size; ++j) {
      temp_min_score = std::min(temp_min_score, matrix[i][j]);
      temp_max_score = std::max(temp_max_score, matrix[i][j]);
    }
    min_score += temp_min_score;
    max_score += temp_max_score;
  }
  return (max_score - min_score)*fd + min_score;
}

ScoringMatrix
ScoringMatrix::StormoScoringMatrix(const Matrix& M, const float *base_comp) {
  ScoringMatrix sm;
  sm.width = M.get_width();
  sm.matrix = new float *[sm.width];
  for (size_t i = 0; i < sm.width; i++) {
    sm.matrix[i] = new float[alphabet_size];
    float count = accumulate(M[i], M[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; j++) 
      sm.matrix[i][j] = cread::log2(M[i][j] + base_comp[j]) - 
	cread::log2(base_comp[j] * (count + 1));
  }
  return sm;
}

ScoringMatrix
ScoringMatrix::StormoScoringMatrix(const Matrix& M, 
				   const vector<float>& base_comp) {
  ScoringMatrix sm;
  sm.width = M.get_width();
  sm.matrix = new float *[sm.width];
  for (size_t i = 0; i < sm.width; i++) {
    sm.matrix[i] = new float[alphabet_size];
    float count = accumulate(M[i], M[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; j++) 
      sm.matrix[i][j] = cread::log2(M[i][j] + base_comp[j]) - 
	cread::log2(base_comp[j] * (count + 1));
  }
  return sm;
}

std::ostream& 
operator<<(std::ostream& s, const ScoringMatrix &sm) {
  return s << sm.tostring(); 
}

float
StadenPValue::normalize_matrix(ScoringMatrix& sm) const {
  
  using std::binder2nd;
  using std::bind2nd;
  using std::transform;
  typedef ScoringMatrix::iterator sm_itr;
  
  using std::multiplies;
  binder2nd<multiplies<float> > scaler(bind2nd(multiplies<float>(), scale));
  for (sm_itr i = sm.begin(); i != sm.end(); ++i)
    transform(*i, *i + alphabet_size, *i, scaler);
  
  float min_entry = std::numeric_limits<float>::max();
  for (sm_itr i = sm.begin(); i != sm.end(); ++i)
    min_entry = std::min(min_entry, *std::min_element(*i, *i + alphabet_size));
  
  using std::plus;
  binder2nd<plus<float> > translater(bind2nd(plus<float>(), -min_entry));
  for (sm_itr i = sm.begin(); i != sm.end(); ++i)
    transform(*i, *i + alphabet_size, *i, translater);
  
  for (sm_itr i = sm.begin(); i != sm.end(); ++i)
    transform(*i, *i + alphabet_size, *i, nearbyint);
  
  return -min_entry;
}

  using std::cerr;
  using std::endl;
  

float 
StadenPValue::get_pvalue(ScoringMatrix sm, const float score) const {
  // get the width once for efficiency
  const size_t width = sm.get_width();

  // do the required translation to make sm non-negative
  const float translation = normalize_matrix(sm);
  
  // determine the maximum score possible for sm
  float max_score = 0;
  for (size_t i = 0; i < width; ++i)
    max_score += *std::max_element(sm[i], sm[i] + alphabet_size);
  const size_t total_discretized_scores = static_cast<int>(max_score) + 1;
  
  // obtain an integer valued scoring matrix by rounding sm
  vector<vector<size_t> > smi(width, vector<size_t>(alphabet_size));
  for (size_t i = 0; i < width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      smi[i][j] = static_cast<size_t>(sm[i][j]);
  
  // make sure there is sufficient space in the working vectors
  if (total_discretized_scores > default_max_score + 1) {
    used.resize(total_discretized_scores);
    prev_used.resize(total_discretized_scores);
    counts.resize(total_discretized_scores);
    prev_counts.resize(total_discretized_scores);
  }
  
  // initialize the working vectors
  fill(used.begin(), used.end(), 0);
  fill(counts.begin(), counts.end(), 0.0);
  for (size_t j = 0; j < alphabet_size; ++j) {
    counts[smi[0][j]] += base_comp[j];
    used[smi[0][j]] = 1;
  }
  copy(used.begin(), used.end(), prev_used.begin());
  copy(counts.begin(), counts.end(), prev_counts.begin());
  
  // do the actual tabulation
  for (size_t i = 1; i < width; ++i) {
    vector<size_t>::const_iterator smi_iter = smi[i].begin();
    for (size_t j = 0; j < alphabet_size; ++j) {
      const float base_comp_multiplier = base_comp[j];
      for (size_t k = 0; k < total_discretized_scores; ++k)
	if (prev_used[k] == i) {
	  const size_t offset = k + smi_iter[j];
	  if (used[offset] < i + 1)
	    counts[offset] = 0;
	  counts[offset] += prev_counts[k]*base_comp_multiplier;
	  used[offset] = i + 1;
	}
    }
    used.swap(prev_used);
    counts.swap(prev_counts);
  }
  
  const size_t cutoff_score = static_cast<size_t>(score*scale + 
						  translation*width);
  double total = 0.0;
  for (size_t i = cutoff_score; i < total_discretized_scores; ++i)
    if (prev_used[i] == width)
      total += prev_counts[i];
  return total;
}

float
StadenPValue::get_score(ScoringMatrix sm, const float pvalue) const {

  // get the width once for efficiency
  const size_t width = sm.get_width();

  // do the required translation to make sm non-negative
  const float translation = normalize_matrix(sm);
  
  // determine the maximum score possible for sm
  float max_score = 0;
  for (size_t i = 0; i < width; ++i)
    max_score += *std::max_element(sm[i], sm[i] + alphabet_size);
  const size_t total_discretized_scores = static_cast<int>(max_score) + 1;
  
  // obtain an integer valued scoring matrix by rounding sm
  vector<vector<size_t> > smi(width, vector<size_t>(alphabet_size));
  for (size_t i = 0; i < width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      smi[i][j] = static_cast<size_t>(sm[i][j]);
  
  // make sure there is sufficient space in the working vectors
  if (total_discretized_scores > default_max_score + 1) {
    used.resize(total_discretized_scores);
    prev_used.resize(total_discretized_scores);
    counts.resize(total_discretized_scores);
    prev_counts.resize(total_discretized_scores);
  }
  
  // initialize the working vectors
  fill(used.begin(), used.end(), 0);
  fill(counts.begin(), counts.end(), 0.0);
  for (size_t j = 0; j < alphabet_size; ++j) {
    counts[smi[0][j]] += base_comp[j];
    used[smi[0][j]] = 1;
  }
  copy(used.begin(), used.end(), prev_used.begin());
  copy(counts.begin(), counts.end(), prev_counts.begin());

  // do the actual tabulation
  for (size_t i = 1; i < width; ++i) {
    vector<size_t>::const_iterator smi_iter = smi[i].begin();
    for (size_t j = 0; j < alphabet_size; ++j) {
      const float base_comp_multiplier = base_comp[j];
      for (size_t k = 0; k < total_discretized_scores; ++k)
	if (prev_used[k] == i) {
	  const size_t offset = k + smi_iter[j];
	  if (used[offset] < i + 1)
	    counts[offset] = 0;
	  counts[offset] += prev_counts[k]*base_comp_multiplier;
	  used[offset] = i + 1;
	}
    }
    used.swap(prev_used);
    counts.swap(prev_counts);
  }
  
  double total = 0.0;
  size_t cutoff_score = 0;
  for (size_t i = 0; i < prev_counts.size() && cutoff_score == 0; ++i)
    if (prev_used[total_discretized_scores - i - 1] == width) {
      total += prev_counts[total_discretized_scores - i - 1];
      if (total > pvalue)
	cutoff_score = total_discretized_scores - i;
    }
  return (cutoff_score - translation*width)/scale;
}

void
StadenPValue::get_pvalues(ScoringMatrix sm, 
			  const vector<float> &scores,
			  vector<float> &pvalues) const {
  // get the width once for efficiency
  const size_t width = sm.get_width();

  // do the required translation to make sm non-negative
  const float translation = normalize_matrix(sm);
  
  // determine the maximum score possible for sm
  float max_score = 0;
  for (size_t i = 0; i < width; ++i)
    max_score += *std::max_element(sm[i], sm[i] + alphabet_size);
  const size_t total_discretized_scores = static_cast<int>(max_score) + 1;
  
  // obtain an integer valued scoring matrix by rounding sm
  vector<vector<size_t> > smi(width, vector<size_t>(alphabet_size));
  for (size_t i = 0; i < width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      smi[i][j] = static_cast<size_t>(sm[i][j]);
  
  // make sure there is sufficient space in the working vectors
  if (total_discretized_scores > default_max_score + 1) {
    used.resize(total_discretized_scores);
    prev_used.resize(total_discretized_scores);
    counts.resize(total_discretized_scores);
    prev_counts.resize(total_discretized_scores);
  }
  
  // initialize the working vectors
  fill(used.begin(), used.end(), 0);
  fill(counts.begin(), counts.end(), 0.0);
  for (size_t j = 0; j < alphabet_size; ++j) {
    counts[smi[0][j]] += base_comp[j];
    used[smi[0][j]] = 1;
  }
  copy(used.begin(), used.end(), prev_used.begin());
  copy(counts.begin(), counts.end(), prev_counts.begin());
  
  // do the actual tabulation
  for (size_t i = 1; i < width; ++i) {
    vector<size_t>::const_iterator smi_iter = smi[i].begin();
    for (size_t j = 0; j < alphabet_size; ++j) {
      const float base_comp_multiplier = base_comp[j];
      for (size_t k = 0; k < total_discretized_scores; ++k)
	if (prev_used[k] == i) {
	  const size_t offset = k + smi_iter[j];
	  if (used[offset] < i + 1)
	    counts[offset] = 0;
	  counts[offset] += prev_counts[k]*base_comp_multiplier;
	  used[offset] = i + 1;
	}
    }
    used.swap(prev_used);
    counts.swap(prev_counts);
  }

  for (size_t i = 0; i < scores.size(); ++i) {
    const size_t cutoff_score = static_cast<size_t>(scores[i]*scale + 
						    translation*width);
    double total = 0.0;
    for (size_t j = cutoff_score; j < total_discretized_scores; ++j)
      if (prev_used[j] == width)
	total += prev_counts[j];
    pvalues.push_back(total);
  }
}

void
StadenPValue::get_scores(ScoringMatrix sm, 
			 const vector<float> &pvalues,
			 vector<float> &scores) const {
  
  // get the width once for efficiency
  const size_t width = sm.get_width();

  // do the required translation to make sm non-negative
  const float translation = normalize_matrix(sm);
  
  // determine the maximum score possible for sm
  float max_score = 0;
  for (size_t i = 0; i < width; ++i)
    max_score += *std::max_element(sm[i], sm[i] + alphabet_size);
  const size_t total_discretized_scores = static_cast<int>(max_score) + 1;
  
  // obtain an integer valued scoring matrix by rounding sm
  vector<vector<size_t> > smi(width, vector<size_t>(alphabet_size));
  for (size_t i = 0; i < width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      smi[i][j] = static_cast<size_t>(sm[i][j]);
  
  // make sure there is sufficient space in the working vectors
  if (total_discretized_scores > default_max_score + 1) {
    used.resize(total_discretized_scores);
    prev_used.resize(total_discretized_scores);
    counts.resize(total_discretized_scores);
    prev_counts.resize(total_discretized_scores);
  }
  
  // initialize the working vectors
  fill(used.begin(), used.end(), 0);
  fill(counts.begin(), counts.end(), 0.0);
  for (size_t j = 0; j < alphabet_size; ++j) {
    counts[smi[0][j]] += base_comp[j];
    used[smi[0][j]] = 1;
  }
  copy(used.begin(), used.end(), prev_used.begin());
  copy(counts.begin(), counts.end(), prev_counts.begin());

  // do the actual tabulation
  for (size_t i = 1; i < width; ++i) {
    vector<size_t>::const_iterator smi_iter = smi[i].begin();
    for (size_t j = 0; j < alphabet_size; ++j) {
      const float base_comp_multiplier = base_comp[j];
      for (size_t k = 0; k < total_discretized_scores; ++k)
	if (prev_used[k] == i) {
	  const size_t offset = k + smi_iter[j];
	  if (used[offset] < i + 1)
	    counts[offset] = 0;
	  counts[offset] += prev_counts[k]*base_comp_multiplier;
	  used[offset] = i + 1;
	}
    }
    used.swap(prev_used);
    counts.swap(prev_counts);
  }

  for (size_t i = 0; i < pvalues.size(); ++i) {
    double total = 0.0;
    size_t cutoff_score = 0;
    for (size_t j = 0; j < prev_counts.size() && cutoff_score == 0; ++j)
      if (prev_used[total_discretized_scores - j - 1] == width) {
	total += prev_counts[total_discretized_scores - j - 1];
	if (total > pvalues[i])
	  cutoff_score = total_discretized_scores - j;
      }
    scores.push_back((cutoff_score - translation*width)/scale);
  }
}
