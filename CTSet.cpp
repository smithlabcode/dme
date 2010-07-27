/*
 * Copyright (C) 2008 Cold Spring Harbor Laboratory and Andrew D Smith
 * Author: Andrew D Smith
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "CTSet.hpp"

using std::vector;
using std::string;
using std::pair;
using std::sort;
using std::max;

float 
CTSet::fixed_matrix[15][4] = {
  { 1.000, 0.000, 0.000, 0.000 },  // A
  { 0.000, 1.000, 0.000, 0.000 },  // C
  { 0.000, 0.000, 1.000, 0.000 },  // G
  { 0.000, 0.000, 0.000, 1.000 },  // T
  { 0.500, 0.500, 0.000, 0.000 },  // M
  { 0.500, 0.000, 0.500, 0.000 },  // R
  { 0.500, 0.000, 0.000, 0.500 },  // W
  { 0.000, 0.500, 0.500, 0.000 },  // S
  { 0.000, 0.500, 0.000, 0.500 },  // Y
  { 0.000, 0.000, 0.500, 0.500 },  // K
  { 0.333, 0.333, 0.333, 0.000 },  // V
  { 0.333, 0.333, 0.000, 0.333 },  // H
  { 0.333, 0.000, 0.333, 0.333 },  // D
  { 0.000, 0.333, 0.333, 0.333 },  // B
  { 0.250, 0.250, 0.250, 0.250 }   // N
};



float
CTSet::get_bits(const vector<float> &col,
		const vector<float> &base_comp) {
  float bits = 0;
  for (size_t i = 0; i < alphabet_size; ++i)
    if (col[i] > 0)
      bits += col[i]*(cread::log2(col[i]) - 
		      cread::log2(base_comp[i]));
  return bits;
}



CTSet::CTSet(const float granularity,
	     const vector<float> &base_comp,
	     const float correction) {
  if (granularity == 0) {
    types = vector<vector<float> >(n_degen_nucs);
    for (size_t i = 0; i < n_degen_nucs; i++)
      copy(fixed_matrix[i], fixed_matrix[i] + alphabet_size,
	   back_inserter(types[i]));
  }
  else {
    size_t volume = static_cast<size_t>(std::ceil(1.0/granularity));
    size_t n_types = (volume + 1)*(volume + 2)*(volume + 3)/6;
    types = vector<vector<float> >(n_types, vector<float>(alphabet_size));
    vector<float> v(alphabet_size);
    generate_column(volume, volume, 0, 0, v);
  }
  sort_types(base_comp);
  build_scoremat(base_comp, correction);
}


CTSet::CTSet(const float *original, const float newgran,
	     const vector<float> &base_comp, 
	     const float correction) {
  vector<float> temp(alphabet_size);
  copy(original, original + alphabet_size, temp.begin());
  types.push_back(temp);
  for (size_t i = 0; i < alphabet_size; ++i) {
    copy(original, original + alphabet_size, temp.begin());
    temp[i] += newgran;
    const float total = accumulate(temp.begin(), temp.end(), 0.0);
    std::transform(temp.begin(), temp.end(), temp.begin(), 
		   std::bind2nd(std::divides<float>(), total));
    types.push_back(temp);
  }
  build_scoremat(base_comp, correction);
}



size_t
CTSet::generate_column(const size_t volume, const size_t max_volume,
		       const size_t depth, size_t index, vector<float>& v) {
  for (size_t i = 0; i <= volume; ++i) {
    v[depth] = i;
    if (depth == 2) {
      v[depth + 1] = volume - i;
      for (size_t j = 0; j < alphabet_size; ++j) 
	types[index][j] = v[j]/max_volume;
      index++;
    }
    else index = generate_column(volume - i, max_volume, depth + 1, index, v);
  }
  return index;
}



void 
CTSet::sort_types(const vector<float> &base_comp) {
  vector<pair<float, size_t> > sorter;
  for (size_t i = 0; i < types.size(); ++i) {
    float bits = 0;
    for (size_t j = 0; j < alphabet_size; j++)
      if (types[i][j] > 0)
	bits += types[i][j]*(log(types[i][j]) - log(base_comp[j]));
    sorter.push_back(std::make_pair(bits, i));
  }
  sort(sorter.begin(), sorter.end(), std::greater<pair<float, size_t> >());
  vector<vector<float> > temp_types;
  for (size_t i = 0; i < types.size(); i++)
    temp_types.push_back(types[sorter[i].second]);
  types.clear();
  unique_copy(temp_types.begin(), temp_types.end(), back_inserter(types));
}



/* convert column types to corresponding scoring matrix columns */
void
CTSet::build_scoremat(const vector<float> &base_comp, const float correction) {
  
  scoremat = types;
  bits.clear();
  for (size_t i = 0; i < types.size(); ++i) {
    float col_bits = 0;
    for (size_t j = 0; j < types[i].size(); j++) {
      scoremat[i][j] = ((types[i][j] > 0.0) ?
			log2(types[i][j]) : log2(correction)) - log2(base_comp[j]);
      if (types[i][j] > 0.0)
	col_bits += types[i][j]*(log2(types[i][j]) - log2(base_comp[j]));
    }
    bits.push_back(col_bits);
  }
}



Matrix
CTSet::path_to_matrix(const DMEPath &path) const {
  float *best_matrix[path.path.size()];
  for (size_t i = 0; i < path.path.size(); ++i) {
    best_matrix[i] = new float[alphabet_size];
    copy(types[path.path[i]].begin(), types[path.path[i]].end(), 
	 best_matrix[i]);
  }
  const Matrix matrix(best_matrix, path.path.size());
  for (size_t i = 0; i < path.path.size(); ++i)
    delete[] best_matrix[i];
  return matrix;
}



Matrix
CTSet::path_to_matrix(const DMEPath &path,
		      const vector<CTSet> &column_types) {
  float *best_matrix[path.path.size()];
  for (size_t i = 0; i < path.path.size(); ++i) {
    best_matrix[i] = new float[alphabet_size];
    copy(column_types[i].types[path.path[i]].begin(),
	 column_types[i].types[path.path[i]].end(),
	 best_matrix[i]);
  }
  const Matrix matrix(best_matrix, path.path.size());
  for (size_t i = 0; i < path.path.size(); ++i)
    delete[] best_matrix[i];
  return matrix;
}
