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

#ifndef CTSET_HPP
#define CTSET_HPP

// Warning: magic number in this file assumes alphabet_size==4

#include <cread.hpp>
#include <Matrix.hpp>

struct DMEPath {
  DMEPath(std::vector<size_t> p, float s) : path(p), score(s) {}
  std::vector<size_t> path;
  float score;
};

class CTSet {
public:
  CTSet(const float *, float, const std::vector<float> &base_comp,
	const float correction = 1e-10);
  CTSet(const float granularity, const std::vector<float> &base_comp,
	const float correction = 1e-10);
  CTSet(const float granularity, const std::vector<float> &base_comp,
	const size_t n_sites);
  
  size_t size() const {return types.size();}
  
  std::vector<std::vector<float> > get_matrix() const {return scoremat;}
  std::vector<float> get_bits() const {return bits;}
  
  Matrix path_to_matrix(const DMEPath &path) const;
  
  static Matrix
  path_to_matrix(const DMEPath &path,
		 const std::vector<CTSet> &column_types);
  
  static float
  get_bits(const std::vector<float> &col, 
	   const std::vector<float> &base_comp);
  
private:
  std::vector<std::vector<float> > types;
  std::vector<std::vector<float> > scoremat;
  std::vector<float> bits;
  
  void sort_types(const std::vector<float> &base_comp);
  size_t generate_column(const size_t volume, const size_t max_volume, 
			 const size_t depth, size_t index, std::vector<float>& v);
  
  static const size_t n_degen_nucs = 15;
  static float fixed_matrix[15][4];
  
  void build_scoremat(const std::vector<float> &base_comp, 
		      const float correction);
  //   void build_scoremat(const std::vector<float> &base_comp, 
  // 		      const size_t n_sites);
  
};

#endif
