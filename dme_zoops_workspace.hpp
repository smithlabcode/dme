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

/**************************************
 * HEADER FOR THE DME ZOOPS WORKSPACE *
 **************************************/

#ifndef DME_ZOOPS_WORKSPACE_HPP
#define DME_ZOOPS_WORKSPACE_HPP

#include <cread.hpp>
#include <Alphabet.hpp>
#include <ScoringMatrix.hpp>

#include "CTSet.hpp"

struct dme_zoops_lextree;

class dme_zoops_workspace {

public:
  dme_zoops_workspace(const std::vector<std::string> &foreground, 
		      const std::vector<std::string> &background, 
		      const size_t motif_width_param,
		      const float adjustment);
  ~dme_zoops_workspace();
  
  DMEPath
  run_dme_zoops(const std::vector<std::vector<float> > &column_types, 
	    const std::vector<float> &bits,
	    const float min_information);
  
  DMEPath 
  run_dme_zoops_local(const std::vector<std::vector<std::vector<float> > > &refined_column_types,
		  const std::vector<std::vector<float> > &refined_bits,
		  const float min_information,
		  const size_t n_changes);
  
  void deactivate(const ScoringMatrix &sm);

private:
  void refined_enumeration(const size_t depth, 
			   const float surplus_information,
			   const size_t remaining_changes);
  void enumeration(const size_t depth, const float surplus_information);

  float get_background_score() const;
  float get_background_refined_score() const;
  
  // These are set initially and not changed
  std::vector<std::vector<std::vector<float> > > fgscore;
  std::vector<std::vector<std::vector<dme_zoops_lextree*> > > fgnodes;
  std::vector<dme_zoops_lextree*> bgtrees;
  
  size_t fgsize;
  size_t bgsize;
  
  size_t motif_width;          // width of motif
  std::vector<size_t> prefix;
  
  // standard search variables
  std::vector<std::vector<float> > scoremat;
  std::vector<float> coltype_bits;
  size_t n_types;
  float max_score;
  
  // local search variables
  std::vector<std::vector<std::vector<float> > > refined_scoremat;
  std::vector<size_t> n_refined_types;
  std::vector<std::vector<float> > refined_coltype_bits;
  float refined_max_score;

  // for both standard and local search (reset right away at each
  // public function call)
  std::vector<size_t> best_path;
  float best_score;
  float lambda;
};

#endif
