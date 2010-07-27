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

#include "dme_tcm_workspace.hpp"

using std::max;
using std::min;
using std::vector;
using std::string;


using std::cerr;
using std::endl;



/* convert column types to corresponding scoring matrix columns */
static float
complement_scoremat(const vector<vector<float> > coltypes, 
		    vector<vector<float> > &scoremat) {
  scoremat = coltypes;
  
  float max_score = 0.0;
  for (size_t i = 0; i < coltypes.size(); ++i)
    for (size_t j = 0; j < coltypes[i].size(); j++)
      max_score = max(max_score, scoremat[i][j]);
  
  for (size_t i = 0; i < coltypes.size(); ++i)
    for (size_t j = 0; j < coltypes[i].size(); j++)
      scoremat[i][j] = max_score - scoremat[i][j];
  
  return max_score;
}



////////////////////////////////////////////////////////////////////////
////////////////////                                ////////////////////
////////////////////    LEXICOGRAPHIC TREE STUFF    ////////////////////
////////////////////                                ////////////////////
////////////////////////////////////////////////////////////////////////

struct dme_tcm_lextree {
  dme_tcm_lextree() : child(0), max_diff(0.0) {}
  void allocate_child_ptrs();
  float set_max_diff();
  void insert(const string::const_iterator seq, const size_t depth, 
	      const size_t motif_width, const float val);
  ~dme_tcm_lextree();
  bool remove_matching(const vector<vector<float> > &score_matrix,
		       const size_t motif_width, const size_t depth, 
		       const float current_score);
  void build(const vector<string> &foreground, 
	     const vector<string> &background,
	     const float lambda, const size_t motif_width);
  
  dme_tcm_lextree **child;
  float max_diff;
};



/* allocate space for the array of children pointers of a dme_tcm_lextree */
void
dme_tcm_lextree::allocate_child_ptrs() {
  child = new dme_tcm_lextree *[alphabet_size];
  std::fill_n(child, alphabet_size, static_cast<dme_tcm_lextree*>(0));
}



/* insert a sequence below a subtree */
void
dme_tcm_lextree::insert(const string::const_iterator seq, const size_t depth, 
		const size_t motif_width, const float val) {
  if (depth < motif_width) {
    const size_t index = base2int(*seq);
    if (child == 0)
      allocate_child_ptrs();
    if (child[index] == 0)
      child[index] = new dme_tcm_lextree;
    child[index]->insert(seq + 1, depth + 1, motif_width, val);
    child[index]->max_diff += val;
  }
}



float
dme_tcm_lextree::set_max_diff() {
  if (child) {
    max_diff = static_cast<float>(0);
    for (size_t i = 0; i < alphabet_size; ++i)
      if (child[i]) {
	child[i]->set_max_diff();
	max_diff += max(child[i]->max_diff, static_cast<float>(0));
      }
  }
  return max_diff;
}



/* build a dme_tcm_lextree from a sets of foreground and background sequences */
void
dme_tcm_lextree::build(const vector<string> &foreground, const vector<string> &background,
	       const float lambda, const size_t motif_width) {
  allocate_child_ptrs();
  for (size_t i = 0; i < foreground.size(); i++) {
    const size_t n_substrings = foreground[i].length() - motif_width + 1;
    for (size_t j = 0; j < n_substrings; ++j) {
      size_t k = 0;
      for (; k < motif_width && toupper(foreground[i][j+k]) != 'N'; k++);
      if (k == motif_width)
	insert(foreground[i].begin() + j, 0, motif_width, 1.0);
    }
  }
  for (size_t i = 0; i < background.size(); i++) {
    const size_t n_substrings = background[i].length() - motif_width + 1;
    for (size_t j = 0; j < n_substrings; ++j) {
      size_t k = 0;
      for (; k < motif_width && toupper(background[i][j + k]) != 'N'; k++);
      if (k == motif_width)
	insert(background[i].begin() + j, 0, motif_width, -lambda);
    }
  }
  set_max_diff();
}



/* recursively free space used by a lexicographic tree */
dme_tcm_lextree::~dme_tcm_lextree() {
  if (child) {
    for (size_t i = 0; i < alphabet_size; ++i) 
      if (child[i]) 
	delete child[i];
    delete[] child;
  }
}

bool
dme_tcm_lextree::remove_matching(const vector<vector<float> > &score_matrix,
			 const size_t motif_width, const size_t depth, 
			 const float current_score) {
  size_t n_inactivated = 0;
  for (size_t i = 0; i < alphabet_size; ++i) {
    if (child[i]) {
      if (current_score  - score_matrix[depth][i] > 0) {
	if (depth == motif_width - 1 ||
	    child[i]->remove_matching(score_matrix, motif_width,
				      depth + 1, current_score -
				      score_matrix[depth][i])) {
	  n_inactivated++;
	  delete child[i];
	  child[i] = 0;
	}
      }
    }
    else n_inactivated++;
  }
  // return true if no more children exist for the current node
  return (n_inactivated == alphabet_size);
}

////////////////////////////////////////////////////////////////////////
/////////////////////                              /////////////////////
/////////////////////   ERASING PREVIOUS MOTIFS    /////////////////////
/////////////////////                              /////////////////////
////////////////////////////////////////////////////////////////////////

/**
 * This function removes all nodes in the tree for which all leaves
 * below them correspond to strings that match the given matrix.
 */
void
dme_tcm_workspace::deactivate(const ScoringMatrix &sm) {
  vector<vector<float> > score_matrix(sm.get_width(),
				      vector<float>(alphabet_size));
  float max_score = 0.0;
  for (size_t i = 0; i < sm.get_width(); ++i) {
    float max_column_score = 0.0;
    for (size_t j = 0; j < alphabet_size; ++j) {
      score_matrix[i][j] = sm[i][j];
      max_column_score = max(score_matrix[i][j], max_column_score);
    }
    for (size_t j = 0; j < alphabet_size; j++)
      score_matrix[i][j] = max_column_score - score_matrix[i][j];
    max_score += max_column_score;
  }
  nodes[0][0]->remove_matching(score_matrix, sm.get_width(), 0, max_score);
  nodes[0][0]->set_max_diff();
}


////////////////////////////////////////////////////////////////////////
//////////////////////////                  ////////////////////////////
//////////////////////////    SEARCH STUFF  ////////////////////////////
//////////////////////////                  ////////////////////////////
////////////////////////////////////////////////////////////////////////


void
dme_tcm_workspace::refined_enumeration(const size_t depth,
				   const size_t prev_frontier,
				   const float surplus_information,
				   const size_t remaining_changes) {
  for (size_t i = 0; i == 0 ||
	 (remaining_changes > 0 &&
	  i < n_refined_types[depth - 1]); ++i)
    if (refined_coltype_bits[depth - 1][i] < surplus_information) {
      float upper_bound = 0.0;
      size_t frontier_size = 0;
      for (size_t j = 0; j < prev_frontier; ++j) {
	const dme_tcm_lextree *n = nodes[depth - 1][j];
	for (size_t k = 0; k < alphabet_size; ++k) {
	  if (n->child[k]) {
	    const float tmp_score = score[depth - 1][j] - 
	      refined_scoremat[depth - 1][i][k];
	    if (tmp_score > 0.0) {
	      score[depth][frontier_size] = tmp_score;
	      nodes[depth][frontier_size] = n->child[k];
	      upper_bound += tmp_score*n->child[k]->max_diff;
	      frontier_size++;
	    }
	  }
	}
      }
      if (upper_bound > best_score) {
	prefix[depth - 1] = i;
	if (depth == motif_width) {
	  best_path = prefix;
	  best_score = upper_bound;
	}
	else
	  refined_enumeration(depth + 1, frontier_size,
			      surplus_information - refined_coltype_bits[depth - 1][i],
			      remaining_changes - (i != 0));
      }
    }
}



DMEPath
dme_tcm_workspace::run_dme_tcm_local(const vector<vector<vector<float> > > &refined_types,
			     const vector<vector<float> > &refined_bits,
			     const float min_information,
			     const size_t n_changes) {
  n_refined_types.resize(refined_types.size());
  for (size_t i = 0; i < refined_types.size(); ++i)
    n_refined_types[i] = refined_types[i].size();
  
  // get the surplus_information content matrix
  float max_col_type_info = 0.0;
  for (size_t i = 0; i < motif_width; ++i)
    for (size_t j = 0; j < n_refined_types[i]; j++)
      max_col_type_info = max(max_col_type_info, refined_bits[i][j]);

  refined_coltype_bits = refined_bits;
  for (size_t i = 0; i < motif_width; ++i)
    for (size_t j = 0; j < n_refined_types[i]; ++j)
      refined_coltype_bits[i][j] = max_col_type_info - 
	refined_coltype_bits[i][j];
  
  // Set the minimum bits/column
  const float surplus_information = (max_col_type_info - min_information)*motif_width;
  
  // initialize the log scoring matrix and get the max score
  refined_scoremat.resize(motif_width);
  score[0][0] = 0.0;
  for (size_t i = 0; i < motif_width; ++i)
    score[0][0] += complement_scoremat(refined_types[i], 
				       refined_scoremat[i]);
  
  // make sure the variables holding the best current motif and score
  // are initialized empty and zero
  best_path.clear();
  best_score = 0.0;
  
  refined_enumeration(1, 1, surplus_information, n_changes);
  
  return DMEPath(best_path, best_score);
}




void 
dme_tcm_workspace::enumeration(const size_t depth, const size_t prev_frontier,
			   const float surplus_information) {
  for (size_t i = 0; i < n_types && coltype_bits[i] <= surplus_information; ++i) {
    float upper_bound = 0.0;
    size_t frontier_size = 0;
    for (size_t j = 0; j < prev_frontier; j++) {
      const dme_tcm_lextree *n = nodes[depth - 1][j];
      for (size_t k = 0; k < alphabet_size; ++k)
	if (n->child[k]) {
	  const float tmp_score = score[depth - 1][j] - scoremat[i][k];
	  if (tmp_score > 0) {
	    score[depth][frontier_size] = tmp_score;
	    nodes[depth][frontier_size] = n->child[k];
	    upper_bound += tmp_score*n->child[k]->max_diff;
	    frontier_size++;
	  }
	}
    }
    if (upper_bound > best_score) {
      prefix[depth - 1] = i;
      if (depth == motif_width) {
	best_path = prefix;
	best_score = upper_bound;
      }
      else
	enumeration(depth + 1, frontier_size, surplus_information - coltype_bits[i]);
    }
  }
}



DMEPath
dme_tcm_workspace::run_dme_tcm(const vector<vector<float> > &column_types,
			       const vector<float > &bits,
			       const float min_information) {
  // set the global alphabet and alphabet size variables
  n_types = column_types.size();
  
  // get the info content matrix
  coltype_bits = bits;
  
  // get the surplus_information content of each column type and the maximum
  // possible surplus_information in a column
  float max_col_type_info = 0.0;
  for (size_t i = 0; i < n_types; ++i)
    max_col_type_info = max(max_col_type_info, coltype_bits[i]);
  for (size_t i = 0; i < n_types; ++i)
    coltype_bits[i] = max_col_type_info - coltype_bits[i];
  
  // set a bound on how much the surplus_information can go under the maximum
  const float surplus_information = (max_col_type_info -
				     min_information)*motif_width;
  
  // get scoring matrix
  const float max_score = complement_scoremat(column_types, scoremat);
  
  // set the value in the first cell in the scores array
  score[0][0] = max_score*motif_width;
  
  best_path.clear();
  best_score = 0.0;
  
  enumeration(1, 1, surplus_information);
  
  return DMEPath(best_path, best_score);
}



dme_tcm_workspace::dme_tcm_workspace(const vector<string> &foreground, 
				     const vector<string> &background, 
				     const int width,
				     const float adjustment) : 
  motif_width(width) {
  
  prefix = vector<size_t>(motif_width, 0);
  
  // get the foreground and background sequence lengths
  size_t fg_seqlen = 0;
  for (size_t i = 0; i < foreground.size(); ++i)
    fg_seqlen += (foreground[i].length() - 
		  std::count(foreground[i].begin(), foreground[i].end(), 'N'));
  
  size_t bg_seqlen = 0;
  for (size_t i = 0; i < background.size(); ++i)
    bg_seqlen += (background[i].length() - 
		  std::count(background[i].begin(), background[i].end(), 'N'));
  
  const float lambda = ((bg_seqlen > 0) ? 
			static_cast<float>(fg_seqlen)/bg_seqlen : 1)*adjustment;

  const size_t total_sites = fg_seqlen + bg_seqlen;

  // allocate the tables of nodes and scores
  score = vector<vector<float> >(motif_width + 1);
  nodes = vector<vector<dme_tcm_lextree*> >(motif_width + 1);
  for (size_t i = 0; i <= motif_width; ++i) {
    const size_t max_frontier_size = 
      min(total_sites, static_cast<size_t>(pow(alphabet_size, i)));
    score[i] = vector<float>(max_frontier_size);
    nodes[i] = vector<dme_tcm_lextree*>(max_frontier_size);
  }
  nodes[0][0] = new dme_tcm_lextree;
  nodes[0][0]->build(foreground, background, lambda, motif_width);
}

dme_tcm_workspace::~dme_tcm_workspace() {
  delete nodes[0][0];
}
