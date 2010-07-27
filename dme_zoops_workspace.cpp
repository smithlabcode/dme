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

#include "dme_zoops_workspace.hpp"
#include "CTSet.hpp"

#include <Matrix.hpp>

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

struct dme_zoops_lextree {
  dme_zoops_lextree() : child(0) {}
  void allocate_child_ptrs();
  void insert(const string::const_iterator seq, const size_t depth,
	      const size_t motif_width);
  ~dme_zoops_lextree();
  bool remove_matching(const vector<vector<float> > &score_matrix,
		       const size_t motif_width, const size_t depth, 
		       const float current_score);
  void build(const string &sequence, const size_t motif_width);

  void get_score(const vector<vector<vector<float> > > &refined_score_matrix,
		 const vector<size_t> &path,
		 size_t depth,
		 float current_score,
		 float &max_score) const;

  void get_score(const vector<vector<float> > &score_matrix,
		 const vector<size_t> &path,
		 size_t depth,
		 float current_score,
		 float &max_score) const;
  
  dme_zoops_lextree **child;
};



/* allocate space for the array of children pointers of a dme_zoops_lextree */
void
dme_zoops_lextree::allocate_child_ptrs() {
  child = new dme_zoops_lextree *[alphabet_size];
  std::fill_n(child, alphabet_size, static_cast<dme_zoops_lextree*>(0));
}



/* insert a sequence below a subtree */
void
dme_zoops_lextree::insert(const string::const_iterator seq, const size_t depth, 
		const size_t motif_width) {
  if (depth < motif_width) {
    const size_t index = base2int(*seq);
    if (child == 0)
      allocate_child_ptrs();
    if (child[index] == 0)
      child[index] = new dme_zoops_lextree;
    child[index]->insert(seq + 1, depth + 1, motif_width);
  }
}



/* build a dme_zoops_lextree from a sets of foreground and background sequences */
void
dme_zoops_lextree::build(const string &sequence, const size_t motif_width) {
  allocate_child_ptrs();
  const size_t n_substrings = sequence.length() - motif_width + 1;
  for (size_t i = 0; i < n_substrings; ++i) {
    size_t j = 0;
    for (; j < motif_width && toupper(sequence[i + j]) != 'N'; ++j);
    if (j == motif_width)
      insert(sequence.begin() + i, 0, motif_width);
  }
}



/* recursively free space used by a lexicographic tree */
dme_zoops_lextree::~dme_zoops_lextree() {
  if (child) {
    for (size_t i = 0; i < alphabet_size; ++i) 
      if (child[i]) {
	delete child[i];
	child[i] = 0;
      }
    delete[] child;
    child = 0;
  }
}

bool
dme_zoops_lextree::remove_matching(const vector<vector<float> > &score_matrix,
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
dme_zoops_workspace::deactivate(const ScoringMatrix &sm) {
  vector<vector<float> > score_matrix(sm.get_width(),
				      vector<float>(alphabet_size));
  float deactivation_max_score = 0.0;
  for (size_t i = 0; i < sm.get_width(); ++i) {
    float max_column_score = 0.0;
    for (size_t j = 0; j < alphabet_size; ++j) {
      score_matrix[i][j] = sm[i][j];
      max_column_score = max(score_matrix[i][j], max_column_score);
    }
    for (size_t j = 0; j < alphabet_size; j++)
      score_matrix[i][j] = max_column_score - score_matrix[i][j];
    deactivation_max_score += max_column_score;
  }
  for (size_t i = 0; i < fgsize; ++i)
    fgnodes.front()[i][0]->remove_matching(score_matrix, sm.get_width(), 
					   0, deactivation_max_score);
  for (size_t i = 0; i < bgsize; ++i)
    bgtrees[i]->remove_matching(score_matrix, sm.get_width(),
				0, deactivation_max_score);
}



////////////////////////////////////////////////////////////////////////
//////////////////////////                           ///////////////////
//////////////////////////    GET BACKGROUND SCORE   ///////////////////
//////////////////////////                           ///////////////////
////////////////////////////////////////////////////////////////////////



void
dme_zoops_lextree::get_score(const vector<vector<vector<float> > > &refined_score_matrix,
		   const vector<size_t> &path,
		   size_t depth,
		   float current_score,
		   float &current_best) const {
  if (child) {
    for (size_t i = 0; i < alphabet_size; ++i)
      if (child[i] && 
	  current_score  - refined_score_matrix[depth][path[depth]][i] > 
	  current_best)
	child[i]->get_score(refined_score_matrix, path,
			    depth + 1,
			    current_score -
			    refined_score_matrix[depth][path[depth]][i],
			    current_best);
  }
  else current_best = max(current_best, current_score);
}



float
dme_zoops_workspace::get_background_refined_score() const {
  float total_score = 0;
  for (size_t i = 0; i < bgsize; ++i) {
    float current_best = 0;
    bgtrees[i]->get_score(refined_scoremat, prefix,
			  0, refined_max_score, current_best);
    if (current_best > 0)
      total_score += current_best;
  }
  return total_score;
}



void
dme_zoops_lextree::get_score(const vector<vector<float> > &score_matrix,
		   const vector<size_t> &path,
		   size_t depth,
		   float current_score,
		   float &current_best) const {
  if (child) {
    for (size_t i = 0; i < alphabet_size; ++i)
      if (child[i] && 
	  current_score  - score_matrix[path[depth]][i] > 
	  current_best)
	child[i]->get_score(score_matrix, path,
			    depth + 1,
			    current_score -
			    score_matrix[path[depth]][i],
			    current_best);
  }
  else current_best = max(current_best, current_score);
}



float
dme_zoops_workspace::get_background_score() const {
  float total_score = 0;
  for (size_t i = 0; i < bgsize; ++i) {
    float current_best = 0;
    bgtrees[i]->get_score(scoremat, prefix,
			  0, max_score, current_best);
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


void 
dme_zoops_workspace::refined_enumeration(const size_t depth,
				     const float surplus_information,
				     const size_t remaining_changes) {
  for (size_t i = 0; i == 0 ||
	 (remaining_changes > 0 &&
	  i < n_refined_types[depth - 1]); ++i)
    if (refined_coltype_bits[depth - 1][i] < surplus_information) {
      float upper_bound = 0.0;
      for (size_t j = 0; j < fgsize; ++j) {
	float current_best = 0.0;
	size_t frontier_size = 0;
	dme_zoops_lextree **lim = &fgnodes[depth - 1][j].front();
	for (size_t k = 0; lim[k] != 0; ++k) {
	  const dme_zoops_lextree *n = fgnodes[depth - 1][j][k];
	  for (size_t l = 0; l < alphabet_size; ++l)
	    if (n->child[l]) {
	      const float score = fgscore[depth - 1][j][k] - 
		refined_scoremat[depth - 1][i][l];
	      if (score > 0.0) {
		fgscore[depth][j][frontier_size] = score;
		fgnodes[depth][j][frontier_size] = n->child[l];
		current_best = max(current_best, score);
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
	  if (upper_bound - background_score*lambda > best_score) {
	    best_path = prefix;
	    best_score = upper_bound - background_score*lambda;
	  }
	}
	else
	  refined_enumeration(depth + 1, 
			      surplus_information - refined_coltype_bits[depth - 1][i],
			      remaining_changes - (i != 0));
      }
    }
}



DMEPath
dme_zoops_workspace::run_dme_zoops_local(const vector<vector<vector<float> > > &refined_types,
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
  const float surplus_information = (max_col_type_info - 
				     min_information)*motif_width;
  
  // initialize the log scoring matrix and get the max score
  refined_scoremat.resize(motif_width);
  
  // float max_score = 0.0;
  refined_max_score = 0.0;
  for (size_t i = 0; i < motif_width; ++i)
    refined_max_score += complement_scoremat(refined_types[i], 
					     refined_scoremat[i]);
  
  // initialize the max score for each fg and bg sequence 
  for (size_t i = 0; i < fgsize; ++i)
    fgscore.front()[i][0] = refined_max_score;
  
  // make sure the variables holding the best current motif and score
  // are initialized empty and zero
  best_path.clear();
  best_score = 0.0;
  
  refined_enumeration(1, surplus_information, n_changes);
  
  return DMEPath(best_path, best_score);
}



void 
dme_zoops_workspace::enumeration(const size_t depth, const float surplus_information) {
  for (size_t i = 0; i < n_types && coltype_bits[i] < surplus_information; ++i) {
    // get the foreground scores
    float upper_bound = 0.0;
    for (size_t j = 0; j < fgsize; ++j) {
      float current_best = 0.0;
      size_t frontier_size = 0;
      dme_zoops_lextree **lim = &fgnodes[depth - 1][j].front();
      for (size_t k = 0; lim[k] != 0; ++k) {
	const dme_zoops_lextree *n = lim[k];
	for (size_t l = 0; l < alphabet_size; ++l)
	  if (n->child[l]) {
	    const float score = fgscore[depth - 1][j][k] - scoremat[i][l];
	    if (score > 0.0) {
	      fgscore[depth][j][frontier_size] = score;
	      fgnodes[depth][j][frontier_size] = n->child[l];
	      current_best = max(current_best, score);
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
	if (upper_bound - background_score*lambda > best_score) {
	  best_path = prefix;
	  best_score = upper_bound - background_score*lambda;
	}
      }
      else
	enumeration(depth + 1, surplus_information - coltype_bits[i]);
    }
  }
}



DMEPath
dme_zoops_workspace::run_dme_zoops(const vector<vector<float> > &column_types,
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
  // const float max_score = complement_scoremat(column_types, scoremat);
  max_score = complement_scoremat(column_types, scoremat)*motif_width;
  
  // set the values in the first cell in the scores arrays
  for (size_t i = 0; i < fgsize; ++i)
    fgscore.front()[i][0] = max_score;
  
  best_path.clear();
  best_score = 0.0;
  
  enumeration(1, surplus_information);
  
  return DMEPath(best_path, best_score);
}



dme_zoops_workspace::dme_zoops_workspace(const vector<string> &foreground, 
					 const vector<string> &background, 
					 const size_t width,
					 const float adjustment) : 
  motif_width(width) {
  
  // initialize the motif prefix
  prefix = vector<size_t>(motif_width, 0);
  
  fgsize = foreground.size();
  bgsize = background.size();
  
  lambda = (!background.empty()) ? static_cast<float>(fgsize)/bgsize : 1;
  
  lambda *= adjustment;
  
  // allocate the tables of nodes and scores
  fgscore = vector<vector<vector<float> > >(width + 1);
  fgnodes = vector<vector<vector<dme_zoops_lextree*> > >(width + 1);
  for (size_t i = 0; i <= width; ++i) {
    fgscore[i] = vector<vector<float> >(fgsize);
    fgnodes[i] = vector<vector<dme_zoops_lextree*> >(fgsize);
    for (size_t j = 0; j < fgsize; ++j) {
      const size_t max_frontier_size =
	min(foreground[j].length(), static_cast<size_t>(pow(alphabet_size, i))) + 1;
      fgscore[i][j] = vector<float>(max_frontier_size);
      fgnodes[i][j] = vector<dme_zoops_lextree*>(max_frontier_size);
    }
  }
  for (size_t i = 0; i < fgsize; ++i) {
    fgnodes.front()[i][0] = new dme_zoops_lextree;
    fgnodes.front()[i][0]->build(foreground[i], width);
    fgnodes.front()[i][1] = 0;
  }
  // cerr << bgsize << endl;
  bgtrees = vector<dme_zoops_lextree*>(bgsize);
  for (size_t i = 0; i < bgsize; ++i) {
    bgtrees[i] = new dme_zoops_lextree;
    bgtrees[i]->build(background[i], width);
  }
}



dme_zoops_workspace::~dme_zoops_workspace() {
  for (size_t i = 0; i < fgnodes.front().size(); ++i)
    if (fgnodes.front()[i][0])
      delete fgnodes.front()[i][0];
  for (size_t i = 0; i < bgtrees.size(); ++i)
    if (bgtrees[i]) {
      delete bgtrees[i];
      bgtrees[i] = 0;
    }
}
