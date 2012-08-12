/*
 * Copyright (C) 2012 University of Southern California and 
 *                    Andrew D Smith
 *
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

#include <iomanip>
#include <numeric>

#include <OptionParser.hpp>

#include <smithlab_os.hpp>
#include <smithlab_utils.hpp>

#include "dme_tcm_workspace.hpp"
#include "dme_zoops_workspace.hpp"
#include "CTSet.hpp"

#include "Matrix.hpp"
#include "ScoringMatrix.hpp"
#include "Motif.hpp"

using std::string;
using std::vector;
using std::pair;
using std::ostream;
using std::ofstream;
using std::ostringstream;
using std::ostream_iterator;
using std::cout;
using std::cerr;
using std::endl;
using std::numeric_limits;
using std::transform;
using std::accumulate;
using std::max;
using std::copy;
using std::ios_base;

using smithlab::alphabet_size;

bool VERBOSE = false;

void
get_seeds_zoops(const bool single_strand,
		const vector<string> &foreground,
		const vector<string> &background,
		const vector<float> &base_comp,
		const size_t motif_width,
		const size_t outputs,
		const float granularity,
		const float bits,
		const float correction,
		const float adjustment,
		vector<Matrix> &seeds) {
  
  static const char *seeds_progress_prefix = "obtaining seeds  ";
  
  dme_zoops_workspace ws(foreground, background, motif_width,
			 adjustment);
  
  CTSet cts(granularity, base_comp, correction);
  
  for (size_t i = 0; i < outputs; ++i) {
    if (VERBOSE)
      cerr << "\r" << seeds_progress_prefix 
	   << i << "/" << outputs;
    const DMEPath path(ws.run_dme_zoops(cts.get_matrix(), cts.get_bits(), bits));
    if (path.score == 0.0)
      break;
    // "erase" the discovered motif from the sequence sets
    const Matrix matrix(cts.path_to_matrix(path));
    const ScoringMatrix sm(matrix, base_comp, correction);
      
    ws.deactivate(sm);
    if (!single_strand)
      ws.deactivate(sm.revcomp());
    seeds.push_back(cts.path_to_matrix(path));
  }
  if (VERBOSE)
    cerr << "\r" << seeds_progress_prefix
	 << outputs << "/" << outputs << endl;
}



void
refine_matrix_zoops(dme_zoops_workspace &ws, const size_t motif_width,
		    const float granularity, const float information,
		    const vector<float> &base_comp, 
		    const size_t n_changes,
		    const size_t n_iterations,
		    const float required_improvement,
		    const string &progress_prefix,
		    Matrix &matrix) {
  
  float score = 0, prev_score = 0, improvement = 1;
  
  for (size_t i = 0; i < n_iterations && 
	 improvement > required_improvement; ++i) {
    
    prev_score = score;
    
    vector<CTSet> refined_cts;
    for (size_t j = 0; j < matrix.get_width(); ++j)
      refined_cts.push_back(CTSet(matrix[j], granularity, base_comp));
    
    vector<size_t> c_count;
    transform(refined_cts.begin(), refined_cts.end(), back_inserter(c_count),
	      std::mem_fun_ref(&CTSet::size));
    
    vector<vector<vector<float> > > matrix_array(refined_cts.size());
    vector<vector<float> > refined_bits(refined_cts.size());
    for (size_t j = 0; j < refined_cts.size(); ++j) {
      matrix_array[j] = refined_cts[j].get_matrix();
      refined_bits[j] = refined_cts[j].get_bits();
    }
    
    const DMEPath refined_path(ws.run_dme_zoops_local(matrix_array, refined_bits,
						      information, n_changes));
    score = refined_path.score;
    improvement = max(static_cast<float>(0), (score - prev_score)/score);
    
    if (VERBOSE) {
      cerr.setf(ios_base::left, ios_base::adjustfield);
      cerr << "\r" << progress_prefix << "    "
	   << "score=" << std::fixed << std::setprecision(2) 
	   << std::setw(10) << refined_path.score << " "
	   << "delta=" 
	   << std::scientific
	   << std::setprecision(1) 
	   << std::setw(10) << improvement;
    }
    
    const Matrix refined(CTSet::path_to_matrix(refined_path, refined_cts));

    if (refined.get_width() > 0) 
      matrix = refined;
  }
}



void
refine_matrices_zoops(const bool single_strand,
		      const vector<string> &foreground,
		      const vector<string> &background,
		      const size_t motif_width,
		      const float refine_granularity,
		      const float bits,
		      const vector<float> &base_comp, 
		      const size_t n_changes,
		      const size_t n_iterations,
		      const float required_improvement,
		      const float correction,
		      const float adjustment,
		      vector<Matrix> &seeds) {
  static const char *refining_progress_prefix = "refining motifs  ";
  
  dme_zoops_workspace ws(foreground, background, motif_width, adjustment);
  for (size_t i = 0; i < seeds.size(); ++i) {
    const string progress_prefix(refining_progress_prefix + 
				 toa(i + 1) + "/" +
				 toa(seeds.size()));
    refine_matrix_zoops(ws, motif_width, refine_granularity,
			bits, base_comp, n_changes,
			n_iterations, required_improvement,
			progress_prefix, seeds[i]);
    
    const ScoringMatrix sm(seeds[i], base_comp, correction);
    ws.deactivate(sm);
    if (!single_strand)
      ws.deactivate(sm.revcomp());
  }
  if (VERBOSE) {
    const string message(refining_progress_prefix + 
			 toa(seeds.size()) + "/" + 
			 toa(seeds.size()));
    cerr << '\r' << std::left << std::setfill(' ') << std::setw(72)
	 << message << endl;
  }
}



void
get_seeds_tcm(const bool single_strand,
	      const vector<string> &foreground,
	      const vector<string> &background,
	      const vector<float> &base_comp,
	      const size_t motif_width,
	      const size_t outputs,
	      const float granularity,
	      const float bits,
	      const float correction,
	      const float adjustment,
	      vector<Matrix> &seeds) {
  
  static const char *seeds_progress_prefix = "obtaining seeds  ";
  
  dme_tcm_workspace ws(foreground, background, motif_width, adjustment);
  // construct the column-type set
  CTSet cts(granularity, base_comp, correction);
  
  for (size_t i = 0; i < outputs; ++i) {
    if (VERBOSE)
      cerr << "\r" << seeds_progress_prefix 
	   << i << "/" << outputs;
    const DMEPath path(ws.run_dme_tcm(cts.get_matrix(), cts.get_bits(), bits));
    if (path.score == 0.0)
      break;
    // "erase" the discovered motif from the sequence sets
    const Matrix matrix(cts.path_to_matrix(path));
    const ScoringMatrix sm(matrix, base_comp, correction);
      
    ws.deactivate(sm);
    if (!single_strand)
      ws.deactivate(sm.revcomp());
    seeds.push_back(cts.path_to_matrix(path));
  }
  if (VERBOSE)
    cerr << "\r" << seeds_progress_prefix
	 << outputs << "/" << outputs << endl;
}



void
refine_matrix_tcm(dme_tcm_workspace &ws, const size_t motif_width,
		  const float granularity, const float information,
		  const vector<float> &base_comp, 
		  const size_t n_changes,
		  const size_t n_iterations,
		  const float required_improvement,
		  const string &progress_prefix,
		  Matrix &matrix) {
  
  float score = 0, prev_score = 0, improvement = 1;
  
  for (size_t i = 0; i < n_iterations && 
	 improvement > required_improvement; ++i) {
    
    prev_score = score;
    
    vector<CTSet> refined_cts;
    for (size_t j = 0; j < matrix.get_width(); ++j)
      refined_cts.push_back(CTSet(matrix[j], granularity, base_comp));
    
    vector<size_t> c_count;
    transform(refined_cts.begin(), refined_cts.end(), back_inserter(c_count),
	      std::mem_fun_ref(&CTSet::size));
    
    vector<vector<vector<float> > > matrix_array(refined_cts.size());
    vector<vector<float> > refined_bits(refined_cts.size());
    for (size_t j = 0; j < refined_cts.size(); ++j) {
      matrix_array[j] = refined_cts[j].get_matrix();
      refined_bits[j] = refined_cts[j].get_bits();
    }
    
    const DMEPath refined_path(ws.run_dme_tcm_local(matrix_array, refined_bits,
						    information, n_changes));
    score = refined_path.score;
    improvement = max(static_cast<float>(0), (score - prev_score)/score);
    
    if (VERBOSE) {
      cerr.setf(ios_base::left, ios_base::adjustfield);
      cerr << "\r" << progress_prefix << "    "
	   << "score=" << std::fixed << std::setprecision(2) 
	   << std::setw(10) << refined_path.score << " "
	   << "delta=" 
	   << std::scientific 
	   << std::setprecision(1) 
	   << std::setw(10) << improvement;
    }
    
    const Matrix refined(CTSet::path_to_matrix(refined_path, refined_cts));
    if (refined.get_width() > 0) 
      matrix = refined;
  }
}



void
refine_matrices_tcm(const bool single_strand,
		    const vector<string> &foreground,
		    const vector<string> &background,
		    const size_t motif_width,
		    const float refine_granularity,
		    const float bits,
		    const vector<float> &base_comp, 
		    const size_t n_changes,
		    const size_t n_iterations,
		    const float required_improvement,
		    const float correction,
		    const float adjustment,
		    vector<Matrix> &seeds) {
  static const char *refining_progress_prefix = "refining motifs  ";
  
  dme_tcm_workspace ws(foreground, background, motif_width,
		       adjustment);
  for (size_t i = 0; i < seeds.size(); ++i) {
    const string progress_prefix(refining_progress_prefix + 
				 toa(i + 1) + "/" +
				 toa(seeds.size()));
    refine_matrix_tcm(ws, motif_width, refine_granularity, 
		      bits, base_comp, n_changes,
		      n_iterations, required_improvement,
		      progress_prefix, seeds[i]);
    
    const ScoringMatrix sm(seeds[i], base_comp, correction);
    ws.deactivate(sm);
    if (!single_strand)
      ws.deactivate(sm.revcomp());
  }
  if (VERBOSE) {
    const string message(refining_progress_prefix + 
			 toa(seeds.size()) + "/" + toa(seeds.size()));
    cerr << '\r' << std::left << std::setfill(' ') << std::setw(72)
	 << message << endl;
  }
}



size_t 
effective_sequence_length(const vector<string>& s) {
  size_t length = 0;
  for (vector<string>::const_iterator i = s.begin(); i != s.end(); ++i)
    length += count_if(i->begin(), i->end(), &isvalid);
  return length;
}



void
mask(size_t width, size_t max_order, vector<string> &seqs) {
  const char masked_base = 'N';
  for (size_t i = 1; i <= max_order; ++i)
    for (size_t j = 0; j < seqs.size(); ++j) {
      size_t rep = i;
      for (size_t k = i; k < seqs[j].length(); ++k)
	if (seqs[j][k] == seqs[j][k - i] && seqs[j][k] != masked_base &&
	    k < seqs[j].length() - 1)
	  rep++;
	else {
	  if (rep > width)
	    std::fill_n(seqs[j].begin() + k - rep, rep, masked_base);
	  rep = i;
	}
    }
}



string
degenerate_consensus(const Matrix &matrix) {
  static const float correction = 0.0000000001;
  static const char degenerate_bases[] = {'A','C','G','T', 'M', 'R', 'W', 'S', 
					  'Y','K','V','H', 'D', 'B', 'N'};

  static const size_t n_degen_nucs = 15;
  static float fixed_matrix[15][4] = {
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
  
  string consensus;
  for (size_t i = 0; i < matrix.get_width(); ++i) {
    float score = numeric_limits<float>::max();
    size_t best = 0;
    for (size_t j = 0; j < n_degen_nucs; ++j) {
      float temp_score = 0.0;
      for (size_t k = 0; k < alphabet_size; ++k) {
	const float matval = (fixed_matrix[j][k] > 0) ? 
	  fixed_matrix[j][k] : correction;
	const float freq = (matrix[i][k] > 0) ? matrix[i][k] : correction;
	temp_score += (matval - freq) * (log2(matval) - log2(freq));
      }
      if (temp_score < score) {
	best = j;
	score = temp_score;
      }
    }
    consensus += degenerate_bases[best];
  }
  return consensus;
}



struct Site {
  size_t seq;
  size_t pos;
  float score;
  bool strand;
  Site(size_t se, size_t p, float sc, bool st) :
    seq(se), pos(p), score(sc), strand(st) {}
  char strand_char() {return (strand) ? 'p' : 'n';}
  static float accumulate_score(float x, const Site& h) {return x + h.score;}
  bool operator>(const Site &s) const {return score > s.score;}
};



Matrix
get_counts_matrix(const vector<Site>& sites, 
		  const vector<string>& seqs,
		  const size_t motif_width) {
  float *counts[motif_width];
  for (size_t i = 0; i < motif_width; ++i) {
    counts[i] = new float[alphabet_size];
    std::fill(counts[i], counts[i] + alphabet_size, 0);
  }
  for (vector<Site>::const_iterator i = sites.begin(); i != sites.end(); ++i) {
    string site = seqs[i->seq].substr(i->pos, motif_width);
    if (!i->strand)
      site = revcomp(site);
    for (size_t j = 0; j < motif_width; ++j)
      counts[j][base2int(site[j])]++;
  }
  const Matrix counts_matrix(counts, motif_width);
  for (size_t i = 0; i < motif_width; ++i)
    delete[] counts[i];
  return counts_matrix;
}



bool
valid_base_id(int c) {
  return (c < static_cast<int>(alphabet_size) && c >= 0);
}

bool 
valid_subsequence(const int* offset, const size_t matwidth) {
  return count_if(offset, offset + matwidth, 
		  std::ptr_fun(&valid_base_id)) == static_cast<int>(matwidth);
}



float 
match_matrix(const ScoringMatrix& sm, const int* offset) {
  const size_t width = sm.get_width();
  float score = 0;
  for (size_t i = 0; i < width; ++i)
    score += sm[i][offset[i]];
  return score;
}



void 
get_sites_zoops(const vector<string>& sequences,
		const ScoringMatrix& sm, 
		const ScoringMatrix& smrc,
		const bool singlestrand, 
		vector<Site>& sites) {
  const size_t matwidth = sm.size();
  for (size_t i = 0; i < sequences.size(); ++i) {
    vector<int> helper(sequences[i].length());
    transform(sequences[i].begin(), sequences[i].end(), 
	      helper.begin(), &base2int);
    const size_t lim = std::max(static_cast<int>(sequences[i].length()) -
				static_cast<int>(matwidth) + 1, 0);
    vector<Site> tied;
    float best_score = std::numeric_limits<float>::min();
    for (size_t j = 0; j < lim; ++j) {
      const int* offset = &helper[j];
      if (valid_subsequence(offset, matwidth)) {
	float score = match_matrix(sm, offset);
	if (score == best_score)
	  tied.push_back(Site(i, j, score, true));
	else if (score > best_score) {
	  tied.clear();
	  tied.push_back(Site(i, j, score, true));
	  best_score = score;
	}
	float score_rc = match_matrix(smrc, offset);
	if (score_rc == best_score)
	  tied.push_back(Site(i, j, score_rc, false));
	else if (score_rc > best_score) {
	  tied.clear();
	  tied.push_back(Site(i, j, score_rc, false));
	  best_score = score_rc;
	}
      }
    }
    copy(tied.begin(), tied.end(), back_inserter(sites));
  }
}



void 
get_sites_tcm(const vector<string>& sequences,
	      const ScoringMatrix& sm, const ScoringMatrix& smrc,
	      const bool singlestrand, vector<Site>& sites) {
  const size_t matwidth = sm.size();
  for (size_t i = 0; i < sequences.size(); ++i) {

    vector<int> helper(sequences[i].length());
    transform(sequences[i].begin(), sequences[i].end(),
	      helper.begin(), &base2int);
    const size_t lim = std::max(static_cast<int>(sequences[i].length()) - 
				static_cast<int>(matwidth) + 1, 0);
    for (size_t j = 0; j < lim; ++j) {
      const int* offset = &helper[j];
      if (valid_subsequence(offset, matwidth)) {
	float score = match_matrix(sm, offset);
	if (score >= 0)
	  sites.push_back(Site(i, j, score, true));
	if (!singlestrand) {
	  score = match_matrix(smrc, offset);
	  if (score >= 0.0)
	    sites.push_back(Site(i, j, score, false));
	}
      }
    }
  }
}



float
get_sites_zoops_score(const vector<Site>& sites) {
  float score = 0;
  float best_in_seq = -numeric_limits<float>::max();
  for (size_t i = 0; i < sites.size(); ++i) {
    best_in_seq = max(sites[i].score, best_in_seq);
    if (i == sites.size() - 1 || sites[i].seq != sites[i + 1].seq) {
      score += best_in_seq;
      best_in_seq = -numeric_limits<float>::max();
    }
  }
  return score;
}



Motif 
prepare_motif_zoops(const string &name, 
		    const Matrix matrix, 
		    const vector<float> base_comp,
		    const float correction,
		    const vector<string> &fgseqs, 
		    const vector<string> &fgnames, 
		    const vector<string> &bgseqs, 
		    const vector<string> &bgnames,
		    const bool singlestrand,
		    const float fg_bg_ratio) {
  
  // convert the pwm to a scoring matrix
  const ScoringMatrix sm(matrix, base_comp, correction);
  const ScoringMatrix smrc(sm.revcomp());
  
  // get sites from the foreground
  vector<Site> fg_sites;
  get_sites_zoops(fgseqs, sm, smrc, singlestrand, fg_sites);
  
  float score = get_sites_zoops_score(fg_sites);
  
  // deal with the background (if it exists)
  vector<Site> bg_sites;
  if (!bgseqs.empty())
    get_sites_zoops(bgseqs, sm, smrc, singlestrand, bg_sites);
  
  score -= get_sites_zoops_score(bg_sites)*fg_bg_ratio;
  
  // Build the motif
  Motif motif(get_counts_matrix(fg_sites, fgseqs, matrix.get_width()));
  motif.set_accession(name);
  motif.set_identifier(degenerate_consensus(matrix)); 
  motif.set_attribute("SCORE", score);
  motif.set_attribute("FGCOUNT", fg_sites.size());
  motif.set_attribute("BGCOUNT", bg_sites.size());
  motif.set_attribute("CORRECTEDBGCOUNT", 
		      static_cast<size_t>(bg_sites.size()*fg_bg_ratio));
  motif.set_attribute("INFO", matrix.info(base_comp)/matrix.get_width());
  
  // Set the sites
  for (vector<Site>::iterator i = fg_sites.begin(); i != fg_sites.end(); ++i) {
    string temp(fgseqs[i->seq].substr(i->pos, matrix.get_width()));
    string site((i->strand) ? temp : revcomp(temp));
    motif.add_site(MotifSite(site, fgnames[i->seq], i->pos, matrix.get_width(),
			     " ", i->strand_char(), i->score));
  }
  return motif;
} // END prepare_motif_zoops()



void
get_base_comp(const vector<string>& sequences, vector<float>& base_comp) {
  vector<size_t> count(alphabet_size, 0);
  for (vector<string>::const_iterator i = sequences.begin();
       i != sequences.end(); ++i)
    for (string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (isvalid(*j)) {
	count[base2int(*j)]++;
      }
  const float total = std::accumulate(count.begin(), count.end(), 0.0);
  base_comp.clear();
  transform(count.begin(), count.end(), back_inserter(base_comp), 
		 std::bind2nd(std::divides<float>(), total));
}



void
preprocess_sequences_zoops(const bool single_strand,
			   const string &fgfilename, 
			   const string &bgfilename, 
			   vector<string> &fgnames, 
			   vector<string> &original_foreground, 
			   vector<string> &foreground,
			   vector<string> &bgnames, 
			   vector<string> &original_background, 
			   vector<string> &background,
			   vector<float> &base_comp, 
			   float &fg_bg_ratio) {
  
  // read fg sequences
  read_fasta_file(fgfilename, fgnames, original_foreground);
  
  static const size_t max_mask_order = 2;
  static const size_t decoy_width = 8;
  mask(decoy_width, max_mask_order, original_foreground);
  
  vector<float> fg_base_comp;
  get_base_comp(original_foreground, fg_base_comp);
  const size_t fg_length = effective_sequence_length(original_foreground);
  
  if (!bgfilename.empty()) {
    // read bg sequences
    read_fasta_file(bgfilename, bgnames, original_background);
    
    mask(decoy_width, max_mask_order, original_background);
    
    vector<float> bg_base_comp;
    get_base_comp(original_background, bg_base_comp);
    const size_t bg_length = effective_sequence_length(original_background);
    for (size_t i = 0; i < alphabet_size; ++i)
      base_comp.push_back((fg_length*fg_base_comp[i] +
			   bg_length*bg_base_comp[i])/(fg_length + bg_length));
    
    foreground = original_foreground;
    if (!single_strand)
      for (size_t i = 0; i < original_foreground.size(); i++)
	foreground[i] += (string("N") + revcomp(foreground[i]));
    
    background = original_background;
    if (!single_strand)
      for (size_t i = 0; i < background.size(); i++)
	background[i] += (string("N") + revcomp(background[i]));
    
    fg_bg_ratio = static_cast<float>(foreground.size())/background.size();
  }
  else {
    base_comp = fg_base_comp;
    foreground = original_foreground;
    if (!single_strand)
      for (size_t i = 0; i < foreground.size(); ++i)
	foreground[i] += "N" + revcomp(foreground[i]);
  }
  
} // END preprocess_sequences_zoops()



Motif 
prepare_motif_tcm(const string &name, 
		  const Matrix matrix, const vector<float> base_comp,
		  const float correction,
		  const vector<string> &fgseqs, 
		  const vector<string> &fgnames, 
		  const vector<string> &bgseqs, 
		  const vector<string> &bgnames,
		  const bool singlestrand,
		  const float length_ratio) {
  
  // convert the pwm to a scoring matrix
  const ScoringMatrix sm(matrix, base_comp, correction);
  const ScoringMatrix smrc(sm.revcomp());
  
  // get sites from the foreground
  vector<Site> fg_sites;
  get_sites_tcm(fgseqs, sm, smrc, singlestrand, fg_sites);
  float score(accumulate(fg_sites.begin(), fg_sites.end(), 0.0,
			 std::ptr_fun(Site::accumulate_score)));
  
  // deal with the background (if it exists)
  vector<Site> bg_sites;
  if (!bgseqs.empty())
    get_sites_tcm(bgseqs, sm, smrc, singlestrand, bg_sites);
  score -= accumulate(bg_sites.begin(), bg_sites.end(), 0.0,
		      std::ptr_fun(Site::accumulate_score))*length_ratio;
  
  // Build the motif
  Motif motif(get_counts_matrix(fg_sites, fgseqs, matrix.get_width()));
  motif.set_accession(name);
  motif.set_identifier(degenerate_consensus(matrix)); 
  motif.set_attribute("SCORE", score);
  motif.set_attribute("FGCOUNT", fg_sites.size());
  motif.set_attribute("BGCOUNT", bg_sites.size());
  motif.set_attribute("CORRECTEDBGCOUNT", 
		      static_cast<size_t>(bg_sites.size()*length_ratio));
  motif.set_attribute("INFO", matrix.info(base_comp)/matrix.get_width());
  
  // Set the sites
  for (vector<Site>::iterator i = fg_sites.begin(); i != fg_sites.end(); ++i) {
    string temp(fgseqs[i->seq].substr(i->pos, matrix.get_width()));
    string site((i->strand) ? temp : revcomp(temp));
    motif.add_site(MotifSite(site, fgnames[i->seq], i->pos, matrix.get_width(),
			     " ", i->strand_char(), i->score));
  }
  return motif;
} // END prepare_motif()



void
preprocess_sequences_tcm(const bool single_strand,
			 const string &fgfilename, 
			 const string &bgfilename, 
			 vector<string> &fgnames, 
			 vector<string> &original_foreground, 
			 vector<string> &foreground,
			 vector<string> &bgnames, 
			 vector<string> &original_background, 
			 vector<string> &background,
			 vector<float> &base_comp, float &length_ratio) {
  
  // read fg sequences
  read_fasta_file(fgfilename, fgnames, original_foreground);
  
  static const size_t max_mask_order = 2;
  static const size_t decoy_width = 8;
  mask(decoy_width, max_mask_order, original_foreground);
    
  vector<float> fg_base_comp;
  get_base_comp(original_foreground, fg_base_comp);
  const size_t fg_length = effective_sequence_length(original_foreground);
	 
  if (!bgfilename.empty()) {
    // read bg sequences
    read_fasta_file(bgfilename, bgnames, original_background);
    
    mask(decoy_width, max_mask_order, original_background);
    
    vector<float> bg_base_comp;
    get_base_comp(original_background, bg_base_comp);
    const size_t bg_length = effective_sequence_length(original_background);
    for (size_t i = 0; i < alphabet_size; ++i)
      base_comp.push_back((fg_length*fg_base_comp[i] +
			   bg_length*bg_base_comp[i])/(fg_length + bg_length));
    
    foreground = original_foreground;
    if (!single_strand)
      for (size_t i = 0; i < original_foreground.size(); i++)
	foreground[i] += (string("N") + revcomp(foreground[i]));
    
    background = original_background;
    if (!single_strand)
      for (size_t i = 0; i < background.size(); i++)
	background[i] += (string("N") + revcomp(background[i]));
    
    length_ratio = static_cast<float>(fg_length)/bg_length;
  }
  else {
    base_comp = fg_base_comp;
    foreground = original_foreground;
    if (!single_strand)
      for (size_t i = 0; i < foreground.size(); ++i)
	foreground[i] += "N" + revcomp(foreground[i]);
  }
} // END preprocess_sequences()



void
validate_parameters(size_t &motif_width, float &bits) {
  
  struct ParamSet {
    int width;
    float bits;
  };
  
  static const size_t max_motif_width = 17;
  const struct ParamSet param_set[] = {
    {  0, 2.000 },
    {  1, 2.000 },
    {  2, 2.000 },
    {  3, 2.000 },
    {  4, 2.000 },
    {  5, 2.000 },
    {  6, 1.900 },
    {  7, 1.900 },
    {  8, 1.800 },
    {  9, 1.675 },
    { 10, 1.600 },
    { 11, 1.550 },
    { 12, 1.500 },
    { 13, 1.450 },
    { 14, 1.400 },
    { 15, 1.350 },
    { 16, 1.300 },
    { 17, 1.250 }
  };
  
  // check to make sure the motif width is reasonable
  if (motif_width > max_motif_width || motif_width == 0)
    throw SMITHLABException("motif width must be at least 1 and less than " +
			    toa(max_motif_width));
  
  if (bits == numeric_limits<float>::max())
    bits = param_set[motif_width].bits;
  
} // END validate_parameters()



int 
main(int argc, const char **argv) {

  try {
    
    // Not a parameter:
    static const char *motif_prep_progress_prefix = "preparing motifs ";
    
    string bgfilename;  // background sequences file
    string outfilename;     // output file

    string accession_prefix("DME"); // names of all
    // output motifs
    // start with this

    static bool ZOOPS = 0;
    static bool TCM = 0;

    static float bits =              // minimum average information content
      numeric_limits<float>::max();  // per position for matrices in
				     // search space
    static bool singlestrand = false; // Indicates if both strands are
				     // to be used
    static float granularity = 1.0;  // granularity of motifs in
				     // search space
    static float refine_granularity = 0.25; // minimum granularity to use
                                            // when refining motifs
    static float correction = 1e-10; // correction value to be added
				     // to each matrix entry
    static float ratio_adjust = 1.0;
    static size_t motif_width = 8;   // minimum width of the motifs to discover
    static size_t outputs = 1;       // number of outputs to print
    static size_t n_changes = 1;
    static size_t n_iterations = 100;
    static const float required_improvement = 1e-6;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "discriminating matrix "
			   "enumeration for motif discovery",
			   "<fasta-sequences>");
    opt_parse.add_opt("zoops", 'z', "use the ZOOPS model (default: hybrid)", 
		      false, ZOOPS);
    opt_parse.add_opt("tcm", 't', "use the TCM model (default: hybrid)", 
		      false, TCM);
    opt_parse.add_opt("background", 'b', "background sequence file (FASTA format)",
		      false, bgfilename);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
		      false, outfilename);
    opt_parse.add_opt("number", 'n', "number of motifs to produce.",
		      false, outputs);
    opt_parse.add_opt("prefix", 'p', "motif accession prefix",
		      false, accession_prefix);
    opt_parse.add_opt("width", 'w', "motif width",
		      false, motif_width);
    opt_parse.add_opt("bits", 'i', "min bits per column (default depends on width)",
		      false, bits);
    // opt_parse.add_opt("granularity", 'g', "see documentation (advanced option)",
    // 		      false, granularity);
    // opt_parse.add_opt("correction", 'c', "correction for 0 in matrices",
    // 		      false, correction);
    // opt_parse.add_opt("refine", 'r', "refinement granularity (default depends on width)",
    // 		      false, refine_granularity);
    // opt_parse.add_opt("adjust", 'a', "adjust contribution of fg and bg",
    // 		      false, ratio_adjust);
    // opt_parse.add_opt("changes", 'c', "changes per refinement",
    // 		      false, n_changes);
    // opt_parse.add_opt("iterations", 'I', "number of refinement iterations",
    // 		      false, n_iterations);
    opt_parse.add_opt("single-strand", '\0', "search only one strand",
		      false, singlestrand);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string fgfilename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    
    // make sure the parameters are sensible
    validate_parameters(motif_width, bits);
    
    vector<string> fgnames, original_foreground, foreground;
    vector<string> bgnames, original_background, background;
    vector<float> base_comp;
    float fg_bg_ratio = 1;

    vector<Matrix> seeds;
    
    /* ZOOPS: Zero or one occurrence per sequence */
    if (ZOOPS) {
      if (TCM) {
	throw SMITHLABException("ZOOPS and TCM options are incompatible");
      }
      preprocess_sequences_zoops(singlestrand,
				 fgfilename, bgfilename,
				 fgnames, original_foreground, foreground,
				 bgnames, original_background, background, 
				 base_comp, fg_bg_ratio);
      get_seeds_zoops(singlestrand,
		      foreground, background, base_comp,
		      motif_width, outputs, granularity, bits,
		      correction, ratio_adjust, seeds);
      
      refine_matrices_zoops(singlestrand,
			    foreground, background, motif_width, 
			    refine_granularity, bits, base_comp, 
			    n_changes, n_iterations, required_improvement, 
			    correction, ratio_adjust, seeds);
      
    }
    /* TCM: Two compenent mixture (0-to-many occurrences per sequence */
    else if (TCM) {
      preprocess_sequences_tcm(singlestrand,
			       fgfilename, bgfilename,
			       fgnames, original_foreground, foreground,
			       bgnames, original_background, background, 
			       base_comp, fg_bg_ratio);

      get_seeds_tcm(singlestrand,
		    foreground, background, base_comp,
		    motif_width, outputs, granularity,
		    bits, correction, ratio_adjust, seeds);
      
      refine_matrices_tcm(singlestrand,
			  foreground, background, motif_width, 
			  refine_granularity, bits, base_comp, 
			  n_changes, n_iterations, required_improvement, 
			  correction, ratio_adjust, seeds);
    }
    else { // HYBRID
      preprocess_sequences_zoops(singlestrand,
				 fgfilename, bgfilename,
				 fgnames, original_foreground, foreground,
				 bgnames, original_background, background, 
				 base_comp, fg_bg_ratio);
      
      get_seeds_tcm(singlestrand,
		    foreground, background, base_comp,
		    motif_width, outputs, granularity,
		    bits, correction, ratio_adjust, seeds);
      
      refine_matrices_zoops(singlestrand,
			    foreground, background, motif_width, 
			    refine_granularity, bits, base_comp, 
			    n_changes, n_iterations, required_improvement, 
			    correction, ratio_adjust, seeds);
      


    }

    vector<Motif> motifs;
    
    /* Turn matrices into motifs when ZOOPS model is used */
    if (!TCM) {
      for (size_t i = 0; i < seeds.size(); ++i) {
	if (VERBOSE)
	  cerr << "\r" << motif_prep_progress_prefix
	       << i + 1 << "/" << seeds.size();
	motifs.push_back(prepare_motif_zoops(accession_prefix + toa(i + 1),
					     seeds[i], base_comp, correction,
					     original_foreground, fgnames,
					     original_background, bgnames,
					     singlestrand, fg_bg_ratio));
      }
    }

    /* Turn matrices into motifs when TCM model is used */
    else {
      for (size_t i = 0; i < seeds.size(); ++i) {
	if (VERBOSE)
	  cerr << "\r" << motif_prep_progress_prefix
	       << i + 1 << "/" << seeds.size();
	motifs.push_back(prepare_motif_tcm(accession_prefix + toa(i + 1),
					   seeds[i], base_comp, correction,
					   original_foreground, fgnames,
					   original_background, bgnames,
					   singlestrand, fg_bg_ratio));
      }
    }
    
    if (VERBOSE)
      cerr << endl << "ranking motifs";
    std::stable_sort(motifs.begin(), motifs.end(), 
		     PatternOrder("SCORE", true, true));
    if (VERBOSE)
      cerr << ". done." << endl;
    
    // prepare the output stream for motifs

    std::ofstream of;
    if (!outfilename.empty()) of.open(outfilename.c_str());
    std::ostream out(outfilename.empty() ? cout.rdbuf() : of.rdbuf());
    copy(motifs.begin(), motifs.end(), ostream_iterator<Motif>(out, "\n"));
  }
  catch (SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
