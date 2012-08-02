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

#ifndef SCORINGMATRIX_HPP
#define SCORINGMATRIX_HPP

/**
  \file ScoringMatrix.hpp

  \brief This header file contains class definitions for ScoringMatrix
  and the related StadenPValue class, along with associated non-member
  function declarations.
*/

#include "cread.hpp"
#include "Matrix.hpp"

/**
   \brief Class to represent log-ratio scoring matrices for motifs.
   \class ScoringMatrix
*/
class ScoringMatrix {
public:
  
  /**
    \brief Default constructor, width-zero scoring matrix
  */
  ScoringMatrix() : matrix(0), width(0) {}
  
  explicit ScoringMatrix(float **m, const size_t w);
  
  /**
    \brief Minimal constructor with non-empty initializations
    
    \param mat Matrix object to convert into the ScoringMatrix
    \param base_comp float vector of base probabilities
    \param correction float to add to matrix entries to prevent log of 0
  */
  ScoringMatrix(const Matrix& mat,
		const std::vector<float> &base_comp,
		float correction = std::numeric_limits<float>::min());
  
  /**
    \brief Minimal constructor with non-empty initializations

    \param mat Matrix object to convert into the ScoringMatrix
    \param base_comp float array of base probabilities
    \param correction float to add to matrix entries to prevent log of 0
  */
  ScoringMatrix(const Matrix& mat, 
		const float base_comp[],
		float correction = std::numeric_limits<float>::min());
  
  /**
    \brief copy constructor; needed because of dynamic allocation
    \param other the other scoring matrix to copy
  */
  ScoringMatrix(const ScoringMatrix &other);
  

  /**
    \brief assignment operator; needed because of dynamic allocation

    \param rhs the scoring matrix being assigned

    \return reference to a scoring matrix
  */
  ScoringMatrix& operator=(const ScoringMatrix &rhs);
  
  /**
    \brief destructor; needed because of dynamic allocation
  */
  ~ScoringMatrix() {delete_matrix(matrix, width);}
  
  /**
    \brief swap method; defined for efficiency

    \param sm the other scoring matrix with which to swap internals
  */
  void swap(ScoringMatrix &sm);

  typedef float** pointer;
  typedef const float*const* const_pointer;
  typedef float*& reference;
  typedef const float*const& const_reference;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  
  /// Iterator to start of scoring matrix columns.
  iterator begin() {return matrix;}

  /// Iterator to end of scoring matrix columns.
  iterator end() {return matrix + width;}

  /// Constant iterator to start of scoring matrix columns.
  const_iterator begin() const {return matrix;}

  /// Constant iterator to end of scoring matrix columns.
  const_iterator end() const {return matrix + width;}

  /// Get the n-th column of the scoring matrix
  reference operator[](size_t n) {return *(begin() + n);}

  /// Get the n-th column of the scoring matrix (as const. ref.)
  const_reference operator[](size_t n) const {return *(begin() + n);}
  
  /// Get string representation of the scoring matrix
  std::string tostring() const;
  
  // accessors

  /// Get the number of columns in the scoring matrix
  size_t size() const {return width;}

  /// Get the number of columns in the scoring matrix
  size_t get_width() const {return width;}

  /// Get the reverse complement of the scoring matrix
  ScoringMatrix revcomp() const;

  /**
     \brief convert a score to a functional depth
  */
  float functional_depth(const float score) const;
  /**
     \brief convert a functional depth to a score
  */
  float functional_depth_to_score(const float fd) const;
  
  static ScoringMatrix StormoScoringMatrix(const Matrix& mat, 
					   const float *base_comp);
  static ScoringMatrix StormoScoringMatrix(const Matrix& mat,
					   const std::vector<float> &basecomp);
  
private:
  float **matrix;
  size_t width;
  
  // Functions to delete matrices (and columns)
  template <class T> static void delete_column(T* p) {delete[] p;}
  template <class T> static void delete_matrix(T **p, int n) {
    std::for_each(p, p + n, delete_column<T>);
    delete[] p;
  }
  
  static float epsilon() {return 0.005;}
  static float default_correction() {return 0.0000000001;}
};

std::ostream& 
operator<<(std::ostream& s, const ScoringMatrix &sm);

/**
   \brief Helper class to organize and speed calculation of Staden
   p-values for scoring matrix matches.
   
   \class StadenPValue
   
   This class is used to implement the calculation of p-values for
   scoring matrix matches as described originally by Staden (ref)
*/
class StadenPValue {
public:
  StadenPValue(const std::vector<float>& bc, const float sc) :
    base_comp(bc), scale(sc) {
    used = std::vector<size_t>(default_max_score + 1);
    prev_used = std::vector<size_t>(default_max_score + 1);
    counts = std::vector<double>(default_max_score + 1);
    prev_counts = std::vector<double>(default_max_score + 1);
  }
  float get_pvalue(ScoringMatrix sm, const float score) const;
  float get_score(ScoringMatrix sm, const float pvalue) const;
  void get_pvalues(ScoringMatrix sm, 
		   const std::vector<float> &scores,
		   std::vector<float> &pvalues) const;
  void get_scores(ScoringMatrix sm, 
		  const std::vector<float> &pvalue,
		  std::vector<float> &scores) const;
  
private:
  static const size_t default_max_score = 10000;
  
  mutable std::vector<size_t> used;
  mutable std::vector<size_t> prev_used;
  mutable std::vector<double> counts;
  mutable std::vector<double> prev_counts;
  
  std::vector<float> base_comp;
  float scale;
  
  float normalize_matrix(ScoringMatrix& sm) const;
};

#endif
