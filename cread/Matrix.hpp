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

#ifndef MATRIX_HPP
#define MATRIX_HPP

/**
  \file Matrix.hpp

  \brief This header file contains class definitions for Matrix, along
  with associated non-member function declarations.
*/

#include "cread.hpp"

/**
   \brief Represents position-weight matrices for DNA sequence motifs.
   
   \class Matrix
*/
class Matrix {
public:
  /**
     \brief Basic constructor
     
     \param mat array of columns, each a float array of size
     equal to alphabet_size (4 for DNA).
     \param w the number of columns in the array.
  */
  Matrix(float **mat = 0, size_t w = 0);
  
  /*!
    \brief Constructor that parses string set representation.
    
    \param lines_from_file an ordered set of lines, one for each
    column of the matrix
  */
  explicit Matrix(std::vector<std::string>& lines_from_file);

  /// Copy constructor
  Matrix(const Matrix& original);

  /// Assignment operator
  Matrix& operator=(const Matrix& lhs);

  /// Destructor
  ~Matrix();

  /// Swap method
  void swap(Matrix &other);

  typedef float** pointer;
  typedef const float *const * const_pointer;
  typedef float*& reference;
  typedef const float *const & const_reference;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  
  /// Iterator to start of matrix columns.
  iterator begin() {return matrix;}

  /// Iterator to end of matrix columns.
  iterator end() {return matrix + width;}

  /// Constant iterator to start of matrix columns.
  const_iterator begin() const {return matrix;}

  /// Constant iterator to end of matrix columns.
  const_iterator end() const {return matrix + width;}

  /// Get the n-th column of the matrix
  reference operator[](size_t n) {return *(begin() + n);}

  /// Get the n-th column of the matrix (const. ref.)
  const_reference operator[](size_t n) const {return *(begin() + n);}

  /// Get string representation of the matrix
  std::string tostring() const;
  
  /// Get the number of columns in the matrix
  size_t get_width() const {return width;}
  
  // mutators
  /*!
    \brief Add a scaled pseudocount to each matrix entry based on
    specified base composition.
    
    The pseudocount at an entry will be determined by the base
    composition (specific to each base), and the weight given to
    pseudocounts in general relative to real counts.
    
    \param base_comp frequencies of bases
    \param weight value multiplied with the base frequency to obtain
    the pseudocount
  */
  void base_comp_pseudocount(const std::vector<float> &base_comp,
			     float weight = 1.0);

  /*!
    \brief Add a unit value to each entry of the matrix
  */
  void Laplace_pseudocount();
  
  // TODO: Whoa!! very bad stuff below...
  const float **at(size_t n) const {return (const float **)(matrix + n);}
  
  /*!
    \brief Copy of matrix, with columns normalized to have unit sum.
    \return Copy of matrix, with columns normalized to have unit sum.
  */
  Matrix freqmat() const;

  /*!
    \brief Copy of matrix, with columns normalized to have unit sum,
    and a base composition pseudocount applied.
    \return Pseudocount-corrected and normalized copy of the matrix.
  */
  Matrix corrected_freqmat(float []) const;
  
  /// Get the reverse complement of the matrix
  Matrix revcomp() const;
  
  /// test if the matrix entries are counts or frequencies
  bool is_count_mat() const;
  
  /// Total information content of the matrix
  float GetInformationContent(float []) const;
  
  /// Total information content of the matrix
  float GetInformationContent(std::vector<float> &f) const;
  
  /// Total information content of the matrix
  float info(const std::vector<float> &) const;
  
  /// Total information content of the matrix
  float info(const float []) const;
  
  static Matrix combine(const Matrix& a, const Matrix& b, 
			int max_overhang = 0);
private:
  float **matrix;      
  size_t width;
  
  friend class MatCompMethods;
};

std::ostream& 
operator<<(std::ostream& s, const Matrix& mat);

/*!
 \exception MatrixException
 \brief Used for reporting errors when handeling matrices.
*/
class MatrixException : public CREADException  {
public:
  MatrixException(std::string m = "") throw() : CREADException(m) {}
}; 

#endif
