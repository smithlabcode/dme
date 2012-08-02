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

#include "Matrix.hpp"

using std::accumulate;
using std::fill;
using std::string;
using std::vector;
using std::transform;
using std::divides;
using std::bind2nd;
using std::reverse_copy;
using std::greater;
using std::endl;

std::ostream& 
operator<<(std::ostream& s, const Matrix& mat) {
  return s << mat.tostring();
}

Matrix::Matrix(float** m, size_t w) : width(w) {
  if (m) {
    matrix = new float*[width];
    for (size_t i = 0; i < width; i++) {
      matrix[i] = new float[alphabet_size];
      std::copy(m[i], m[i] + alphabet_size, matrix[i]);
    }
  }
  else matrix = 0;
}

Matrix::Matrix(const Matrix& M) : width(M.width) {
  if (M.matrix) {
    matrix = new float *[width];
    for (size_t i = 0; i < width; i++) {
      matrix[i] = new float[alphabet_size];
      std::copy(M.matrix[i], M.matrix[i] + alphabet_size, matrix[i]);
    }
  }
  else matrix = 0;
}

Matrix::Matrix(vector<string> &s) {
  vector<float*> columns;
  for (size_t i = 0; i < s.size(); ++i) {
    if (isdigit(s[i][0])) {
      float* column = new float[alphabet_size];
      if (sscanf(s[i].c_str(), "%*d %f %f %f %f %*c", 
		 column, &column[1], &column[2], &column[3]) - 
	  alphabet_size != 0) {
	delete column;
	throw MatrixException("Bad Matrix column line: " + s[i]);
      }
      columns.push_back(column);
    }
  }
  width = columns.size();
  if (width > 0) {
    matrix = new float*[width];
    for (size_t i = 0; i < width; ++i) {
      matrix[i] = new float[alphabet_size];
      std::copy(columns[i], columns[i] + alphabet_size, matrix[i]);
    }
    for (vector<float*>::iterator i = columns.begin(); 
	 i != columns.end(); ++i)
      delete[] *i;
  }
  else matrix = 0;
}

void
Matrix::swap(Matrix &other) {
  std::swap(matrix, other.matrix);
  std::swap(width, other.width);
}

Matrix& 
Matrix::operator=(const Matrix& M) {
  Matrix tmp(M);
  swap(tmp);
  return *this;
}


Matrix::~Matrix() {
  if (matrix) {
    std::for_each(matrix, matrix + width, cread::delete_ptr<float>);
    delete[] matrix;
  }
}


string 
Matrix::tostring() const {
  std::ostringstream s;
  if (matrix != static_cast<float **>(0)) {
    s << "P0       A       C       G       T\n";
    for (size_t i = 0; i < width; i++) {
      s.width(2);
      s.fill('0');
      s << i + 1;
      for (size_t j = 0; j < alphabet_size; j++) {
	s.width(8); 
	s.fill(' '); 
	s.setf(std::ios_base::right);
	s.setf(std::ios_base::fixed, std::ios_base::floatfield);
	// NEED TO USE SOME is_count_mat() TYPE THING HERE
	if (matrix[i][j] - floor(matrix[i][j]) > 0.00001)
	  s.precision(3);
	else s.precision(0);
	s << matrix[i][j];
      }
      if (i != width - 1) s << endl;
    }
  }
  return s.str();
}

Matrix
Matrix::combine(const Matrix &matrix1, const Matrix &matrix2, int offset) {
  Matrix newmat;
  newmat.width = std::max(matrix1.width,
			  matrix2.width + offset);
  newmat.matrix = new float*[newmat.width];
  for (size_t i = 0; i < newmat.width; ++i) {
    newmat.matrix[i] = new float[alphabet_size];
    fill(newmat.matrix[i], newmat.matrix[i] + alphabet_size, 0);
  }
  for (size_t i = 0; i < matrix1.width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      newmat.matrix[i][j] += matrix1.matrix[i][j];
  for (size_t i = 0; i < matrix2.width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      newmat.matrix[i + offset][j] += matrix2.matrix[i][j];
  for (size_t i = 0; i < newmat.width; ++i) {
    const float total = accumulate(newmat.matrix[i],
				   newmat.matrix[i] + alphabet_size, 0.0);
    transform(newmat.matrix[i], newmat.matrix[i] + alphabet_size,
	      newmat.matrix[i], bind2nd(divides<float>(), total));
  }
  return newmat;
}

Matrix
Matrix::freqmat() const {
  float **r = new float*[width];
  for (size_t i = 0; i < width; i++) {
    float total = std::accumulate(matrix[i], matrix[i] + alphabet_size, 0.0);
    r[i] = new float[alphabet_size];
    transform(matrix[i], matrix[i] + alphabet_size, r[i],
	      std::bind2nd(std::divides<float>(), total));
  }
  return Matrix(r, width);
}

Matrix
Matrix::corrected_freqmat(float base_comp[]) const {
  if (!is_count_mat()) 
    return Matrix(*this);
  float **r = new float*[width];
  for (size_t i = 0; i < width; i++) {
    float total = std::accumulate(matrix[i],
				  matrix[i] + alphabet_size, 0.0) + 1;
    r[i] = new float[alphabet_size];
    transform(matrix[i], matrix[i] + alphabet_size, base_comp, r[i],
	      std::plus<float>());
    transform(r[i], r[i] + alphabet_size, r[i],
	      std::bind2nd(std::divides<float>(), total));
  }
  return Matrix(r, width);
}

float
Matrix::GetInformationContent(float f[]) const {
  float information = 0;
  for (size_t i = 0; i < width; i++){
    float sum = std::accumulate(matrix[i], matrix[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; j++)
      if (matrix[i][j] > 0)
        information += matrix[i][j]/sum*(cread::log2(matrix[i][j]/sum) -
					 cread::log2(f[j]));
  }
  return information;
}

float
Matrix::GetInformationContent(vector<float> &f) const {
  float information = 0;
  for (size_t i = 0; i < width; i++){
    float sum = std::accumulate(matrix[i], matrix[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; j++)
      if (matrix[i][j] > 0)
        information += matrix[i][j]/sum*(cread::log2(matrix[i][j]/sum) -
					 cread::log2(f[j]));
  }
  return information;
}




float
Matrix::info(const float f[]) const {
  float information = 0;
  for (size_t i = 0; i < width; i++){
    float sum = std::accumulate(matrix[i], matrix[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; j++)
      if (matrix[i][j] > 0)
        information += matrix[i][j]/sum*(cread::log2(matrix[i][j]/sum) -
					 cread::log2(f[j]));
  }
  return information;
}

float
Matrix::info(const vector<float> &f) const {
  float information = 0;
  for (size_t i = 0; i < width; i++){
    float sum = std::accumulate(matrix[i], matrix[i] + alphabet_size, 0.0);
    for (size_t j = 0; j < alphabet_size; j++)
      if (matrix[i][j] > 0)
        information += matrix[i][j]/sum*(cread::log2(matrix[i][j]/sum) -
					 cread::log2(f[j]));
  }
  return information;
}

Matrix
Matrix::revcomp() const {
  Matrix m;
  m.width = width;
  m.matrix = new float *[width];
  for (size_t i = 0; i < width; ++i) {
    m.matrix[i] = new float[alphabet_size];
    reverse_copy(matrix[width - i - 1],
		 matrix[width - i - 1] + alphabet_size, m.matrix[i]);
  }
  return m;
}

bool
Matrix::is_count_mat() const {
  // TODO: this needs to be much more robust, or maybe we need to have
  // a 3 way split: matrices can be freq or count, or just plain wrong
  for (Matrix::const_iterator i = begin(); i != end(); ++i)
    if (find_if(*i, *i + alphabet_size, bind2nd(greater<float>(), 1)) 
	!= *i + alphabet_size)
      return true;
  return false;
}

void 
Matrix::base_comp_pseudocount(const std::vector<float> &base_comp,
			      const float weight) {
  for (size_t i = 0; i < width; ++i)
    for (size_t j = 0; j < alphabet_size; ++j)
      matrix[i][j] += base_comp[j]*weight;
}

void 
Matrix::Laplace_pseudocount() {
  for (size_t i = 0; i < width; ++i)
    transform(matrix[i], matrix[i] + alphabet_size, matrix[i],
	      std::bind2nd(std::plus<float>(), 1.0));
}
