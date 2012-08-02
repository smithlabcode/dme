/*
 * Copyright (C) 2007 Cold Spring Harbor Laboratory
 * Authors: Andrew D. Smith
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

#ifndef ALPHABET_HPP
#define ALPHABET_HPP

/*! \file Alphabet.hpp
    \brief Contains commonly used functions that operate on representations
	   of bases and sequences and their representations.
*/ 



#include <string>
#include <vector>

struct Alpha {
  static const size_t invalid = static_cast<size_t>(-1);
};

/// Returns the reverse complement of a string ACGTN are legal characters
std::string
reverse_complement(const std::string& seq);

/// Returns the reverse complement of a string ACGTN are legal characters
std::string 
revcomp(const std::string&);

/*!
  \brief Compute the base composition of a set of sequences.
  \param[in] sequences
      A set of sequences.
  \param[out] base_comp
      An array of size alphabet_size containing the frequencies
      for each nucleotide (sum to 1).
  \warning
      Space for base_comp should be allocated before calling this function.
      Only valid bases will be counted.
 */
void 
get_base_comp(const std::vector<std::string>& sequences, 
	      float* base_comp);
/*!
  \brief Compute the base composition of a set of sequences.
  \param[in] sequences
      A set of sequences.
  \param[out] base_comp
      The frequency vector for each nucleotide (sum to 1).
  \warning
      Only valid bases will be counted.
 */
void 
get_base_comp(const std::vector<std::string>& sequences,
	      std::vector<float>& bc);

void 
get_base_comp(const std::vector<std::string>& sequences,
	      std::vector<double>& bc);

/// Returns the number of valid bases in the input string
size_t 
count_valid_bases(const std::string &s);

/// Returns the number of valid bases in the string vector
size_t 
count_valid_bases(const std::vector<std::string> &s);

/// Convert from character to integer
int
base2int(char b);

int
base2int_rc(char b);


/// Convert from integer to character
char
int2base(int i);

/// Returns true if c is a valid base
bool
valid_base(char c);

/// Returns true if c is a valid integer representation for a base
bool
valid_base_id(int c);

bool
isgap(char c);

/// Returns the identity of the complement of the base corresponding to i
char
complement(int i);

/*!
  \brief Convert the n-length word embedded in array s to an unsigned integer.
  \parm[in] s character array representation for a word
  \parm[in] n length of the word
  \return unsigned integer representation for a word
 */
size_t
mer2index(const char *s, size_t n);

size_t
mer2index_rc(const char *s, size_t n);

/*!
  \brief Count the number of occurrences of each word of length k in a vector
	 of sequences.
  \parm[in] seqs A vector of DNA of sequences.
  \param[in] k Word size.
  \param[out] counts
      A Vector of counts for unsigned integer representations for each word
      of length k.
  \return Total number of k-length words made of valid bases.
 */
size_t
kmer_counts(const std::vector<std::string> &seqs, 
	    std::vector<size_t> &counts, 
	    size_t k);
#endif
