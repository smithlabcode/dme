/*
 * Copyright (C) 2006 Cold Spring Harbor Laboratory
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

#include "Alphabet.hpp"
#include "cread.hpp"

using std::string;
using std::vector;
using std::transform;

// private static members
static const int b2i_size = 20;
static const int b2i[] = {
//A, b, C, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, T
  0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3
};
static const int b2i_rc[] = {
//A, b, C, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, T
  3,-1, 2,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0
};
static const int b2c_size = 20;
static const char b2c[] = {
 //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
  'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
};
static const char *i2b = "ACGTN";


int
base2int(char b) {
  b = std::toupper(b);
  if (b - 'A' >= 0 && b - 'A' < b2i_size) 
    return b2i[b - 'A'];
  else
    return -1;
}


int
base2int_rc(char b) {
  b = std::toupper(b);
  if (b - 'A' >= 0 && b - 'A' < b2i_size) 
    return b2i_rc[b - 'A'];
  else
    return -1;
}


char
int2base(int i) {
  if (i < static_cast<int>(alphabet_size) && i >= 0)
    return i2b[i];
  else return i2b[alphabet_size];
}


bool
valid_base_id(int c) {
  return (c < static_cast<int>(alphabet_size) && c >= 0);
}


char
complement(int i) {
  if (i - 'A' >= 0 && i - 'A' < b2c_size)
    return b2c[i - 'A'];
  else return 'N';
}


bool
valid_base(char c) {
  char i = std::toupper(c);
  return (i == 'A' || i == 'C' || i == 'G' || i == 'T');
}

bool
isgap(char c) {
  return c == '-';
}

void
get_base_comp(const vector<string>& sequences, float *base_comp) {
  std::fill(base_comp, base_comp + alphabet_size, 0.0);
  int total = 0;
  for (vector<string>::const_iterator i = sequences.begin();
       i != sequences.end(); ++i)
    for (string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (valid_base(*j)) {
	base_comp[base2int(*j)]++;
	total++;
      }
  transform(base_comp, base_comp + alphabet_size, base_comp, 
		 bind2nd(std::divides<float>(), total));
}


void
get_base_comp(const vector<string>& sequences, vector<float>& base_comp) {
  vector<size_t> count(alphabet_size, 0);
  for (vector<string>::const_iterator i = sequences.begin();
       i != sequences.end(); ++i)
    for (string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (valid_base(*j)) {
	count[base2int(*j)]++;
      }
  const float total = std::accumulate(count.begin(), count.end(), 0.0);
  base_comp.clear();
  transform(count.begin(), count.end(), back_inserter(base_comp), 
		 std::bind2nd(std::divides<float>(), total));
}


void
get_base_comp(const vector<string>& sequences, vector<double>& base_comp) {
  vector<size_t> count(alphabet_size, 0);
  for (vector<string>::const_iterator i = sequences.begin();
       i != sequences.end(); ++i)
    for (string::const_iterator j = i->begin(); j != i->end(); ++j)
      if (valid_base(*j)) {
	count[base2int(*j)]++;
      }
  const double total = std::accumulate(count.begin(), count.end(), 0.0);
  base_comp.clear();
  transform(count.begin(), count.end(), back_inserter(base_comp), 
		 std::bind2nd(std::divides<double>(), total));
}


string
reverse_complement(const string& s) {
  string r;
  transform(s.begin(), s.end(), back_inserter(r), complement);
  reverse(r.begin(), r.end());
  return r;
}


string
revcomp(const string& s) {
  string r;
  transform(s.begin(), s.end(), back_inserter(r), complement);
  reverse(r.begin(), r.end());
  return r;
}


size_t
count_valid_bases(const string& s) {
  return count_if(s.begin(), s.end(), &valid_base);
}


size_t
count_valid_bases(const vector<string>& s) {
  size_t n_valid = 0;
  for (vector<string>::const_iterator i = s.begin(); i != s.end(); ++i)
    n_valid += count_valid_bases(*i);
  return n_valid;
}


size_t
kmer_counts(const vector<string> &seqs, 
	    vector<size_t> &counts, size_t k) {
  counts.clear();
  size_t nwords = static_cast<size_t>(pow(static_cast<float>(alphabet_size), 
					  static_cast<int>(k)));
  counts.resize(nwords, 0);
  size_t total = 0;
  for (size_t i = 0; i < seqs.size(); ++i) {
    char seq[seqs[i].length() + 1];
    seq[seqs[i].length()] = '\0';
    copy(seqs[i].begin(), seqs[i].end(), seq);
    for (size_t j = 0; j < seqs[i].length() - k + 1; ++j)
      if (std::count_if(seq + j, seq + j + k, &valid_base) == 
	  static_cast<int>(k)) {
	counts[mer2index(seq + j, k)]++;
	++total;
      }
  }
  return total;
}


size_t
mer2index(const char *s, size_t n) {
  size_t multiplier = 1, index = 0;
  do {
    --n;
    index += base2int(s[n])*multiplier;
    multiplier *= alphabet_size;
  } while (n > 0);
  return index;
}

size_t
mer2index_rc(const char *s, size_t n) {
  size_t multiplier = 1, index = 0;
  size_t k = 0;
  do {
    index += base2int_rc(s[k])*multiplier;
    multiplier *= alphabet_size;
  } while (++k < n);
  return index;
}
