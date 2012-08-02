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

#ifndef CREAD_H
#define CREAD_H
/*!
  \file cread.hpp
  \brief Contains the cread namespace, constants that are used project wide,
	 and commonly used functions for string parsing and manipulation.
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <set>
#include <map>
#include <ctime>
#include <cctype>
#include <functional>
#include <utility>
#include <limits>
#include <queue>
#include <list>
#include <numeric>
#include <cassert>

// TODO: put this into the cread namespace:
/// The number of DNA bases
const size_t alphabet_size = 4;

/*!
   \namespace cread
   \brief Contains commonly used functions.
 */
namespace cread {
  /// Split a string into a vector of strings using a character array delimiter
  std::vector<std::string> split(std::string, const char *, 
				 bool get_empty_fields = false);
  /**
     \brief Return a vector of strings resulting from splitting
     the input string to_split using white space delimiters.
   */
  std::vector<std::string> split_whitespace_quoted(std::string to_split);
  /*! \brief Store a vector of strings in v resulting from splitting
   	     the input using white space delimiters.
   */
  void split_whitespace(const std::string& s, std::vector<std::string> &v);

  /// Produce a string with striped leading and trailing white spaces
  std::string strip(const std::string& s);

  /// Translate to a string
  template <class T> std::string toa(T);

  /// Log base 2
  template <class T> T log2(T t) {return log(t)/log(2.0);}

  /// Delete a pointer (used with STL algorithms to process containers)
  template <class T> T* delete_ptr(T* p) {delete p; return 0;}

}

/**
   \exception CREADException
   \brief Provides a base class for all cread exceptions.
*/
class CREADException : public std::exception {
public:
  CREADException() throw() {}
  CREADException(std::string m) throw(): message(m) {}
  virtual ~CREADException() throw() {}
  const char* what() const throw() { return message.c_str(); }
protected:
  std::string message;
};

template <class T> std::string 
cread::toa(T t) {
  std::ostringstream s;
  s << t;
  return s.str();
}

inline std::vector<std::string> 
cread::split(std::string s, const char *delim, bool get_empty_fields) {
  std::vector <std::string> parts;
  size_t i = 0, j = 0, dlen = strlen(delim);
  while (i < s.length()) {
    bool prev_in_delim = false;
    if (i > 0)
      for (size_t k = 0; k < dlen; ++k) 
	if (s[i-1] == delim[k]) 
	  prev_in_delim = true;
    bool curr_in_delim = false;
    for (size_t k = 0; k < dlen; ++k) 
      if (s[i] == delim[k]) 
	curr_in_delim = true;
    if (curr_in_delim && (get_empty_fields || !prev_in_delim)) {
      parts.push_back(s.substr(j, i - j));
      j = i + 1;
    }
    i++;
  }
  if (i > j || (get_empty_fields && i == j)) 
    parts.push_back(s.substr(j, i - j));
  return parts;
}

inline std::string
cread::strip(const std::string& s) {
  const size_t len = s.length();
  size_t i = 0;
  while (i < len && isspace(s[i])) {
    i++;
  }
  size_t j = len;
  do {
    j--;
  } while (j >= i && isspace(s[j]));
  j++;
  if (i == 0 && j == len)
    return s;
  else return s.substr(i, j - i);
}

inline void
cread::split_whitespace(const std::string& s, std::vector<std::string> &v) {
  size_t i = 0, len = s.length();
  while (i < len) {
    while (i < len && isspace(s[i])) ++i;
    size_t j = i;
    while (i < len && !isspace(s[i])) ++i;
    if (j < i)
      v.push_back(s.substr(j, i - j));
  }
}

inline std::vector<std::string>
cread::split_whitespace_quoted(std::string to_split) {
  static const char *non_word_chars = " \t\"'";
  to_split = cread::strip(to_split);
  
  std::vector<std::string> words;
  size_t start_pos = 0, end_pos = 0;
  const size_t length_of_to_split = to_split.length();
  while (start_pos < length_of_to_split) {
    /** find next position that is not a word character */
    end_pos = to_split.find_first_of(non_word_chars, end_pos);
    if (end_pos > to_split.length()) { /** If we hit the end: done */
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      break;
    }
    /** unescaped, unquoted white space: definitely a word delimiter */
    if (to_split[end_pos] == ' ' || to_split[end_pos] == '\t') { 
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      end_pos = to_split.find_first_not_of(" \t", end_pos);
      start_pos = end_pos;
    }
    /** preserve whatever is being escaped; will become part of the
	current word */
    else if (to_split[end_pos] == '\\')
      end_pos = to_split.find_first_not_of(non_word_chars, end_pos + 2);
    else {
      const std::string current_delim = "\\" + to_split.substr(end_pos, 1);
      do { // slurp doubly- or singly-quoted string
	end_pos = to_split.find_first_of(current_delim, end_pos + 1);
	if (end_pos == std::string::npos) {
	  end_pos = length_of_to_split;
	  break;
	}
	if (to_split[end_pos] == '\\')
	  ++end_pos;
	else break;
      } while (true);
      ++end_pos;
    }
    if (end_pos >= length_of_to_split) {
      words.push_back(to_split.substr(start_pos,
				      end_pos - start_pos));
      break;
    }
  }
  return words;
}


#endif // CREAD_H
