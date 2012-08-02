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

#ifndef FASTA_FILE
#define FASTA_FILE

#include "cread.hpp"

//!  A class for reading large FASTA format files one chunk at a time.
/*!  This class is for those really big FASTA files (e.g. chromosomes)
that may be better read into memory one chunk at a time.
*/
class BigFastaFile {
  static const size_t input_buffer_size = 1000000;
  
  struct SequenceInfo {
    std::vector<size_t> newlines;
    std::string name;
    size_t start;
    size_t end;
    explicit SequenceInfo(std::string n, size_t s, size_t e, 
			  std::vector<size_t> &nl) :
      newlines(nl), name(n), start(s), end(e) {}
    bool operator<(size_t val) const {
      return end < val;
    }
    size_t count_newlines(size_t range_start, size_t range_end) const;
    std::string tostring() const;
  };
  
  typedef std::vector<SequenceInfo>::const_iterator SeqInfoIter;
  
  std::vector<SequenceInfo> SequenceTable;
  std::string filename;
  size_t filesize;
  SeqInfoIter find_containing_sequence(size_t position) const;
  void make_index(std::string filename);
  
public:
  explicit BigFastaFile(std::string fn);
  std::vector<std::string> get_names(const size_t chunk_start,
				     const size_t chunk_end) const;
  std::vector<std::string> get_sequences(const size_t chunk_start,
					 const size_t chunk_size) const;
  size_t get_first_sequence_offset(const size_t chunk_start) const;
  size_t get_filesize() const {return filesize;}
  std::string tostring() const;
};

class FastaFile {
  // PRIVATE STATIC MEMBERS
  const static size_t default_input_buffer_size = 1000000;
  
  // INSTANCE VARIABLES
  std::string filename;
  size_t input_buffer_size;

  // all this stuff is mutable so if can be initialized only if
  // and when needed:
  mutable size_t filesize;
  mutable std::vector<std::string> names;
  mutable std::vector<std::string> sequences;

  // This is used to keep track of where the sequences (and their
  // names) begin or end when we are dealing with huge files that are
  // not resident in memory all at once.
  struct faindex {
    std::string name;
    size_t start;
    size_t end;
    explicit faindex(std::string n, size_t s, size_t e) : 
      name(n), start(s), end(e) {}
    bool operator<(size_t val) const {
      return start < val;
    }
  };
  // mutable for just-in-time initialization
  mutable std::vector<faindex> index;
  
  // PRIVATE MEMBER FUNCTIONS
  void read(void) const;  // these can be const because
  void clean() const;
  void toupper() const;
public:
  
  // CONSTRUCTORS
  FastaFile() : input_buffer_size(default_input_buffer_size),
		filesize(std::numeric_limits<size_t>::max()) {}
  FastaFile(std::string &fn) : filename(fn),
			       input_buffer_size(default_input_buffer_size),
			       filesize(std::numeric_limits<size_t>::max()) {}
  FastaFile(const char *fn) : filename(fn),
			      input_buffer_size(default_input_buffer_size),
			      filesize(std::numeric_limits<size_t>::max()) {}
  // accessors
  size_t get_filesize() const;
  size_t size() const {return sequences.size();}
  std::vector<std::string> get_sequences() const {
    if (sequences.empty()) read();
    return sequences;
  }
  std::vector<std::string> get_names() const {
    if (names.empty()) read();
    return names;
  }
  std::vector<std::string> get_ids() const;
  std::vector<std::string> get_descriptions() const;
  
  // static class functions
  static std::vector<std::string> read_seqs_dir(const std::string& dir, 
						const std::string& suff);
  static void base_comp_from_file(std::string fn, float *bc, 
				  size_t buffer_size = default_input_buffer_size);
  static void base_comp_from_files(const std::vector<std::string>& fns, 
				   std::vector<float>& bc,
				   size_t buffer_size = default_input_buffer_size);
  static void clean(std::vector<std::string>&);
  static void toupper(std::vector<std::string>&);
};

class FastaFileException : public CREADException {
public:
  FastaFileException(std::string m = "") : CREADException(m) {}
}; 

#endif
