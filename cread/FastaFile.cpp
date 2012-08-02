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

#include "FastaFile.hpp"
#include "Alphabet.hpp"

#include <dirent.h>
#include <errno.h>

using std::string;
using std::vector;
using std::pair;
using std::cerr;
using std::ostringstream;
using std::endl;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////////// BigFastaFile stuff  ///////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

string 
BigFastaFile::SequenceInfo::tostring() const {
  ostringstream buff;
  buff << name << "\t" << start << "\t" 
       << end << "\t" << newlines.size();
  return buff.str();  
}

size_t 
BigFastaFile::SequenceInfo::count_newlines(size_t range_start,
					   size_t range_end) const {
  vector<size_t>::const_iterator lower =
    lower_bound(newlines.begin(), newlines.end(), range_start - start);
  vector<size_t>::const_iterator upper = 
    lower_bound(newlines.begin(), newlines.end(), range_end - start);
  return upper - lower;
}


void
BigFastaFile::make_index(string filename) {
  std::ifstream in(filename.c_str());
  if (!in)
    throw FastaFileException("cannot open input file " + filename);
  char buffer[input_buffer_size + 1];
  size_t offset = 0;
  size_t sequence_start = 0;
  string name;
  bool inside_name = true;
  vector<size_t> newlines;
  size_t char_count = 0;
  while (!in.eof()) {
    in.clear();
    in.read(buffer, input_buffer_size);
    const size_t characters_read = in.gcount();
    for (size_t position = 0; position < characters_read; ++position) {
      if (!inside_name) {
	if (buffer[position] == '>') {
	  SequenceTable.push_back(SequenceInfo(name, sequence_start, 
					       offset + position, newlines));
	  inside_name = true;
	  name.clear();
	  newlines.clear();
	  char_count = 0;
	}
	else {
	  if (buffer[position] == '\n') {
	    newlines.push_back(char_count);
	  }
	  char_count++;
	}
      }
      else if (inside_name) {
	if (buffer[position] == '\n') {
	  inside_name = false;
	  // we need the '+ 1' below because it is the next position
	  // where the sequence actually starts; the current one is
	  // the preceding newline.
	  sequence_start = offset + position + 1;
	}
	else if (buffer[position] != '>')
	  name += buffer[position];
      }
    }
    offset += characters_read;
  }
  SequenceTable.push_back(SequenceInfo(name, sequence_start,
				       offset, newlines));
  in.close();
}


typedef std::vector<BigFastaFile::SequenceInfo>::const_iterator SeqInfoIter;


SeqInfoIter
BigFastaFile::find_containing_sequence(size_t position) const {
  return lower_bound(SequenceTable.begin(), SequenceTable.end(), position);
}


vector<string>
BigFastaFile::get_names(const size_t chunk_start, 
			const size_t chunk_size) const {
  SeqInfoIter start_index = find_containing_sequence(chunk_start);
  //   SeqInfoIter end_index = (chunk_start + chunk_size < filesize) ?
  //     find_containing_sequence(chunk_start + chunk_size) + 1 :
  //     SequenceTable.end();

  SeqInfoIter end_index = find_containing_sequence(chunk_start + chunk_size);
  // This check is needed because the end could be between sequences:
  // in the name part
  if (end_index != SequenceTable.end() && 
      chunk_start + chunk_size > end_index->start)
    ++end_index;
  
  vector<string> names;
  for (SeqInfoIter i = start_index; i < end_index; ++i)
    names.push_back(i->name);
  return names;
}


/* The constructor for BigFastaFile sets the filename, builds an index
 * of the locations of sequences within the file, and sets the size of
 * the file, which is needed to make sure nobody uses the object to
 * access parts of the file that don't exist. 
 */
BigFastaFile::BigFastaFile(string fn) : filename(fn) {
  make_index(filename);
  filesize = SequenceTable.back().end;
}


vector<string>
BigFastaFile::get_sequences(const size_t chunk_start, 
			    const size_t chunk_size) const {
  
  std::ifstream in(filename.c_str());
  if (!in)
    throw FastaFileException("cannot open input file " + filename);
  
  in.seekg(chunk_start);
  char *buffer = new char[chunk_size];
  in.read(buffer, chunk_size);
  const size_t characters_read = in.gcount();
  in.close();
  
  // The index of the sequence (from the file) that contains the
  // position where the chunk starts.
  SeqInfoIter start_index = find_containing_sequence(chunk_start);

  // The index of the sequence that contains the final position of the
  // chunk
  //   SeqInfoIter end_index = (chunk_start + chunk_size < filesize) ?
  //     find_containing_sequence(chunk_start + chunk_size) + 1 :
  //     SequenceTable.end();
  
  SeqInfoIter end_index = find_containing_sequence(chunk_start + chunk_size);
  // This check is needed because the end could be between sequences:
  // in the name part
  if (end_index != SequenceTable.end() &&
      chunk_start + chunk_size > end_index->start)
    ++end_index;
  
  // actually read in the chunk of the file
  vector<string> sequences;
  for (SeqInfoIter i = start_index; i < end_index; ++i) {
    char *sequence_start, *sequence_end;
    if (i->start > chunk_start) {
      sequence_start = buffer + (i->start - chunk_start);
    }
    else sequence_start = buffer;
    if (i->end > chunk_start + characters_read) {
      sequence_end = buffer + characters_read;
    }
    else 
      sequence_end = buffer + (i->end - chunk_start);
    sequences.push_back(string());
    remove_copy_if(sequence_start, sequence_end, 
 		   back_inserter(sequences.back()),
  		   bind2nd(std::equal_to<char>(), '\n'));
  }
  // done with the buffer, let's delete it
  delete buffer;
  // turn anything not in {a,A,c,C,g,G,t,T} to an N
  FastaFile::clean(sequences);
  // make every base upper case
  FastaFile::toupper(sequences);
  return sequences;
}


/* get_first_sequence_offset(chunk_start) returns the true offset of
 * the first sequence in a chunk. That is the number of characters
 * after the start of that sequence inside the file where the chunk
 * begins. This is needed to correct for locations that are specified
 * relative to the start of the first sequence in a chunk, which may
 * have started part-way into a sequence within the file, in which
 * case a previous chunk probably contained the preceding part of that
 * sequence.
 */
size_t
BigFastaFile::get_first_sequence_offset(const size_t chunk_start) const {
  SeqInfoIter seq_id = find_containing_sequence(chunk_start);
  if (chunk_start < seq_id->start)
    return 0;
  else {
    const size_t start_of_sequence = seq_id->start;
    const size_t newlines_before_chunk_start =
      seq_id->count_newlines(start_of_sequence, chunk_start);
    return chunk_start - start_of_sequence -
      newlines_before_chunk_start;
  }
}


string
BigFastaFile::tostring() const {
  ostringstream buff;
  buff << filename << "\t" << filesize << endl;
  for (size_t i = 0; i < SequenceTable.size() - 1; ++i)
    buff << SequenceTable[i].tostring() << endl;
  if (!SequenceTable.empty())
    buff << SequenceTable.back().tostring();
  return buff.str();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////////////////////// FastaFile stuff  ///////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void 
FastaFile::read() const {
  std::ifstream in(filename.c_str());
  if (!in) throw FastaFileException("cannot open input file " + filename);
  string s, name;
  bool first_line = true;
  while (!in.eof()) {
    char buffer[input_buffer_size + 1];
    in.getline(buffer, input_buffer_size);
    if (in.gcount() == static_cast<int>(input_buffer_size))
      throw FastaFileException("Line in " + name + "\nexceeds max length: " +
			       cread::toa(input_buffer_size));
    // CORRECT FOR DOS CARRIAGE RETURNS BEFORE NEWLINES
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    if (buffer[0] == '>') {
      if (first_line == false && s.length() > 0) {
	names.push_back(name);
	sequences.push_back(s);
      }
      else first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("> "));
      s = "";
    }
    else s += buffer;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
  }
  clean();   // turn anything not in {a,A,c,C,g,G,t,T} to an N
  toupper(); // make every base upper case
}


vector<string> 
FastaFile::get_ids() const {
  if (names.empty()) read();
  vector<string> ids;
  for (vector<string>::iterator i = names.begin(); i != names.end(); ++i)
    ids.push_back(i->substr(0, i->find_first_of(" \t")));
  return ids;
}


vector<string> 
FastaFile::get_descriptions() const {
  if (names.empty()) read();
  vector<string> descriptions;
  for (vector<string>::iterator i = names.begin(); i != names.end(); ++i) {
    const size_t desc_offset = i->find_first_of(" \t");
    if (desc_offset != string::npos)
      descriptions.push_back(i->substr(desc_offset));
    else
      descriptions.push_back(string());
  }
  return descriptions;
}


void
FastaFile::clean() const {
  for (vector<string>::iterator i = sequences.begin(); 
       i != sequences.end(); ++i)
    replace_if(i->begin(), i->end(), 
	       std::not1(std::ptr_fun(valid_base)), 'N');
}


void
FastaFile::clean(vector<string>& sequences) {
  for (vector<string>::iterator i = sequences.begin(); 
       i != sequences.end(); ++i)
    replace_if(i->begin(), i->end(), 
	       std::not1(std::ptr_fun(valid_base)), 'N');
}


void
FastaFile::toupper() const {
  for (vector<string>::iterator i = sequences.begin(); 
       i != sequences.end(); ++i)
    std::transform(i->begin(), i->end(), i->begin(), ::toupper);
}


void
FastaFile::toupper(vector<string>& sequences) {
  for (vector<string>::iterator i = sequences.begin(); 
       i != sequences.end(); ++i)
    std::transform(i->begin(), i->end(), i->begin(), ::toupper);
}


size_t 
FastaFile::get_filesize() const {
  if (filesize == std::numeric_limits<size_t>::max()) {
    std::ifstream f(filename.c_str());
    if (!f.good()) {return 0;}
    size_t begin_pos = f.tellg();
    f.seekg(0, std::ios_base::end);
    size_t end_pos = f.tellg();
    f.close();
    filesize = end_pos - begin_pos;
  }
  return filesize;
}

void
FastaFile::base_comp_from_file(string fn, float *base_comp, 
			       size_t buffer_size) {
  std::ifstream fin(fn.c_str());
  bool reading_sequence = true;
  char buffer[buffer_size + 1];
  std::fill(base_comp, base_comp + alphabet_size, 0);
  while (!fin.eof()) {
    fin.clear();
    fin.get(buffer, buffer_size + 1, '\0');
    const size_t limit = static_cast<size_t>(fin.gcount());
    for (size_t i = 0; i < limit; ++i)
      if (buffer[i] == '>') 
	reading_sequence = false;
      else if (!reading_sequence && buffer[i] == '\n')
	reading_sequence = true;
      else if (reading_sequence && buffer[i] != '\n' &&
	       valid_base(buffer[i]))
	base_comp[base2int(buffer[i])]++;
  }
  fin.close();
  // get the total number of valid bases
  const float total = std::accumulate(base_comp, 
				      base_comp + alphabet_size, 0.0);
  // normalize the base composition
  std::transform(base_comp, base_comp + alphabet_size, base_comp,
		 bind2nd(std::divides<float>(), total));
}


void
FastaFile::base_comp_from_files(const vector<string>& fns, 
				vector<float>& base_comp, 
				size_t buffer_size) {
  vector<size_t> count(alphabet_size, 0);
  char buffer[buffer_size + 1];
  for (size_t i = 0; i < fns.size(); ++i) {
    std::ifstream fin(fns[i].c_str());
    if (!fin) 
      throw FastaFileException("cannot open input file " + fns[i]);
    bool reading_sequence = false;
    while (!fin.eof()) {
      fin.clear();
      fin.read(buffer, buffer_size);
      const size_t limit = static_cast<size_t>(fin.gcount());
      for (size_t j = 0; j < limit; ++j)
	if (buffer[j] == '>') 
	  reading_sequence = false;
	else if (!reading_sequence && buffer[j] == '\n') 
	  reading_sequence = true;
	else if (reading_sequence && buffer[j] != '\n' && 
		 valid_base(buffer[j])) {
	  count[base2int(buffer[j])]++;
	}
    }
    fin.close();
  }
  // get the total number of valid bases
  const double total = std::accumulate(count.begin(), count.end(), 0.0);
  base_comp.clear();
  transform(count.begin(), count.end(), back_inserter(base_comp),
	    std::bind2nd(std::divides<double>(), total));
}


vector<string>
FastaFile::read_seqs_dir(const string& dirname, const string& suffix) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw FastaFileException("could not open directory: " + dirname);
  const string dotted_suffix('.' + suffix);
  const size_t sufflen = dotted_suffix.length();
  vector<string> filenames;
  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    const string d_name_str(ent->d_name);
    if (d_name_str.length() > sufflen &&
	!d_name_str.substr(d_name_str.length() - sufflen,
			   sufflen).compare(dotted_suffix))
      filenames.push_back(dirname + "/" + d_name_str);
    errno = 0;
  }
  if (errno)
    throw FastaFileException("readdir failed with error: " + 
			     cread::toa(errno) +
			     "for directory: " + dirname);
  return filenames;
}
