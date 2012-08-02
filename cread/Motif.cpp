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

#include "Motif.hpp"

using std::string;
using std::vector;
using std::map;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::endl;

const char *
Motif::type_id = "Motif";

const size_t
Motif::type_id_size = 5;

const char *
Motif::matrix_start = "P0";

Motif::Motif(const Motif& m) : 
  Pattern(m), matrix(m.matrix), sites(m.sites) {
  type = type_id;
}

Motif::Motif(const Matrix &m) : matrix(m) {
  accession = "(none)";
  type = type_id;
}

Motif& 
Motif::operator=(const Motif& m) {
  if (this != &m) {
    Pattern::operator=(m);
    matrix = m.matrix;
    type = type_id;
    sites = m.sites;
  }
  return *this;
}

void
Motif::format_representation(ostream& os) const {
  os << matrix << endl 
     << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Motif::format_sites(ostream& os) const {
  string sep("\n");
  sep += PatternID::BINDING_SITE_START + string("  ");
  if (sites.size() > 0) {
    os << PatternID::BINDING_SITE_START << "  ";
    std::copy(sites.begin(), sites.end() - 1,
	      std::ostream_iterator<MotifSite>(os, sep.c_str()));
    os << sites.back() << endl;
    os << PatternID::BLANK_PATTERN_LINE << endl;
  }
}

// TODO: in the event of an exception in the constructor below, make
// sure the state is reset so that there are no leaks, and strange
// state.

// TODO: need to verify that the information is valid -- read_accession
// isn't even actually being used.
Motif::Motif(vector<string>& lines) : Pattern(lines) {
  bool read_matrix = false, read_sites = false;
  for (size_t i = 0; i < lines.size(); ++i) {
    string line = lines[i];
    if (line_type(line, matrix_start)) {
      vector<string> matrix_lines;
      matrix_lines.push_back(line);
      for (size_t j = i + 1; j < lines.size() && isdigit(lines[j][0]); ++j)
	matrix_lines.push_back(lines[j]);
      matrix = Matrix(matrix_lines);
      read_matrix = true;
      i += matrix_lines.size();
    }
    else if (!read_sites &&
	     line_type(line, PatternID::BINDING_SITE_START)) {
      while (line_type(lines[i], PatternID::BINDING_SITE_START))
	sites.push_back(MotifSite(remove_line_id(lines[i++])));
      read_sites = true;
    }

    // TODO: need to verify format better in here
  }
  if (type.empty())
    type = type_id;
  if (!read_matrix)
    throw MotifFormatException();
}

vector<Motif>
Motif::ReadMotifVector(string file_name) {
  vector<vector<string> > motif_lines;
  ReadPatternLines(file_name, motif_lines);
  vector<Motif> motifs;
  for (size_t i = 0; i < motif_lines.size(); ++i) {
    bool version_record = false;
    for (size_t j = 0; j < motif_lines[i].size(); ++j) {
      if (line_type(motif_lines[i][j], version_line)) 
	version_record = true;
    }
    if (!version_record)
      motifs.push_back(Motif(motif_lines[i]));
  }
  return motifs;
}
