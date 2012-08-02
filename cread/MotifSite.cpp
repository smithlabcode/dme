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

#include "MotifSite.hpp"

using std::string;
using std::ostringstream;
using std::vector;


bool
MotifSite::is_valid_orientation(char c) {
  // static const char* MotifSite::valid_orientations = "pn";
  // static const size_t n_valid_orientations = 2;
  return c == 'p' || c == 'n';
}


MotifSite::MotifSite(std::string s, std::string sn, std::string st, 
		     std::string l, std::string g, std::string o, std::string sc) :
  site(s), seq_name(sn), gaps(g) {
  orientation = (is_valid_orientation(o[0])) ? o[0] : 'p';
  start = atoi(st.c_str());
  score = atof(sc.c_str());
  int temp_length = atoi(l.c_str());
  if (temp_length > 0)
    length = static_cast<size_t>(temp_length);
  else throw InvalidMotifSiteException(
	      static_cast<string>("non-positive site length"));
}

MotifSite::MotifSite(string s) {
  vector<string> parts = cread::split(s, ";", true);
  if (parts.size() < 6) {
    string message("invalid site format: ");
    message += s;
    throw InvalidMotifSiteException(message);
  }
  site = cread::strip(parts[0]);
  seq_name = cread::strip(parts[1]);
  start = atoi(parts[2].c_str());
  // TODO: make sure the length is sensible
  length = static_cast<size_t>(atoi(parts[3].c_str()));
  gaps = cread::strip(parts[4]);  // this is always empty
  orientation = cread::strip(parts[5])[0];
  if (parts.size() > 6)
    score = atof(parts[6].c_str());
  else score = -std::numeric_limits<float>::max();
}

string
MotifSite::tostring() const {
  ostringstream os;
  os << site << "; " << seq_name << "; " << start << "; " 
     << length << "; " << gaps << "; " << orientation << "; " << score;
  return os.str();
}

bool
MotifSite::operator<(const MotifSite& ms) const {
  return ((seq_name < ms.seq_name) || 
	  (seq_name == ms.seq_name && start < ms.start) || 
	  (seq_name == ms.seq_name && start == ms.start && 
	   length < ms.length) || 
	  (seq_name == ms.seq_name && start == ms.start && 
	   length == ms.length && orientation < ms.orientation) || 
	  (seq_name == ms.seq_name && start == ms.start && 
	   length == ms.length && orientation == ms.orientation &&
	   score < ms.score));
}

bool 
MotifSite::operator==(const MotifSite& ms) const {
  return site == ms.site && seq_name == ms.seq_name && 
    start == ms.start && length == ms.length && 
    orientation == ms.orientation && score == ms.score;
}
