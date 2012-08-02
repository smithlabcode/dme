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

#include "PatternFactory.hpp"

using std::vector;
using std::string;

Pattern *
PatternFactory::CreatePattern(vector<string>& pattern_lines, string type_id) {
  Pattern *r;
  if (type_id.empty()) 
    type_id = Pattern::extract_type(pattern_lines);
  if (Motif::is_type(type_id)) 
    r = new Motif(pattern_lines);
  else if (Word::is_type(type_id))
    r = new Word(pattern_lines);
  else if (Module::is_type(type_id))
    r = new Module(pattern_lines);
  else 
    throw PatternFactoryException(string("invalid Pattern type: ") + type_id);
  return r;
}

vector<Pattern*>
PatternFactory::ReadPatternVector(string file_name, string type_id) {
  vector<vector<string> > pattern_lines;
  Pattern::ReadPatternLines(file_name, pattern_lines);
  vector<Pattern*> patterns;
  for (size_t i = 0; i < pattern_lines.size(); ++i)
    patterns.push_back(CreatePattern(pattern_lines[i], type_id));
  return patterns;
}
