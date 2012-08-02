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

#include "Pattern.hpp"

using std::string;
using std::vector;
using std::map;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::endl;

// TODO: in the format_whatever() functions, get rid of the hardcoded
// strings and replace them with the constants defined here:
const size_t Pattern::PatternID::id_width = 2;
const char *Pattern::PatternID::ACCESSION_START = "AC";
const char *Pattern::PatternID::TYPE_START = "TY";
const char *Pattern::PatternID::BASIS_START = "BA";
const char *Pattern::PatternID::BINDING_FACTOR_START = "BF";
const char *Pattern::PatternID::BINDING_SITE_START = "BS";
// NEED TO FIX THIS:
const char *Pattern::PatternID::COMMENT_START = "CO";
const char *Pattern::PatternID::DESCRIPTION_START = "DE";
const char *Pattern::PatternID::ATTRIBUTE_START = "AT";
const char *Pattern::PatternID::IDENTIFIER_START = "ID";
const char *Pattern::PatternID::FACTOR_NAME_START = "NA";
const char *Pattern::PatternID::REFERENCE_AUTHOR_START = "RA";
const char *Pattern::PatternID::REFERENCE_TITLE_START = "RT";
const char *Pattern::PatternID::REFERENCE_START = "RN";
const char *Pattern::PatternID::VERSION_LINE = "VV";
const char *Pattern::PatternID::BLANK_PATTERN_LINE = "XX";
const char *Pattern::PatternID::PATTERN_TERMINATOR = "//";
const char *Pattern::PatternID::ORGANIZATION_LABEL = "OR";


const char *
Pattern::version_line = "VV";

// TODO: verify that the format is correct
Pattern::Pattern(vector<string>& lines) {
  bool read_accession = false, read_type = false;
  for (size_t i = 0; i < lines.size(); ++i) {
    string line = lines[i];
    if (line_type(line, PatternID::ACCESSION_START)) {
      set_accession(remove_line_id(line));
      read_accession = true;
    }
    else if (read_accession) {
      if (line_type(line, PatternID::TYPE_START)) {
	type = remove_line_id(line);
	read_type = true;
      }
      if (line_type(line, PatternID::COMMENT_START)) 
	set_comment(remove_line_id(line));
      if (line_type(line, PatternID::ATTRIBUTE_START)) {
	pair<string, string> att = parse_attribute_line(line);
	set_attribute(att.first, att.second);
      }
      if (line_type(line, PatternID::IDENTIFIER_START)) 
	identifier = remove_line_id(line);
      if (line_type(line, PatternID::FACTOR_NAME_START)) 
	factor_names = remove_line_id(line);
      if (line_type(line, PatternID::DESCRIPTION_START)) 
	factor_description = remove_line_id(line);
      if (line_type(line, PatternID::BASIS_START)) 
	basis = remove_line_id(line);

      if (line_type(line, PatternID::ORGANIZATION_LABEL))
	organization = remove_line_id(line);

    }
    else {
      string bad_line = line;
      if (!line_type(line, PatternID::BLANK_PATTERN_LINE) &&
	  !line_type(line, PatternID::VERSION_LINE)) 
	throw PatternFileException("Accession must come first in a pattern, "
				   "not: " + bad_line, i);
    }
  }
  if (!read_type)
    throw PatternFileException("Missing pattern type");
}

Pattern::Pattern(const Pattern& p) {
  accession = p.accession;
  type = p.type;
  comment = p.comment;

  attributes = p.attributes;
  // binding_sites = p.binding_sites;

  basis = p.basis;
  binding_factors = p.binding_factors;
  factor_description = p.factor_description;
  identifier = p.identifier;
  factor_names = p.factor_names;

  organization = p.organization;

}

Pattern& 
Pattern::operator=(const Pattern& p) {
  if (this != &p) {
    accession = p.accession;
    type = p.type;
    // binding_sites = p.binding_sites;
    comment = p.comment;

    attributes = p.attributes;

    basis = p.basis;
    binding_factors = p.binding_factors;
    factor_description = p.factor_description;
    identifier = p.identifier;
    factor_names = p.factor_names;
    organization = p.organization;

    //   reference_authors = p.reference_authors;
    //   reference_data = p.reference_data;
    //   reference_titles = p.reference_titles;
    //   references = p.references;
  }
  return *this;
}

string
Pattern::tostring() const {
  ostringstream os;
  format_accession(os);
  format_type(os);
  format_identifier(os);
  format_factor_names(os);

  format_organization(os);

  format_factor_description(os);
  
  format_representation(os);
  
  format_comment(os);
  format_attributes(os);
  
  format_basis(os);
  format_binding_factors(os);
  
  format_sites(os);
  os << PatternID::PATTERN_TERMINATOR;
  return os.str();
}


pair<string, string>
Pattern::parse_attribute_line(string s) {
  s = remove_line_id(s);
  size_t equal_sign_offset = s.find_first_of("=");
  string key = s.substr(0, equal_sign_offset);
  size_t value_offset = s.find_first_of(" ", equal_sign_offset);
  string value = s.substr(equal_sign_offset + 1, 
			  value_offset - equal_sign_offset);
  return pair<string, string>(key, value);
}

string
Pattern::get_attribute(string s) const {
  std::map<string, string>::const_iterator a = attributes.find(s);
  if (a != attributes.end()) 
    return a->second;
  else return "";
}


void 
Pattern::format_attributes(ostream& os) const {
  for (map<string, string>::const_iterator i = attributes.begin();
       i != attributes.end(); ++i)
    os << "AT  " << i->first << "=" << i->second << endl;
  if (!attributes.empty())
    os << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Pattern::format_accession(ostream& os) const {
  if (!accession.empty())
    os << "AC  " << accession << endl 
       << PatternID::BLANK_PATTERN_LINE << endl;
  else throw PatternFormatterException(string("bad accession format") +
				       accession);
}

void
Pattern::format_type(ostream& os) const {
  if (!type.empty())
    os << PatternID::TYPE_START << "  " << type << endl
       << PatternID::BLANK_PATTERN_LINE << endl;
  else throw PatternFormatterException(string("bad pattern type: ") +
				       type);
}

void
Pattern::format_comment(ostream& os) const {
  if (!comment.empty())
    os << "CO  " << comment << endl 
       << PatternID::BLANK_PATTERN_LINE << endl;
}


void
Pattern::format_basis(ostream& os) const {
  if (!basis.empty())
    os << PatternID::BASIS_START << "  " << basis << endl 
       << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Pattern::format_binding_factors(ostream& os) const {
  if (!binding_factors.empty())
    os << PatternID::BINDING_FACTOR_START << "  " << binding_factors << endl 
       << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Pattern::format_factor_description(ostream& os) const {
  if (!factor_description.empty())
    os << PatternID::DESCRIPTION_START << "  " 
       << factor_description << endl 
       << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Pattern::format_identifier(ostream& os) const {
  if (!identifier.empty())
    os << PatternID::IDENTIFIER_START << "  " << identifier 
       << endl << PatternID::BLANK_PATTERN_LINE << endl;
}

void
Pattern::format_factor_names(ostream& os) const {
  if (!factor_names.empty())
    os << PatternID::FACTOR_NAME_START << "  " << factor_names << endl 
       << PatternID::BLANK_PATTERN_LINE << endl;
}


void
Pattern::format_organization(ostream& os) const {
  if (!organization.empty())
    os << PatternID::ORGANIZATION_LABEL << "  " << organization << endl
       << PatternID::BLANK_PATTERN_LINE << endl;
}


// template <class T> void
// Motif<T>::format_reference_authors(ostream& os)
// {
//   if (!reference_authors.empty()) 
//     os << reference_authors;
// }

// template <class T> void
// Motif<T>::format_reference_titles(ostream& os)
// {
//   os << reference_titles;
// }

// template <class T> void
// Motif<T>::format_references(ostream& os)
// {
//   os << references;
// }


void
Pattern::separate_pattern_lines(vector<string> &in, 
				vector<vector<string> > &out) {
  out.clear();
  out.push_back(vector<string>());
  size_t i;
  for (i = 0; i < in.size() - 1; ++i) {
    if (line_type(in[i], PatternID::PATTERN_TERMINATOR))
      out.push_back(vector<string>());
    else out.back().push_back(in[i]);
  }
  if (!line_type(in[i], PatternID::PATTERN_TERMINATOR))
    throw PatternFileException("Bad pattern terminator: " + string(in[i]), i);
}

bool
Pattern::line_type(string& line, const char *type) {
  return !line.compare(0, PatternID::id_width, type);
}

string
Pattern::remove_line_id(string& s) {
  return s.substr(s.find_first_not_of("\t ", PatternID::id_width));
}

// TODO: make this more efficient by using the fact that only blank
// lines may come between the accession and the type, and the
// accession must come first in the file.
string
Pattern::extract_type(vector<string>& lines) {
  for (size_t i = 0; i < lines.size(); ++i)
    if (line_type(lines[i], PatternID::TYPE_START))
      return remove_line_id(lines[i]);
  
  return "";
}

bool
PatternOrder::operator()(const Pattern *m1, const Pattern *m2) const {
  if (m1->attributes.find(key) == m1->attributes.end() ||
      m2->attributes.find(key) == m1->attributes.end())
    throw BadKeyException(key, m1->get_accession());
  string val1 = m1->attributes.find(key)->second;
  string val2 = m2->attributes.find(key)->second;
  if (reverse && numeric) return atof(val1.c_str()) > atof(val2.c_str());
  else if (reverse) return val1 > val2;
  else if (numeric) return atof(val1.c_str()) < atof(val2.c_str());
  else return val1 < val2;
}

bool
PatternOrder::operator()(const Pattern &m1, const Pattern &m2) const {
  if (m1.attributes.find(key) == m1.attributes.end() ||
      m2.attributes.find(key) == m1.attributes.end())
    throw BadKeyException(key, m1.get_accession());
  string val1 = m1.attributes.find(key)->second;
  string val2 = m2.attributes.find(key)->second;
  if (reverse && numeric) return atof(val1.c_str()) > atof(val2.c_str());
  else if (reverse) return val1 > val2;
  else if (numeric) return atof(val1.c_str()) < atof(val2.c_str());
  else return val1 < val2;
}

bool
PatternCutoff::operator()(const Pattern &p) const {
  if (p.attributes.find(key) == p.attributes.end())
    throw BadKeyException(key, p.get_accession());
  string val = p.attributes.find(key)->second;
  if (reverse && numeric) 
    return atof(val.c_str()) > atof(value.c_str());
  else if (reverse) 
    return val > value;
  else if (numeric) 
    return atof(val.c_str()) < atof(value.c_str());
  else return val < value;
}

bool
PatternCutoff::operator()(const Pattern *p) const {
  if (p->attributes.find(key) == p->attributes.end())
    throw BadKeyException(key, p->get_accession());
  string val = p->attributes.find(key)->second;
  if (reverse && numeric) 
    return atof(val.c_str()) > atof(value.c_str());
  else if (reverse) 
    return val > value;
  else if (numeric) 
    return atof(val.c_str()) < atof(value.c_str());
  else return val < value;
}

void
Pattern::ReadPatternLines(string file_name,
			  vector<vector<string> > &pattern_lines) {
  std::ifstream in(file_name.c_str());
  if (!in)
    throw PatternFileException("cannot open pattern file " +
			       string(file_name));
  vector<string> lines;
  while (!in.eof()) {
    char buffer[pattern_line_size+1];
    in.getline(buffer, pattern_line_size);
    if (strlen(buffer) > 0) 
      lines.push_back(buffer);
  }
  if (lines.empty())
    throw PatternFileException("could not read lines in pattern file " +
			       string(file_name));
  Pattern::separate_pattern_lines(lines, pattern_lines);
}
