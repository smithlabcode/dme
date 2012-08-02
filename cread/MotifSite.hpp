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

#ifndef MOTIFSITE_HPP
#define MOTIFSITE_HPP

#include "cread.hpp"
#include "Pattern.hpp"

class MotifSite {
private:
  
  // static const char* valid_orientations = "pn";
  // static const size_t n_valid_orientations = 2;
  
  std::string site;
  std::string seq_name;
  int start;
  size_t length;
  std::string gaps;
  char orientation;
  float score;

  bool is_valid_orientation(char c);
  
public:
  MotifSite(std::string s, std::string sn, int st, 
	    size_t l, std::string g, char o, float sc) : 
    site(s), seq_name(sn), start(st), length(l), gaps(g), 
    orientation(o), score(sc) {}
  
  MotifSite(std::string s, std::string sn, std::string st, 
	    std::string l, std::string g, std::string o, std::string sc);
  MotifSite(std::string);
  std::string tostring() const;
  friend std::ostream& operator<<(std::ostream &s, const MotifSite &bs) {
    return s << bs.tostring();
  }
  
  std::string get_site() const {return site;}
  std::string get_seq_name() const {return seq_name;}
  std::string get_gaps() const {return gaps;}
  char get_orientation() const {return orientation;}
  bool negstrand() const {return orientation == 'n';}
  bool posstrand() const {return orientation == 'p';}
  float get_score() const{return score;}

  int get_start() const {return start;}
  size_t get_length() const {return length;}

  void set_score(float s) {score = s;}
  bool operator<(const MotifSite& ms) const;
  bool operator==(const MotifSite& ms) const;
  
};

/*!
   \class InvalidMotifSiteException
   \brief Exception class for handling invalid motif site exceptions
*/
class InvalidMotifSiteException : public PatternSiteException {
public:
  //! Constructor that initializes the message
  InvalidMotifSiteException(const std::string m) : PatternSiteException(m) {}
};

#endif
