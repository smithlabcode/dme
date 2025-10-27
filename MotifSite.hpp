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

#include <cstdint>
#include <ostream>
#include <string>

struct MotifSite {
  std::string site;
  std::string seq_name;
  int start;
  std::uint32_t length;
  std::string gaps;
  char orientation;
  float score;

  bool
  is_valid_orientation(char c);

  MotifSite(std::string s, std::string sn, int st, std::uint32_t l,
            std::string g, char o, float sc) :
    site(s), seq_name(sn), start(st), length(l), gaps(g), orientation(o),
    score(sc) {}

  std::string
  get_site() const {
    return site;
  }
  std::string
  get_seq_name() const {
    return seq_name;
  }
  std::string
  get_gaps() const {
    return gaps;
  }
  char
  get_orientation() const {
    return orientation;
  }
  bool
  negstrand() const {
    return orientation == 'n';
  }
  bool
  posstrand() const {
    return orientation == 'p';
  }
  float
  get_score() const {
    return score;
  }

  int
  get_start() const {
    return start;
  }

  std::uint32_t
  get_length() const {
    return length;
  }

  void
  set_score(float s) {
    score = s;
  }
  // bool
  // operator<(const MotifSite &ms) const;
  // bool
  // operator==(const MotifSite &ms) const;
};

inline std::ostream &
operator<<(std::ostream &o, const MotifSite &s) {
  // clang-format off
  return o << s.site << "; "
          << s.seq_name << "; "
          << s.start << "; "
          << s.length << "; "
          << s.gaps << "; "
          << s.orientation << "; "
          << s.score;
  // clang-format on
}

#endif
