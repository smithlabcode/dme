/* Copyright (C) 2025 Andrew D. Smith
 *
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#ifndef MOTIF_HPP
#define MOTIF_HPP

#include "Matrix.hpp"
#include "MotifSite.hpp"

#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

struct Motif {
  void
  add_site(const MotifSite &bs) {
    sites.push_back(bs);
  }

  std::string accession;
  std::string identifier;

  float score{};
  std::uint32_t fgcount{};
  std::uint32_t bgcount{};
  std::uint32_t correctedbgcount{};
  float info{};

  Matrix matrix;
  std::vector<MotifSite> sites;

  [[nodiscard]] std::string
  tostring() const {
    static constexpr auto blank_line = "XX";
    std::ostringstream os;

    os << "AC  " << accession << '\n' << blank_line << '\n';
    os << "ID  " << identifier << '\n' << blank_line << '\n';
    os << matrix << '\n' << blank_line << '\n';
    os << "AT  FGCOUNT=" << fgcount << '\n';
    os << "AT  BGCOUNT=" << bgcount << '\n';
    os << "AT  CORRECTEDBGCOUNT=" << correctedbgcount << '\n';
    os << "AT  INFO=" << info << '\n';
    os << "AT  SCORE=" << score << '\n';
    os << blank_line << '\n';

    for (const auto &s : sites)
      os << "BS  " << s << '\n';
    if (!sites.empty())
      os << blank_line << '\n';
    os << "//";
    return os.str();
  }
};

inline std::ostream &
operator<<(std::ostream &o, const Motif &m) {
  return o << m.tostring();
}

#endif
