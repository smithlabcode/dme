/* MIT License
 *
 * Copyright (c) 2025 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
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
