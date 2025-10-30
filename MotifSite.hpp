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

#ifndef MOTIFSITE_HPP_
#define MOTIFSITE_HPP_

#include <cstdint>
#include <ostream>
#include <string>

struct MotifSite {
  std::string site;
  std::string seq_name;
  int start{};
  std::uint32_t length{};
  std::string gaps;
  char orientation{};
  float score{};

  auto
  is_valid_orientation(const char c) -> bool;

  MotifSite(const std::string &site, const std::string &seq_name,
            const int start, const std::uint32_t length,
            const std::string &gaps, const char orientation,
            const float score) :
    site{site}, seq_name{seq_name}, start{start}, length{length}, gaps{gaps},
    orientation{orientation}, score{score} {}

  [[nodiscard]] auto
  get_site() const -> const std::string & {
    return site;
  }

  [[nodiscard]] auto
  get_seq_name() const -> const std::string & {
    return seq_name;
  }

  [[nodiscard]] auto
  get_gaps() const -> const std::string & {
    return gaps;
  }

  [[nodiscard]] auto
  get_orientation() const -> char {
    return orientation;
  }

  [[nodiscard]] auto
  negstrand() const -> bool {
    return orientation == 'n';
  }

  [[nodiscard]] auto
  posstrand() const -> bool {
    return orientation == 'p';
  }

  [[nodiscard]] auto
  get_score() const -> float {
    return score;
  }

  [[nodiscard]] auto
  get_start() const -> int {
    return start;
  }

  [[nodiscard]] auto
  get_length() const -> std::uint32_t {
    return length;
  }

  auto
  set_score(float s) {
    score = s;
  }
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

#endif  // MOTIFSITE_HPP_
