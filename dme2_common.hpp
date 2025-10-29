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

#ifndef DME2_COMMON_HPP_
#define DME2_COMMON_HPP_

#include <array>
#include <cstddef>
#include <cstdint>

inline constexpr std::uint32_t alphabet_size = 4;  // DNA!

typedef std::array<float, alphabet_size> Column;

static constexpr std::uint8_t encode_base[128] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 0
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 17
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 33
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 49
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  //@,A-O
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // P-Z
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  //`,a-o
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // p-z
};

static constexpr char encoding_offset = 64;

// clang-format off
static constexpr char complement[] = {
  'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 't', 'N', 'g', 'N', 'N', 'N', 'c',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'a', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
};
// clang-format on

[[nodiscard]] inline auto
get_complement(const char c) -> char {
  return c > encoding_offset
           ? complement[static_cast<std::uint8_t>(c - encoding_offset)]
           : 4;
}

// clang-format off
static constexpr bool valid_encoding[] = {
  false, true,  false, true,  false, false, false, true,
  false, false, false, false, false, false, false, false,
  false, false, false, false, true,  false, false, false,
  false, false, false, false, false, false, false, false,
  false, true,  false, true,  false, false, false, true,
  false, false, false, false, false, false, false, false,
  false, false, false, false, true,  false, false, false,
  false, false, false, false, false, false, false, false,
};
// clang-format on

[[nodiscard]] inline auto
isvalid(const char c) -> char {
  return c > encoding_offset
           ? valid_encoding[static_cast<std::uint8_t>(c - encoding_offset)]
           : false;
}

#endif  // DME2_COMMON_HPP_
