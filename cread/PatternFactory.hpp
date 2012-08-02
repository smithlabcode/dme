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

#include "cread.hpp"
#include "Pattern.hpp"
#include "Motif.hpp"
#include "Module.hpp"
#include "Word.hpp"

class PatternFactory {
public:
  static Pattern *CreatePattern(std::vector<std::string>&, std::string = "");
  static std::vector<Pattern*> ReadPatternVector(std::string, std::string = "");
private:
  PatternFactory() {}
};

class PatternFactoryException : public CREADException {
public:
  PatternFactoryException(const std::string m) : CREADException(m) {}
  PatternFactoryException() throw() {}
};
