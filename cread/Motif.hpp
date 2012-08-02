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

#ifndef MOTIF_HPP
#define MOTIF_HPP
/**
  \file Motif.hpp

  \brief This header file contains class definitions for Motif class,
  a subclass of Pattern.
*/

#include "cread.hpp"
#include "Matrix.hpp"
#include "MotifSite.hpp"
#include "Pattern.hpp"

/**
  \brief Class to represent DNA sequence (e.g. TF binding site)
  motifs.

  \class Motif
 */
class Motif : public Pattern {
public:
  /// This constructor parses a vector of strings to populate the motif pattern
  Motif(std::vector<std::string>&);
  // The copy constructor
  Motif(const Motif&);
  /*!
   \brief This constructor creates a minimal motif pattern from a matrix;
	  the accession is set to (none).
  */
  Motif(const Matrix&);
  /// Assignment operator
  Motif& operator=(const Motif&);

  // ACCESSORS
  /// Returns the matrix
  Matrix get_matrix() const {return matrix;}
  /// Returns the MotifSites
  std::vector<MotifSite> get_sites() const {return sites;}

  /// Add a site that can not be modified
  void add_site(const MotifSite& bs) {sites.push_back(MotifSite(bs));}
  /// Add a site
  void add_site(MotifSite& bs) {sites.push_back(MotifSite(bs));}
  /// Parse a site from a string and add it
  void add_site(std::string bs) {sites.push_back(MotifSite(bs));}
  /// Delete all sites
  void clear_sites() {sites.clear();}
  
  /// Returns a constant matrix by reference
  const Matrix& const_get_matrix() const {return matrix;}
  /// Returns the width of the matrix
  size_t get_width() const {return matrix.get_width();}
  /// Returns true if the type is a motif type
  static bool is_type(std::string s) {
    return !s.compare(0, type_id_size, type_id);
  }

  /// Read a vector of motifs from a file
  static std::vector<Motif> ReadMotifVector(std::string file_name);

private:
  /// Stores the type id
  static const char *type_id;
  /// Stores the size of the type-id label
  static const size_t type_id_size;
  /// Stores the label of the matrix start
  static const char *matrix_start;
  
  void format_representation(std::ostream& os) const;
  void format_sites(std::ostream& os) const;
  
  Matrix matrix;
  std::vector<MotifSite> sites;
  
};

/*!
  \exception ModuleFormatterException
  \brief Reports exceptions when writing a motif.
*/
class MotifFormatterException : public PatternFormatterException {};
/*!
  \exception ModuleFormatterException
  \brief Reports exceptions when parsing a Motif.
*/
class MotifFormatException : public PatternFormatException {};

#endif
