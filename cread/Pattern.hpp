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

#ifndef PATTERN_HPP
#define PATTERN_HPP

/*!
   \file Pattern.hpp
   \brief Definition for the Pattern base class.
 */

#include "cread.hpp"

/*!
   \class Pattern
   \brief This base class provides a standard interface, exception handling
	  and functionality that is common to patterns.
 */
class Pattern {
public:

  //! Empty constructor
  Pattern() {}
  /*!
     \brief This constructor reads a vector of strings to populate the pattern
            with standard data types included in PatternID.
   */
  Pattern(std::vector<std::string>&);
  //! Copy constructor
  Pattern(const Pattern &);
  //! Assignment operator
  Pattern& operator=(const Pattern &);
  //! Empty destructor
  virtual ~Pattern() {}
  
  // ACCESSORS
  //! Returns pattern accession
  std::string get_accession() const {return accession;}
  //! Returns pattern type
  std::string get_type() const {return type;}
  //! Returns pattern comment
  std::string get_comment() const {return comment;}

  /*!
    \brief Returns the value of the attribute if it exists, empty string otherwise.
    \param[in] s Attribute.
    \return Attribute value or empty string if the attribute is not set.
   */
  std::string get_attribute(std::string s) const;
  
  //! Returns the pattern basis
  std::string get_basis() const {return basis;}
  //! Returns the binding factors
  std::string get_binding_factors() const {return binding_factors;}
  //! Returns the factor description
  std::string get_factor_description() const {return factor_description;}
  //! Returns the pattern identifier
  std::string get_identifier() const {return identifier;}
  //! Returns the factor names
  std::string get_factor_names() const {return factor_names;}
  //! Returns the pattern organization
  std::string get_organization() const {return organization;}
  
  //! Returns the entire pattern attributes map
  std::map<std::string, std::string> get_attributes() const {
    return attributes;
  }
  //! String representation for the pattern
  virtual std::string tostring() const;
  
  // MUTATORS
  //! Sets the pattern accession
  void set_accession(std::string s) {accession = s;}
  //! Sets the pattern comment
  void set_comment(std::string s) {comment = s;}
  //! Sets the pattern identifier
  void set_identifier(std::string s) {identifier = s;}
  
  /*!
    \brief Sets the attribute key to value.
    \param[in] key Attribute name.
    \param[in] value Attribute value.
   */
  template <class U, class V> void set_attribute(U key, V value) {
    attributes[cread::toa(key)] = cread::toa(value);
  }

  //! Extract the pattern type from a set of strings
  static std::string extract_type(std::vector<std::string>&);
  /*!
    \brief Split a vector of strings containing multiple patterns into
	   a vector of string vectors such that each string vector contains
	   exactly one pattern.
    \param[in] in
	A vector of strings containing multiple patterns.
    \param[out] out
	A vector of string vectors such that each string vector contains
	data on one pattern.
  */
  static void separate_pattern_lines(std::vector<std::string> & in, 
				     std::vector<std::vector<std::string> >
								      & out);
  
  //! Print the pattern
  friend std::ostream& operator<<(std::ostream& s, const Pattern &p) {
    return s << p.tostring();
  }

protected:

  /*!
   \struct PatternID
   \brief Standard pattern data type labels.
   */
  struct PatternID {
    static const size_t id_width;
    static const char *ACCESSION_START;
    static const char *TYPE_START;
    static const char *BASIS_START;
    static const char *BINDING_FACTOR_START;
    static const char *BINDING_SITE_START;
    static const char *COMMENT_START;
    static const char *DESCRIPTION_START;
    static const char *ATTRIBUTE_START;
    static const char *IDENTIFIER_START;
    static const char *FACTOR_NAME_START;
    static const char *REFERENCE_AUTHOR_START;
    static const char *REFERENCE_TITLE_START;
    static const char *REFERENCE_START;
    static const char *VERSION_LINE;
    static const char *BLANK_PATTERN_LINE;
    static const char *PATTERN_TERMINATOR;
    static const char *ORGANIZATION_LABEL;
  };


  /*
    \brief Parse the input string to identify an attribute and its value.
    \param[in] s String containing an attribute and its value.
    \return Pair of strings, where first contains the attribute key
	    and the second contains its value.
    \warning The input string must contain an attribute and its value;
	     value is empty string otherwise.
   */
  static std::pair<std::string, std::string> parse_attribute_line(std::string s);
  //! Returns true if the line starts with the label type
  static bool line_type(std::string& line, const char *type);
  //! Returns a string where the type level has been removed
  static std::string remove_line_id(std::string& s);

  std::string accession;
  std::string type;
  std::string comment;
  std::map<std::string, std::string> attributes;

  std::string basis;
  std::string binding_factors;
  std::string factor_description;
  std::string identifier;
  std::string factor_names;
  std::string organization;

  /*! \brief This pure virtual function provides an interface for each pattern
      to print its representation.
   */
  virtual void format_representation(std::ostream &os) const = 0;
  /*! \brief This pure virtual function provides an interface for each pattern
      to print its sites.
   */
  virtual void format_sites(std::ostream& os) const = 0;
  
  //! If pattern accession exists, print it into ostream os
  void format_accession(std::ostream& os) const;
  //! If pattern type exists, print it into ostream os
  void format_type(std::ostream& os) const;
  //! If pattern comment exists, print it into ostream os
  void format_comment(std::ostream& os) const;
  //! For each attribute, print its key and value into ostream os
  void format_attributes(std::ostream& os) const;

  //! If pattern basis exists, print it into ostream os
  void format_basis(std::ostream& os) const;
  //! If pattern binding factors information exists, print it into ostream os
  void format_binding_factors(std::ostream& os) const;
  //! If pattern factor description exists, print it into ostream os
  void format_factor_description(std::ostream& os) const;
  //! If pattern identifier exists, print it into ostream os
  void format_identifier(std::ostream& os) const;
  //! If pattern factor names exists, print it into ostream os
  void format_factor_names(std::ostream& os) const;
  //! If pattern organization exists, print it into ostream os
  void format_organization(std::ostream& os) const;

  //! \brief This class provides comparison functionality
  friend class PatternOrder;
  //! This class provides cutoff comparison functionality
  friend class PatternCutoff;
  //! This container class needs access to type information for error checking
  friend class PatternFactory;

  /*! 
    \brief Retrieve lines from a file.
    \param[in] file_name File name.
    \param[out] pattern_lines A Vector of strings set to the line in file_name.
    \warning Line maximum length is pattern_line_size.
   */
  static void ReadPatternLines(std::string file_name, 
			       std::vector<std::vector<std::string> >&
							pattern_lines);

  //! Version is a special non-type information
  static const char *version_line;
  //! Maximum length per line
  const static size_t pattern_line_size = 1024;
};

/*!
  \exception BadKeyException
  \brief Used to handle erros of uninitialized or missing order
  attributes.
 */
class BadKeyException : public CREADException {
  std::string key;
  std::string accession;
public:
  /// The constructor takes an attribute and a pattern name
  BadKeyException(std::string k, std::string a) throw() : 
     CREADException(), key(k), accession(a) {}
  const char * what() const throw() {
    std::ostringstream s;
    s << message <<"\n"
      << "key:\t" << key << "\n"
      << "pattern:\t" << accession;
    return s.str().c_str();
  }
  virtual ~BadKeyException() throw() {}
};

/*!
  \class PatternOrder
  \brief Provides operators to compare patters.
 */
class PatternOrder : public std::binary_function<Pattern*, Pattern*, bool> {
protected:
  //! The attribute used for comparison
  std::string key;
  //! Set for descending order, ascending by default
  bool reverse;
  //! Set for numeric (default), string comparison otherwise
  bool numeric;
public:
  //! The initialization constructor
  explicit PatternOrder(std::string k, bool r = false, bool n = true) : 
    key(k), reverse(r), numeric(n) {}
  //! Comparison function for pattern pointers
  virtual bool operator()(const Pattern *m1, const Pattern *m2) const;
  //! Comparison function
  virtual bool operator()(const Pattern &m1, const Pattern &m2) const;
  //! Empty destructor
  virtual ~PatternOrder() {}
};

/*!
  \class PatternCutoff
  \brief Provides operators to compare values to a predetermined cutoff
 */
class PatternCutoff : public std::binary_function<Pattern*, std::string, bool> {
protected:
  //! The attribute used for comparison
  std::string key;
  //! The cutoff value
  std::string value;
  //! Set for descending order, ascending by default
  bool reverse;
  //! Set for numeric (default), string comparison otherwise
  bool numeric;
public:
  //! The initialization constructor
  explicit PatternCutoff(std::string k, std::string v, 
			 bool r = false, bool n = true) : 
    key(k), value(v), reverse(r), numeric(n) {}
  //! Comparison function for a pattern pointer
  virtual bool operator()(const Pattern *p) const;
  //! Comparison function
  virtual bool operator()(const Pattern &p) const;
  //! Empty destructor
  virtual ~PatternCutoff() {}
};

/*!
  \exception PatternFormatterException
  \brief Used for reporting errors when writing patterns.
 */
class PatternFormatterException : public CREADException {
public:
  /// Constructor that initializes the message
  PatternFormatterException(std::string m = std::string()) throw() : CREADException (m) {}
};

/*!
  \exception PatternFormatException
  \brief Used for reporting errors when reading patterns.
 */
struct PatternFormatException {
  //! Constructor that initializes the message
  PatternFormatException(std::string m) : message(m) {}
  //! Empty constructor
  PatternFormatException() {}
  //! The message
  std::string message;
};

/*!
  \exception PatternFileException
  \brief Used for reporting errors when reading patterns from a file.
  \param message The message.
  \param line_number The line number where the error was detected.
 */
class PatternFileException : public CREADException {
  //! The line number where the error was detected
  int line_number;
public:
  //! The initialization constructor takes optional message and line number
  PatternFileException(std::string m = "", int l = -1) : 
    CREADException(m), line_number(l) {}
  /*! 
    \brief Return a string composed of the message and the line number
    indication position in the file in which the error was detected.
   */
  const char * what() const throw() {
    std::ostringstream s;
    s << message;
    if (line_number > 0) s << " (line " << line_number << ")";
    return s.str().c_str();
  }
};

/*!
  \exception PatternSiteException
  \brief Used for reporting errors in pattern sites.
  \param message The message.
 */
class PatternSiteException : public CREADException {
public:
  //! The initialization constructor takes optional message
  PatternSiteException(std::string m = std::string()) : CREADException(m) {}
};

#endif
