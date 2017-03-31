/** \file debug.h
 *
 * Routines for error reporting.
 */
/*
 * Copyright (c) 2007-2011 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _DEBUG_H
#define _DEBUG_H

#include <exception>
#include <string>

class swps3_exception: public std::exception
{
public:
  swps3_exception(std::string message) {
	  this->str = message;
  }

  ~swps3_exception() throw() {}

  virtual const char* what() const throw()
  {
    return this->str.c_str();
  }

private:
  std::string str;
};

extern void error(const char* str, ...);
extern void warning(const char* str, ...);

#ifdef DEBUG
#	include <cstdio>
#	define debug(...) { std::fprintf(std::stderr,__VA_ARGS__); std::fflush(std::stderr); }
#else
#	define debug(...)
#endif

#endif /* _DEBUG_H */
