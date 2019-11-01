/** \file fasta.h
 *
 * Routines for reading FASTA databases.
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

#ifndef FASTA_H
#define FASTA_H

#include "main.h"
#include "debug.h"
#include <fstream>
#include <exception>
#include <string>
#include <map>
#include <vector>

class fasta_exception: public swps3_exception {
public:
	fasta_exception(std::string str) : swps3_exception(str) {};
};

typedef struct seq {
	std::string seq;
	std::string name;
	bool valid;

	operator bool() const {
		return valid;
	}
} seq_t;

class FastaLib {
public:
	FastaLib( const char * filename );
	~FastaLib();
	seq_t nextSequence();
	std::map<std::string,std::string> readAll();
	std::map<std::string,std::string> readAll(std::vector<std::string> &order);

private:
	std::fstream *fp;
};

void write_fasta(const std::map<std::string,std::string > &alignment, std::ostream &out);
void write_fasta(const std::map<std::string,std::string > &alignment, const std::vector<std::string> &order, std::ostream &out);

#endif /* FASTA_H */
