/** \file fasta.c
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

#include "Fasta.h"
#include "main.h"
#include "debug.h"
#include <fstream>
#include <iostream>
#include <string>

FastaLib::FastaLib ( const char * filename ) {
	this->fp = new std::fstream(filename, std::fstream::in);

	if(!*this->fp) {
		throw fasta_exception(std::string("error opening file"));
	}

	if (this->fp->peek() != '>') {
		throw fasta_exception(std::string("format error"));
	}
}

static std::string strip(const std::string &s) {
	size_t start = s.find_first_not_of(" \t\f\v\n\r");
	if(start == std::string::npos) return "";

	size_t end = s.find_last_not_of(" \t\f\v\n\r");
	return s.substr(start, end-start+1);
}

seq_t FastaLib::nextSequence( ) {
	seq_t out;

	if (this->fp->eof()) {
		out.valid = false;
		return out;
	}

	out.valid = true;

	this->fp->ignore(); // ignore '>'
	std::getline(*this->fp, out.name);
	out.name = strip(out.name);

	std::string line = "";
	while (*this->fp && this->fp->peek() != '>') {
		line = "";
		std::getline(*this->fp, line);
		line = strip(line);
		out.seq += line;
	}

	return out;
}

std::map<std::string,std::string> FastaLib::readAll( ) {
	std::map<std::string,std::string> result;
	seq_t seq;

	while((seq = this->nextSequence())) {
		if(result.find(seq.name) != result.end()) {
			error("duplicate sequence name \"%s\"",seq.name.c_str());
		}
		result[seq.name] = seq.seq;
	}

	return result;
}

std::map<std::string,std::string> FastaLib::readAll(std::vector<std::string> &order) {
	std::map<std::string,std::string> result;
	seq_t seq;

	while((seq = this->nextSequence())) {
		if(result.find(seq.name) != result.end()) {
			error("duplicate sequence name \"%s\"",seq.name.c_str());
		}
		result[seq.name] = seq.seq;
		order.push_back(seq.name);
	}

	return result;
}


FastaLib::~FastaLib() {
	delete this->fp;
}

void write_fasta(const std::map<std::string,std::string> &alignment, std::ostream &out) {
	for(std::map<std::string,std::string>::const_iterator it=alignment.begin(); it != alignment.end(); ++it) {
		out << ">" << it->first << std::endl << it->second << std::endl;
	}
}

void write_fasta(const std::map<std::string,std::string> &alignment, const std::vector<std::string> &order, std::ostream &out) {
	for(std::vector<std::string>::const_iterator it=order.begin(); it != order.end(); ++it) {
		out << ">" << *it << std::endl << alignment.at(*it) << std::endl;
	}
}
