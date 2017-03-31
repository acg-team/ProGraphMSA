#include <fstream>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <string>

#include "TempFile.h"

std::vector<std::string> TempFile::getDirs(const char* directory) {
	std::vector<std::string> dirs;

	char *tmpdir = getenv("TMPDIR");
	if(tmpdir) dirs.push_back(std::string(tmpdir));

	tmpdir = getenv("TEMP");
	if(tmpdir) dirs.push_back(std::string(tmpdir));

	if(directory) dirs.push_back(std::string(directory));

	dirs.push_back(std::string(P_tmpdir));

	dirs.push_back(std::string("/tmp"));

	return dirs;
}

TempFile::TempFile(const char* prefix, const char* directory) {
	std::vector<std::string> dirs = getDirs(directory);

	if(!prefix) prefix = "tmpfile-";

	for(std::vector<std::string>::const_iterator dir=dirs.begin(); dir != dirs.end(); ++dir) {
		std::string fname = *dir + delim + prefix + XXX;

		char *file = strdup(fname.c_str());
		int fd = mkstemp(file);
		if(fd >= 0) {
			this->filename = std::string(file);
			free(file);
			this->open(this->filename.c_str(), std::fstream::trunc | std::fstream::in | std::fstream::out);
			::close(fd);

			if(*this) {
				return;
			} else {
				unlink(this->filename.c_str());
			}
		} else {
			free(file);
		}
	}
	this->filename = "";
	this->setstate(std::fstream::failbit);
}

TempFile::~TempFile() {
	this->close();
	if(this->filename != "") {
		unlink(this->filename.c_str());
	}
}

std::string TempFile::getFilename() const {
	return this->filename;
}

const std::string TempFile::delim = std::string("/");
const std::string TempFile::XXX = std::string("XXXXXX");
