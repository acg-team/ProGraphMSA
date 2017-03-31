#ifndef TEMPFILE_H_
#define TEMPFILE_H_

#include <fstream>
#include <vector>
#include <string>

class TempFile : public std::fstream {
private:
	std::string filename;
	std::vector<std::string> getDirs(const char* directory=NULL);
	static const std::string delim;
	static const std::string XXX;

public:
	TempFile(const char* prefix=NULL, const char* directory=NULL);
	virtual ~TempFile();
	std::string getFilename() const;
};

#endif /* TEMPFILE_H_ */
