#include "RepeatDetectionTReks.h"
#include "Alphabet.h"
#include <assert.h>
#include <fstream>

#ifdef WITH_CODON
template <>
std::map< std::string, std::vector<repeat_t> > detect_repeats<Codon>(const std::map<std::string,sequence_t<Codon> > &seqs)
{
	std::string filename;
	TempFile tmp("tmpseqrep-");

	if(!tmp) {
		throw treks_exception("could not create temporary file");
	}

	std::map<std::string,std::string> seqs2;
	for (std::map<std::string, sequence_t<Codon> >::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
		seqs2[(*it).first] = stringFromSequence<AA>(translateCodons((*it).second));
	}
	write_fasta(seqs2,tmp);
	tmp.flush();

	if(!tmp) {
		throw treks_exception("could not write temporary file");
	}

	std::map< std::string, std::vector<repeat_t> > result;
	result = do_detect_repeats(tmp.getFilename(),seqs2);

	return result;
}
#endif

static std::string strip(const std::string &s) {
	size_t start = s.find_first_not_of(" \t\f\v\n\r");
	if(start == std::string::npos) return "";

	size_t end = s.find_last_not_of(" \t\f\v\n\r");
	return s.substr(start, end-start+1);
}

/**
 * Parses repeats from a T-REKS (.trd) file.
 *
 * Checks the repeats against the full-length sequence for the protein.
 *
 * Example TREKS:
 *
 *     >sp|Q10567|AP1B1_HUMAN AP-1 complex subunit beta-1 OS=Homo sapiens OX=9606 GN=AP1B1 PE=1 SV=2
 *     Length: 4 residues - nb: 5  from  625 to 645 - Psim:0.7 region Length:21
 *     G-D-LL
 *     G-D-LL
 *     NLD-L-
 *     G-PPVS
 *     G-P-PL
 * @param ss File stream
 * @param seqs Map from IDs to the sequence string (parsed earlier from FASTA or other source)
 * @return map from the ID to a list of repeats
 */
static std::map< std::string, std::vector<repeat_t> > parse_treks_output(std::istream *ss, const std::map<std::string,std::string> &seqs)
{
	std::map< std::string, std::vector<repeat_t> > map;

	index_t n_sequences = 0;
	index_t n_repeats = 0;

	std::string name;
	std::string line;

	while(std::getline(*ss,line)) {
		if(line[0] == '>') {
			// identifier
			name = strip(line.substr(1));
			++n_sequences;
		} else if (line.compare(0,7,"Length:") == 0) {
			// header line
			size_t from = line.find("from");
			if(from == line.npos) {
				throw treks_exception("format error (from)");
			}

			size_t to = line.find("to",from);
			if(to == line.npos) {
				throw treks_exception("format error (to)");
			}

			repeat_t repeat;
			++n_repeats;

			std::stringstream ssstart;
			ssstart << line.substr(from+4,to-from-4);
			ssstart >> repeat.start;
			if(!ssstart || repeat.start <= 0) {
				throw treks_exception("format error (number)");
			}
			repeat.start -= 1;

			std::map<std::string,std::string>::const_iterator orig_entry = seqs.find(name);
			if(orig_entry == seqs.end()) {
				throw treks_exception("unknown sequence name: " + name);
			}
			std::string::const_iterator orig = orig_entry->second.begin() + repeat.start;

			repeat.len = (index_t)-1;
			int line_no = 0;
			while(std::getline(*ss,line)) {
				line = strip(line);
				++line_no;
				if(line.compare(0,22,"**********************") == 0) break;

				for(index_t i=0; i<line.size(); ++i) {
					char c = line[i];
					if(c == '-' || c == ' ' || c == '\n' || c == '\t' || c == '\r') {
						line[i] = '_';
					}
				}

				if(repeat.len != (index_t)-1 && line.size() != repeat.len) {
					throw treks_exception("repeat unit lengths differ");
				}
				repeat.len = line.size();

				for(index_t i=0; i<repeat.len; ++i) {
					if(line[i] != '_') {
						repeat.tr_hom.push_back(i);
						if(*orig != line[i]) {
							int orig_index =  orig - orig_entry->second.begin() + 1;
							int line_index = i + 1;
							std::stringstream ss;
							ss << "character mismatch (repeat " << n_repeats <<
									", seq \"" << name << "\", orig pos " <<
									orig_index << " char '" << line[i] <<
									"', trmsa line " << line_no << " pos " <<
									line_index << " char '" << *orig << "')";
							throw treks_exception(ss.str());
						}
						++orig;
					}
				}
			}

			map[name].push_back(repeat);
		}
		// ignore comments, blank lines, etc
	}

	std::cerr << "found " << n_repeats << " repeats in " << n_sequences << " sequences" << std::endl;

	return(map);
}

std::map< std::string, std::vector<repeat_t> > do_read_repeats(const std::string &filename, const std::map<std::string,std::string> &seqs)
{
	std::ifstream f(filename.c_str());
	return parse_treks_output(&f,seqs);
}

std::map< std::string, std::vector<repeat_t> > do_detect_repeats(const std::string &filename, const std::map<std::string,std::string> &seqs)
{
	char buffer[BUFFERSIZE];
	std::iostream *io;

	if(cmdlineopts.trdout_file == "") {
		io = new std::stringstream(std::stringstream::in | std::stringstream::out | std::stringstream::app);
	} else {
		io = new std::fstream(cmdlineopts.trdout_file.c_str(), std::fstream::in | std::fstream::out | std::fstream::trunc);
	}
	std::string cmd = std::string(java_path) + " -jar " + treks_jar_path + " -infile=\"" + filename + "\"";
	if(cmdlineopts.customtr_cmd != "") {
		cmd = cmdlineopts.customtr_cmd + " \"" + filename + "\"";
	}

	std::FILE *fd = popen(cmd.c_str(),"r");
	if(fd == NULL) {
		throw treks_exception("could not launch Java/T-Reks");
	}

	index_t ret;
	index_t total = 0;
	while((ret = std::fread(buffer, 1, BUFFERSIZE-1, fd)) > 0) {
		total += ret;
		buffer[ret] = '\0';
		*io << buffer;
	}

	pclose(fd);
	io->flush();
	io->seekg(0);

	std::map< std::string, std::vector<repeat_t> > repeats;
	repeats = parse_treks_output(io,seqs);

	delete io;

	return repeats;
}
