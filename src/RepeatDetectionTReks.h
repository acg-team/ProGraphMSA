#ifndef REPEATDETECTIONTREKS_H_
#define REPEATDETECTIONTREKS_H_

#include <exception>
#include <string>
#include <map>
#include <vector>
#include "debug.h"
#include "Alphabet.h"
#include "Repeat.h"

class treks_exception: public swps3_exception {
public:
	treks_exception(std::string str) : swps3_exception(str) {};
};

template <class ALPHABET>
typename std::map< std::string, std::vector<repeat_t> > read_repeats(const std::string &filename, const std::map<std::string,sequence_t<ALPHABET> > &seqs);
template <class ALPHABET>
typename std::map< std::string, std::vector<repeat_t> > detect_repeats(const std::map<std::string,sequence_t<ALPHABET> > &seqs, const std::string &trd_output_file);
#ifdef WITH_CODON
template <>
std::map< std::string, std::vector<repeat_t> > detect_repeats<Codon>(const std::map<std::string,sequence_t<Codon> > &seqs, const std::string &trd_output_file);
#endif

/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include <iostream>
#include <sstream>
#include <cstdio>
#include "main.h"
#include "TempFile.h"
#include "Fasta.h"
#include "external_programs.h"

#define BUFFERSIZE 1024

std::map< std::string, std::vector<repeat_t> > do_detect_repeats(const std::string &filename, const std::map<std::string,std::string> &seqs);
std::map< std::string, std::vector<repeat_t> > do_read_repeats(const std::string &filename, const std::map<std::string,std::string> &seqs);

template <class ALPHABET>
typename std::map< std::string, std::vector<repeat_t> > detect_repeats(const std::map<std::string,sequence_t<ALPHABET> > &seqs)
{
	std::string filename;
	TempFile tmp("tmpseqrep-");

	if(!tmp) {
		throw treks_exception("could not create temporary file");
	}

	std::map<std::string,std::string> seqs2;
	for (typename std::map<std::string, sequence_t<ALPHABET> >::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
		seqs2[(*it).first] = stringFromSequence<ALPHABET>((*it).second);
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

template <class ALPHABET>
typename std::map< std::string, std::vector<repeat_t> > read_repeats(const std::string &filename, const std::map<std::string,sequence_t<ALPHABET> > &seqs)
{
	std::map<std::string,std::string> seqs2;
	for (typename std::map<std::string, sequence_t<ALPHABET> >::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
		seqs2[(*it).first] = stringFromSequence<ALPHABET>((*it).second);
	}

	std::map< std::string, std::vector<repeat_t> > result;
	result = do_read_repeats(filename, seqs2);

	return result;
}

#include "ProgressiveAlignment.h"
#include "TreeNJ.h"
#include "Model.h"
#include "CSProfile.h"

template <class ALPHABET>
std::map< std::string, std::vector<repeat_t> > align_repeats(const typename std::map<std::string, sequence_t<ALPHABET> > &seqs, const std::map< std::string, std::vector<repeat_t> > &reps, const CSProfile *csprofile, const ModelFactory<ALPHABET> *model_factory) {
	std::map< std::string, std::vector<repeat_t> > dummy_reps;
	std::map< std::string, std::vector<repeat_t> > new_reps;

	for(typename std::map< std::string, std::vector<repeat_t> >::const_iterator it = reps.begin(); it != reps.end(); ++it) {
		const std::string &seq_name = it->first;
		const sequence_t<ALPHABET> &seq = seqs.at(seq_name);
		const std::vector<repeat_t> &cur_repeat_list = it->second;
		std::vector<repeat_t> new_repeat_list;

		for(typename std::vector<repeat_t>::const_iterator it2 = cur_repeat_list.begin(); it2 != cur_repeat_list.end(); ++it2) {
			std::map<std::string,sequence_t<ALPHABET> > seqs2;
			std::vector<std::string> units;
			const repeat_t &cur_repeat = *it2;

			index_t unit=0;
			index_t start=0;
			while(start < cur_repeat.tr_hom.size()) {
				std::stringstream ss;
				ss << unit;
				std::string sunit = ss.str();
				units.push_back(sunit);

				index_t end = start+1;
				while(end < cur_repeat.tr_hom.size() && cur_repeat.tr_hom[end] > cur_repeat.tr_hom[end-1]) ++end;

				seqs2[sunit] = seq.substr(cur_repeat.start+start,end-start);

				++unit;
				start = end;
			}

			PhyTree *tree = TreeNJ<ALPHABET>(seqs2,false,model_factory);
			std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > alignment_cache;
			ProgressiveAlignmentResult<ALPHABET> result = progressive_alignment<ALPHABET>(seqs2,*tree,dummy_reps,csprofile,*model_factory,alignment_cache);
			delete tree;

			repeat_t new_repeat;
			new_repeat.start = cur_repeat.start;
			new_repeat.len = result.aligned_sequences[units.front()].size();
			new_repeat.tr_hom.reserve(cur_repeat.tr_hom.size());
			for(std::vector<std::string>::const_iterator it = units.begin(); it != units.end(); ++it) {
				const sequence_t<ALPHABET> &s = result.aligned_sequences[*it];
				for(index_t i=0; i < new_repeat.len; ++i) {
					if(!s[i].isGap()) {
						new_repeat.tr_hom.push_back(i);
					}
				}
			}
			assert(cur_repeat.tr_hom.size() == new_repeat.tr_hom.size());

			new_repeat_list.push_back(new_repeat);
		}
		new_reps[seq_name] = new_repeat_list;
	}

	return new_reps;
}

#endif /* REPEATDETECTIONTREKS_H_ */
