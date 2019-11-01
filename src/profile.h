#ifndef PROFILE_H_
#define PROFILE_H_

#include <iostream>
#include <map>
#include <string>

#include "main.h"
#include "Alphabet.h"
#include "Model.h"

template <class ALPHABET>
void write_profile(const std::map<std::string,typename Model<ALPHABET>::Profile> &profiles, std::ostream &out)
{
	for(typename std::map<std::string,typename Model<ALPHABET>::Profile>::const_iterator i = profiles.begin(); i != profiles.end(); ++i) {
		const std::string name = i->first;
		const typename Model<ALPHABET>::Profile &profile = i->second;

		out << '>' << name << std::endl;
		for(int j=0; j < ALPHABET::DIM; ++j) {
			out << ALPHABET(j).asString();

			for(index_t k=0; k < (index_t)profile.cols(); ++k) {
				out << '\t' << profile(j,k);
			}

			out << std::endl;
		}
	}
}

#endif /* PROFILE_H_ */
