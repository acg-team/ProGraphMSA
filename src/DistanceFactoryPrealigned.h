#ifndef DISTANCEFACTORYPREALIGNED_H_
#define DISTANCEFACTORYPREALIGNED_H_

#include "main.h"
#include "DistanceFactoryML.h"
#include "ModelFactory.h"

template <class ALPHABET>
class DistanceFactoryPrealigned : public DistanceFactoryML<ALPHABET> {
private:
	typedef typename DistanceFactoryML<ALPHABET>::CountMatrix CountMatrix;

public:
	DistanceFactoryPrealigned();
	DistanceFactoryPrealigned(const ModelFactory<ALPHABET> *model_factory) : DistanceFactoryML<ALPHABET>(model_factory) {};
	virtual ~DistanceFactoryPrealigned() {};
	virtual DistanceMatrix computePwDistances(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::vector<std::string> &order);
};


/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include <iostream>
#include <fstream>

template <class ALPHABET>
DistanceMatrix DistanceFactoryPrealigned<ALPHABET>::computePwDistances(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::vector<std::string> &order)
{
	DistanceMatrix distances(order.size());
	Eigen::Matrix<int, Eigen::Dynamic, 1> seq_len(order.size());

	for(index_t i=0; i < order.size(); ++i) {
		const std::string &ii = order[i];
		typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator iii = sequences.find(ii);
		assert(iii != sequences.end());
		const sequence_t<ALPHABET> &seq1 = iii->second;

        seq_len(i) = seq1.size();

		for(index_t j=i+1; j < order.size(); ++j) {
			const std::string &jj = order[j];
			typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator jjj = sequences.find(jj);
			assert(jjj != sequences.end());
			const sequence_t<ALPHABET> &seq2 = jjj->second;

			CountMatrix counts = CountMatrix::Zero();
			index_t gaps = 0;

			assert(seq1.length() == seq2.length());

			bool gap_opened1 = false;
			bool gap_opened2 = false;

			for(index_t k=0; k<seq1.length(); ++k) {
				if(!seq1[k].isGap() && !seq2[k].isGap()) {
					int c1 = seq1[k].value();
					int c2 = seq2[k].value();
					if(c1 >= 0 && c1 < 20 && c2 >= 0 && c2 < 20)
						++counts(c1,c2);
					gap_opened1 = false;
					gap_opened2 = false;
				} else if (seq1[k].isGap() && seq2[k].isGap()) {
					// skip
				} else if(!seq1[k].isGap() && !gap_opened1) {
					++gaps;
					gap_opened1 = true;
					gap_opened2 = false;
				} else if(!seq2[k].isGap() && !gap_opened2) {
					++gaps;
					gap_opened1 = false;
					gap_opened2 = true;
				}
			}

			distvar_t distvar = this->computeDistance(counts,gaps,(seq1.length()+seq2.length())/2.0);

			distances.distances(i, j) = distances.distances(j, i) = distvar.dist;
			distances.variances(i, j) = distances.variances(j, i) = distvar.var;
		}
	}

	return distances;
}

#endif /* DISTANCEFACTORYPREALIGNED_H_ */
