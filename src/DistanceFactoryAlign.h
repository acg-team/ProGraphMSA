#ifndef DISTANCEFACTORYALIGN_H_
#define DISTANCEFACTORYALIGN_H_

#include "main.h"
#include "DistanceFactoryML.h"
#include "Alphabet.h"

template <class ALPHABET>
class DistanceFactoryAlign : public DistanceFactoryML<ALPHABET> {
private:
	typedef int dp_int_t;
	typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> DynProgMatrix;
	typedef typename DistanceFactoryML<ALPHABET>::CountMatrix CountMatrix;

public:
	DistanceFactoryAlign() { initMatrix(); };
	DistanceFactoryAlign(const ModelFactory<ALPHABET> *model_factory) : DistanceFactoryML<ALPHABET>(model_factory) { initMatrix(); };
	virtual ~DistanceFactoryAlign() {};
	virtual DistanceMatrix computePwDistances(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::vector<std::string> &order);

protected:
	distvar_t alignPair(const sequence_t<ALPHABET> &seq1, const sequence_t<ALPHABET> &seq2);
	void initMatrix();
	Eigen::Matrix<int,ALPHABET::DIM+1,ALPHABET::DIM+1> scoring_matrix;
	int gap_open;
	int gap_extend;
};

template <class ALPHABET>
DistanceMatrix DistanceFactoryAlign<ALPHABET>::computePwDistances(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::vector<std::string> &order)
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

			distvar_t distvar = this->alignPair(seq1,seq2);
			distances.distances(i,j) = distances.distances(j,i) = distvar.dist;
			distances.variances(i,j) = distances.variances(j,i) = distvar.var;
		}
	}

	return distances;
}

template <class ALPHABET>
distvar_t DistanceFactoryAlign<ALPHABET>::alignPair(const sequence_t<ALPHABET> &seq1, const sequence_t<ALPHABET> &seq2) {
	const int minfty = -10000;
	DynProgMatrix W((int)seq2.size()+1,(int)seq1.size()+1);
	DynProgMatrix X((int)seq2.size()+1,(int)seq1.size()+1);
	DynProgMatrix Y((int)seq2.size()+1,(int)seq1.size()+1);

	Eigen::Matrix<int,Eigen::Dynamic,1> s1((int)seq1.size()+1);
	Eigen::Matrix<int,Eigen::Dynamic,1> s2((int)seq2.size()+1);

	W(0,0) = 0;

	for(index_t x=1; x<seq1.size()+1; ++x) {
		s1(x) = seq1[x-1].value();
		if(s1(x) < 0) s1(x) = 20;
		X(0,x) = W(0,x) = this->gap_open + (x-1)*this->gap_extend;
		Y(0,x) = minfty;
	}

	for(index_t y=1; y<seq2.size()+1; ++y) {
		s2(y) = seq2[y-1].value();
		if(s2(y) < 0) s2(y) = 20;
		Y(y,0) = W(y,0) = this->gap_open + (y-1)*this->gap_extend;
		X(y,0) = minfty;
	}

	for(index_t y=1; y<seq2.size()+1; ++y) {
		for(index_t x=1; x<seq1.size()+1; ++x) {
			W(y,x) = W(y-1,x-1)+this->scoring_matrix(s2(y),s1(x));
			X(y,x) = std::max(X(y,x-1)+this->gap_extend,W(y,x-1)+this->gap_open);
			Y(y,x) = std::max(Y(y-1,x)+this->gap_extend,W(y-1,x)+this->gap_open);
			W(y,x) = std::max(std::max(X(y,x),Y(y,x)),W(y,x));
		}
	}

	CountMatrix counts = CountMatrix::Zero();
	index_t gaps = 0;

	bool gap_opened1 = false;
	bool gap_opened2 = false;

	// backtrack and count identical characters in matches
	for (index_t y = seq2.size(), x = seq1.size(); y != 0 && x != 0;) {
		if (W(y, x) == W(y - 1, x - 1) + this->scoring_matrix(s2(y), s1(x))) {
			if(s1(x) < ALPHABET::DIM && s2(y) < ALPHABET::DIM)
				++counts(s1(x),s2(y));
			gap_opened1 = false;
			gap_opened2 = false;
			--x;
			--y;
		} else if (W(y, x) == X(y, x)) {
			if(!gap_opened1)
				++gaps;
			gap_opened1 = true;
			gap_opened2 = false;
			--x;
		} else if (W(y, x) == Y(y, x)) {
			if(!gap_opened2)
				++gaps;
			gap_opened1 = false;
			gap_opened2 = true;
			--y;
		} else {
			error("error while backtracking");
		}
	}

	distvar_t distvar = this->computeDistance(counts,gaps,(seq1.length()+seq2.length())/2.0);
	return distvar;
}

template <>
void DistanceFactoryAlign<AA>::initMatrix();
#ifdef WITH_CODON
template <>
void DistanceFactoryAlign<Codon>::initMatrix();
#endif
#ifdef WITH_DNA
template <>
void DistanceFactoryAlign<DNA>::initMatrix();
#endif

#endif /* DISTANCEFACTORYALIGN_H_ */
