#ifndef DISTANCEFACTORYANGLE_H_
#define DISTANCEFACTORYANGLE_H_

#include "main.h"
#include "DistanceFactory.h"
#include "Alphabet.h"

template<int X, int K>
class Pow {
public:
    enum{
    	value = X*Pow<X, K-1>::value
    };
};

template<int X>
class Pow<X, 0> {
public:
    enum{
    	value = 1
    };
};

template <class ALPHABET,int K>
class DistanceFactoryAngle : public DistanceFactory<ALPHABET> {
private:
	typedef Eigen::Matrix<int, Eigen::Dynamic, Pow<ALPHABET::DIM,K>::value> CountMatrix;
	typedef Eigen::Matrix<double, Pow<ALPHABET::DIM,K>::value, Eigen::Dynamic> CountMatrixSVD;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> SingularValues;

public:
	DistanceFactoryAngle();
	virtual ~DistanceFactoryAngle() {};
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

#include "debug.h"
#include <Eigen/Dense>

template <class ALPHABET,int K>
DistanceFactoryAngle<ALPHABET,K>::DistanceFactoryAngle()
{
}

template <class ALPHABET,int K>
DistanceMatrix DistanceFactoryAngle<ALPHABET,K>::computePwDistances(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::vector<std::string> &order)
{
	DistanceMatrix distances(order.size());
	CountMatrix counts = CountMatrix::Zero((int)order.size(),Pow<ALPHABET::DIM,K>::value);
	Eigen::Matrix<int, Eigen::Dynamic, 1> seq_len(order.size());

	for(index_t i=0; i < order.size(); ++i) {
		const std::string &ii = order[i];
		typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator iii = sequences.find(ii);
		assert(iii != sequences.end());
		const sequence_t<ALPHABET> &seq = iii->second;

		seq_len(i) = seq.size();

		int chars[K];
		for(index_t k=0; k<K; ++k) {
			chars[k] = -1;
		}

		for(index_t j=0; j<seq.size(); ++j) {
			for(index_t k=1; k<K; ++k) {
				chars[k-1] = chars[k];
			}
			chars[K-1] = seq[j].value();
			if(chars[K-1] < 0 || chars[K-1] >= ALPHABET::DIM) chars[K-1] = -1;

			int index = 0;
			for(index_t k=0; k<K; ++k) {
				if(chars[k] == -1) {
					index = -1;
					break;
				}

				index *= ALPHABET::DIM;
				index += chars[k];
			}

			if(index != -1) {
				counts(i,index) += 1;
			}
		}
	}

#if 0
	Eigen::JacobiSVD<CountMatrixSVD> svd(counts.template cast<double>().transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);


	SingularValues sv = svd.singularValues();
	double sv_sum = sv.sum();
	double sv_sum_part = 0;
	for(index_t i=0; i < (index_t)sv.rows(); ++i) {
		if(sv_sum_part > .95*sv_sum) sv(i) = 0;
		sv_sum_part += sv(i);
	}

	CountMatrixSVD counts2 = svd.matrixU().leftCols(sv.rows()) * sv.asDiagonal() * svd.matrixV().transpose();
#else
	CountMatrixSVD counts2 = counts.template cast<double>().transpose();
#endif

	distances.distances = counts2.colwise().norm().array().inverse().matrix().asDiagonal() * counts2.transpose() * counts2 * counts2.colwise().norm().array().inverse().matrix().asDiagonal();
	distances.distances = -1.0 * ((distances.distances.array().square() + 0.4) / 1.4).log();
	if(!(cmdlineopts.mldist_flag) && !(cmdlineopts.mldist_gap_flag)) {
		distances.distances = distances.distances.array().exp();
		distances.distances = - 0.5 * (5.0*distances.distances.array() - (45.0*distances.distances.array().square() - 20.0*distances.distances.array()).sqrt()) * distances.distances.array().inverse();
	}

	distances.variances.setConstant(1);
	distances.variances = distances.variances * seq_len.cast<distance_t>().asDiagonal();
	distances.variances = ((distances.variances + distances.variances.transpose())/2).eval().array().inverse();
	distances.variances.array() *= distances.distances.array();

	//adjust minimum variances
	distances.variances.array() = distances.variances.array().max(1e-5 * Eigen::Array<distance_t, Eigen::Dynamic, Eigen::Dynamic>::Ones(distances.variances.rows(),distances.variances.cols()));

	return distances;
}

#endif /* DISTANCEFACTORYANGLE_H_ */
