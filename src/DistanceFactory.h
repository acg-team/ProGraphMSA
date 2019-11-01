#ifndef DISTANCEFACTORY_H_
#define DISTANCEFACTORY_H_

#include <Eigen/Core>
#include <map>
#include <vector>
#include <string>

#include "Alphabet.h"
#include "ModelFactory.h"

struct DistanceMatrix {
	DistanceMatrix(int dim) {
		distances = Eigen::Matrix<distance_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(dim,dim);
		variances = Eigen::Matrix<distance_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(dim,dim);
	}
	Eigen::Matrix<distance_t, Eigen::Dynamic, Eigen::Dynamic> distances;
	Eigen::Matrix<distance_t, Eigen::Dynamic, Eigen::Dynamic> variances;

	DistanceMatrix reduce(int i) {
		int dim = this->distances.rows();
		DistanceMatrix reduced_dist(dim - 1);

		if(i == 0) {
			reduced_dist.distances = this->distances.bottomRightCorner(dim-1,dim-1);
			reduced_dist.variances = this->variances.bottomRightCorner(dim-1,dim-1);
		} else if (i == dim - 1) {
			reduced_dist.distances = this->distances.topLeftCorner(dim-1,dim-1);
			reduced_dist.variances = this->variances.topLeftCorner(dim-1,dim-1);
		} else {
			reduced_dist.distances << this->distances.topLeftCorner(i,i), this->distances.topRightCorner(i,dim-i-1),
					this->distances.bottomLeftCorner(dim-i-1,i), this->distances.bottomRightCorner(dim-i-1,dim-i-1);

			reduced_dist.variances << this->variances.topLeftCorner(i,i), this->variances.topRightCorner(i,dim-i-1),
					this->variances.bottomLeftCorner(dim-i-1,i), this->variances.bottomRightCorner(dim-i-1,dim-i-1);
		}

		return reduced_dist;
	}

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <class ALPHABET>
class DistanceFactory {
public:
	virtual DistanceMatrix computePwDistances(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::vector<std::string> &order) = 0;
	virtual ~DistanceFactory() {};

	static DistanceFactory<ALPHABET> *getDefault(const ModelFactory<ALPHABET> *model_factory=NULL, bool prealigned=false);
};

#endif /* DISTANCEFACTORY_H_ */
