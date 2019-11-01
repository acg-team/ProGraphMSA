#ifndef MODEL_H_
#define MODEL_H_

#include <Eigen/Core>
#include "main.h"
#include "Alphabet.h"

template <class ALPHABET>
struct Model {
	typedef Eigen::Matrix<score_t,ALPHABET::DIM,ALPHABET::DIM> Subst;
	typedef Eigen::Matrix<score_t,ALPHABET::DIM,1> Freqs;
	typedef Eigen::Matrix<score_t,ALPHABET::DIM,Eigen::Dynamic> Profile;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Subst M;
	Subst P;
	Subst Q;
	Freqs pi;
	double delta;
	double epsilon;

	distance_t distance;
	distance_t divergence;
};

#endif /* MODEL_H_ */
