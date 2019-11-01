#ifndef MODELFACTORYCUSTOM_H_
#define MODELFACTORYCUSTOM_H_

#include <iostream>
#include "ModelFactory.h"
#include "main.h"
#include "debug.h"

class custom_model_exception: public swps3_exception {
public:
	custom_model_exception(std::string str) : swps3_exception(str) {};
};

template <class ALPHABET>
class ModelFactoryCustom : public ModelFactory<ALPHABET> {
public:
	ModelFactoryCustom(std::istream &input) throw (custom_model_exception);
	virtual ~ModelFactoryCustom() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include <sstream>
#include <cmath>
#include <Eigen/Dense>

template <class ALPHABET>
ModelFactoryCustom<ALPHABET>::ModelFactoryCustom(std::istream &input) throw (custom_model_exception) {
   for(index_t i=1; i<ALPHABET::DIM; ++i) {
	   for(index_t j=0; j<i; ++j) {
		   score_t value = 0;
		   input >> value;

		   if(!input) throw custom_model_exception("error reading exchangeability matrix from file");
		   if(!(value > 0 && value < INFINITY)) throw custom_model_exception("negative/infinity/zero value in exchangeability matrix");

		   this->Q(j,i) = this->Q(i,j) = value;
	   }
   }

   for(index_t i=0; i<ALPHABET::DIM; ++i) {
	   score_t value = 0;
	   input >> value;

	   if(!input) throw custom_model_exception("error reading amino acid frequencies");
	   if(!(value > 0 && value < INFINITY)) throw custom_model_exception("negative/infinity/zero value in amino acid frequencies");

	   this->freqs(i) = value;
   }

   this->freqs /= this->freqs.sum();

   // normalize rate
   this->Q.diagonal().setZero();
   this->Q.diagonal() = -this->Q.rowwise().sum().eval();
   this->Q /= -(this->freqs.transpose() * this->Q.diagonal())(0,0);

   Eigen::EigenSolver<typename Model<ALPHABET>::Subst> solver(this->Q);
   this->sigma = solver.eigenvalues().real();
   this->V = solver.eigenvectors().real();
   this->Vi = this->V.inverse();
}

#endif /* MODELFACTORYCUSTOM_H_ */
