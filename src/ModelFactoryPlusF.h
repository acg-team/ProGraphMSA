#ifndef MODELFACTORYPLUSF_H_
#define MODELFACTORYPLUSF_H_

#include <vector>
#include <string>

#include "ModelFactory.h"
#include "main.h"

template <class ALPHABET>
class ModelFactoryPlusF : public ModelFactory<ALPHABET> {
public:
	ModelFactoryPlusF(const ModelFactory<ALPHABET> *base_model);
	virtual ~ModelFactoryPlusF();
	virtual double getEpsilon(distance_t distance) const;
	virtual double getDelta(distance_t distance) const;
	void estimateFreqs(const std::vector<sequence_t<ALPHABET> > sequences);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	typename Model<ALPHABET>::Freqs freqsOld;
	const ModelFactory<ALPHABET> *base_model;
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
#include <sstream>
#include <cmath>
#include "debug.h"
#include <Eigen/Dense>

template <class ALPHABET>
ModelFactoryPlusF<ALPHABET>::ModelFactoryPlusF(const ModelFactory<ALPHABET> *base_model) {
	Model<ALPHABET> model = base_model->getModel(1);
	this->Q = model.Q;
	this->freqs = this->freqsOld = model.pi;

	Eigen::EigenSolver<typename Model<ALPHABET>::Subst> solver(this->Q);
	this->sigma = solver.eigenvalues().real();
	this->V = solver.eigenvectors().real();
	this->Vi = this->V.inverse();

	this->base_model = base_model;
}

template <class ALPHABET>
ModelFactoryPlusF<ALPHABET>::~ModelFactoryPlusF() {
	delete this->base_model;
}

template <class ALPHABET>
double ModelFactoryPlusF<ALPHABET>::getEpsilon(distance_t distance) const
{
	return this->base_model->getEpsilon(distance);
}

template <class ALPHABET>
double ModelFactoryPlusF<ALPHABET>::getDelta(distance_t distance) const
{
	return this->base_model->getDelta(distance);
}

template <class ALPHABET>
void ModelFactoryPlusF<ALPHABET>::estimateFreqs(const std::vector<sequence_t<ALPHABET> > sequences)
{
	Model<ALPHABET> model = base_model->getModel(1);
	this->Q = model.Q;
	this->freqsOld = model.pi;

	typename Model<ALPHABET>::Freqs freqs;
	freqs = this->freqsOld;
	freqs *= cmdlineopts.pseudo_count;

	for(typename std::vector<sequence_t<ALPHABET> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
		const sequence_t<ALPHABET> &seq = *it;

		for(typename sequence_t<ALPHABET>::const_iterator it2 = seq.begin(); it2 != seq.end(); ++it2) {
			if(it2->isValid()) {
				freqs(it2->value()) += 1;
			}
		}
	}

	this->freqs = freqs/freqs.sum();

	// exchange AA freqs
	this->Q *= (this->freqs.array() / this->freqsOld.array()).matrix().asDiagonal();

	// normalize rate
	this->Q.diagonal().setZero();
	this->Q.diagonal() = -this->Q.rowwise().sum().eval();
	this->Q /= -(this->freqs.transpose() * this->Q.diagonal())(0,0);

	Eigen::EigenSolver<typename Model<ALPHABET>::Subst> solver(this->Q);
	this->sigma = solver.eigenvalues().real();
	this->V = solver.eigenvectors().real();
	this->Vi = this->V.inverse();
}

#endif /* MODELFACTORYPLUSF_H_ */
