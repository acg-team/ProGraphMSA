#ifndef MODELFACTORY_H_
#define MODELFACTORY_H_

#include <string>
#include <map>

#include "Model.h"
#include "main.h"

template <class ALPHABET>
class ModelFactory {
private:
	typedef typename Model<ALPHABET>::Profile Profile;
	typedef typename Model<ALPHABET>::Freqs Freqs;
	typedef typename Model<ALPHABET>::Subst Subst;

public:
	virtual Model<ALPHABET> getModel(distance_t distance) const;
	virtual Model<ALPHABET> getModel(distance_t distance, distance_t gap_distance) const;
	virtual ~ModelFactory() {}
	virtual double getEpsilon(distance_t distance) const;
	virtual double getDelta(distance_t distance) const;

	static ModelFactory<ALPHABET> *getDefault(const std::map<std::string,sequence_t<ALPHABET> > &seqs);

protected:
	ModelFactory();
	static void parseDistance(distance_t distance, Model<ALPHABET> &model);
	Freqs freqs;
	Subst Q;
	Subst V;
	Subst Vi;
	Freqs sigma;
};


/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include <cmath>

template <class ALPHABET>
Model<ALPHABET> ModelFactory<ALPHABET>::getModel(distance_t distance) const {
	Model<ALPHABET> model;

	parseDistance(distance,model);

	model.P = this->V * (this->sigma * model.distance).array().exp().matrix().asDiagonal() * this->Vi;
	model.M = this->freqs.asDiagonal() * model.P;

	if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
		model.divergence = 1.0 - model.M.diagonal().sum();
	}
	model.pi = this->freqs;
	model.Q = this->Q;

	model.epsilon = getEpsilon(model.distance);
	model.delta = getDelta(model.distance);

	return model;
}

template <class ALPHABET>
Model<ALPHABET> ModelFactory<ALPHABET>::getModel(distance_t distance, distance_t gap_distance) const {
	Model<ALPHABET> model;

	parseDistance(gap_distance,model);

	model.epsilon = getEpsilon(model.distance);
	model.delta = getDelta(model.distance);

	parseDistance(distance,model);

	model.P = this->V * (this->sigma * model.distance).array().exp().matrix().asDiagonal() * this->Vi;
	model.M = this->freqs.asDiagonal() * model.P;

	if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
		model.divergence = 1.0 - model.M.diagonal().sum();
	}
	model.pi = this->freqs;
	model.Q = this->Q;

	return model;
}

template <class ALPHABET>
double ModelFactory<ALPHABET>::getEpsilon(distance_t distance) const {
	(void)distance;
	return cmdlineopts.gapext_prob;
}

template <class ALPHABET>
double ModelFactory<ALPHABET>::getDelta(distance_t distance) const {
	return (1.0 - std::exp(-distance*cmdlineopts.indel_rate))/2.0;
}

template <class ALPHABET>
void ModelFactory<ALPHABET>::parseDistance(distance_t distance, Model<ALPHABET> &model) {
	distance = std::max(0.0,distance);

	if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
		if(std::isnan(distance)) {
			distance = 5.2;
		}
		model.distance = distance;
		double ed = std::exp(model.distance);
		model.divergence = - 0.5 * (5.0*ed - std::sqrt(45.0*ed*ed - 20.0*ed)) / ed;
	} else {
		if(std::isnan(distance)) {
			distance = 1.0;
		}
		if(distance > 0.85) {
			model.distance = 5.2;
		} else {
			model.distance = - std::log(1.0 - distance - 0.2 * distance * distance);
		}
		model.divergence = distance;
	}
	model.distance = std::max(std::min(model.distance,cmdlineopts.max_dist),cmdlineopts.min_dist);
	model.divergence = std::max(std::min(model.divergence,cmdlineopts.max_pdist),cmdlineopts.min_pdist);
}

template <class ALPHABET>
ModelFactory<ALPHABET>::ModelFactory() {}

#endif /* MODELFACTORY_H_ */
