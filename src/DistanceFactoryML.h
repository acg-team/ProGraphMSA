#ifndef DISTANCEFACTORYML_H_
#define DISTANCEFACTORYML_H_

#include "DistanceFactory.h"
#include "ModelFactory.h"
#include "main.h"
#include <Eigen/Dense>

struct distvar_t {
	distvar_t(distance_t dist, distance_t var) {
		this->dist = dist;
		this-> var = var;
	}
	distance_t dist;
	distance_t var;
};

template <class ALPHABET>
class DistanceFactoryML: public DistanceFactory<ALPHABET> {
public:
	DistanceFactoryML(const ModelFactory<ALPHABET> *model_factory);
	virtual ~DistanceFactoryML();
protected:
	typedef Eigen::Matrix<int, ALPHABET::DIM, ALPHABET::DIM> CountMatrix;
	distvar_t computeDistance(const CountMatrix &counts, index_t gaps, double seqlen);
	const ModelFactory<ALPHABET> *model_factory;

private:
	distvar_t computeMLDist(const CountMatrix &counts, index_t gaps, double seqlen, double dist0, double var0);
	static const double DIST_MAX;
	static const double VAR_MAX;
	static const double VAR_MIN;
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
#include <string>
#include <iostream>

#include "Alphabet.h"
#include "debug.h"
#include "ModelFactoryWag.h"

#define MAXITER 20
#define EPSILON 1e-5

template <class ALPHABET>
DistanceFactoryML<ALPHABET>::DistanceFactoryML(const ModelFactory<ALPHABET> *model_factory)
{
	this->model_factory = model_factory;
}

template <class ALPHABET>
DistanceFactoryML<ALPHABET>::~DistanceFactoryML()
{
}

template <class ALPHABET>
distvar_t DistanceFactoryML<ALPHABET>::computeMLDist(const CountMatrix &counts, index_t gaps, double seqlen, double dist0, double var0) {
	double dist_min = 0;
	double dist_max = INFINITY;

	double dist=dist0, var=var0;

	double delta = 1;
	index_t iteration = 0;

	while(std::abs(delta) > EPSILON) {
		if(iteration > MAXITER) {
#ifdef DEBUG
			std::cerr << "maximum number of iterations reached: " << dist << " / " << dist0 << std::endl;
#endif
			if(dist_max == INFINITY) {
				dist = DIST_MAX;
				var = VAR_MAX;
			} else {
				dist = dist0;
				var = var0;
			}

			break;
		}
		Model<ALPHABET> model = this->model_factory->getModel(dist);
		typename Model<ALPHABET>::Subst p = model.P;
		typename Model<ALPHABET>::Subst pp = model.Q * p;
		typename Model<ALPHABET>::Subst ppp = model.Q * pp;

		double f,ff;

		if(cmdlineopts.mldist_gap_flag) {
			double grate = cmdlineopts.indel_rate * seqlen * dist;
			//double g0 = std::pow(grate,gaps) * std::exp(-grate); // XXX / gaps!
			double g = (-grate + gaps)/dist;
			double gg = -gaps/(dist*dist);

			//lf = (p.array().log() * counts.cast<double>().array()).sum() + gaps * std::log(g0);
			f = (counts.template cast<double>().array() * pp.array() / p.array()).sum() + g;
			ff = ((counts.template cast<double>().array() * (ppp.array() * p.array() - pp.array().square())) / p.array().square()).sum() + gg;
		} else {
			//lf = (p.array().log() * counts.cast<double>().array()).sum();
			f = (counts.template cast<double>().array() * pp.array() / p.array()).sum();
			ff = ((counts.template cast<double>().array() * (ppp.array() * p.array() - pp.array().square())) / p.array().square()).sum();
		}

		/* variance estimated by Fisher information */
		var = -1.0/ff;

		/* compute next estimate of distance by either newton or bisection */
		if(f > 0) {
			dist_min = std::max(dist_min,dist);
		} else {
			dist_max = std::min(dist_max,dist);
		}

		double new_dist = dist - f/ff;
		if(!(new_dist < dist_max && new_dist > dist_min)) {
			double upper = (dist_max == INFINITY) ? dist*3 : dist_max;
			double lower = dist_min;
			new_dist = (upper+lower)/2.0;
		}
		delta = 1.0 - new_dist/dist;
		dist = new_dist;

		++iteration;
	}

	return distvar_t(dist,var);
}

template <class ALPHABET>
distvar_t DistanceFactoryML<ALPHABET>::computeDistance(const CountMatrix &counts, index_t gaps, double seqlen)
{
	double ident = counts.diagonal().sum();
	double total = counts.sum();

	double dist0 = 1.0 - ident / total;
	double dist;
	double var;

	if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
		if (total == 0 || dist0 > 0.85) {
			dist = dist0 = DIST_MAX;
			var = VAR_MAX;
		} else {
			dist = dist0 = -std::log(1.0 - dist0 - 0.2 * dist0 * dist0);
			var = dist/total;
		}

		if(total > 0 && ident != total) {
			distvar_t dv = this->computeMLDist(counts, gaps, seqlen, dist, var);
			dist = dv.dist;
			var = dv.var;
		}
	} else {
		if(total == 0) {
			dist = dist0 = 1.0;
			var = VAR_MAX;
		} else {
			dist = dist0;
			var = dist0/total;
		}
	}

	if(!(dist < DIST_MAX)) {
		dist = DIST_MAX;
		var = VAR_MAX;
	}

	if(dist > cmdlineopts.cutoff_dist) {
		dist = cmdlineopts.cutoff_dist;
	}

	if(var < VAR_MIN) {
		var = VAR_MIN;
	}

	if(!(var < VAR_MAX)) {
		var = VAR_MAX;
	}

	return distvar_t(dist,var);
}

#undef EPSILON
#undef MAXITER

#endif /* DISTANCEFACTORYML_H_ */
