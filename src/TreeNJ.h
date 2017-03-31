#ifndef TREENJ_H_
#define TREENJ_H_

#include <map>
#include <string>

#include "PhyTree.h"
#include "ModelFactory.h"
#include "Alphabet.h"
#include "LeastSquares.h"

PhyTree* buildNJTree(std::vector<std::string> seqs_order, DistanceMatrix dist, const PhyTree *topo);

/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include "DistanceFactoryAlign.h"
#include "DistanceFactoryAngle.h"
#include "DistanceFactoryPrealigned.h"

template <class ALPHABET>
PhyTree* TreeNJ(const std::map<std::string,sequence_t<ALPHABET> > &seqs, bool prealigned, const ModelFactory<ALPHABET> *model_factory, const PhyTree *topo=NULL)
{
	if(seqs.size() < 2) {
		error("cannot construct tree from < 2 sequences");
	}

	std::vector<std::string> seqs_order(seqs.size());

	index_t i=0;
	for(typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator it = seqs.begin(); it != seqs.end(); ++it, ++i) {
		seqs_order[i] = it->first;
	}

	DistanceFactory<ALPHABET> *dist_factory = DistanceFactory<ALPHABET>::getDefault(model_factory,prealigned);
	DistanceMatrix dist = dist_factory->computePwDistances(seqs,seqs_order);
	delete dist_factory;

	dist.distances.diagonal().setConstant(0);
	dist.variances.diagonal().setConstant(0);

	assert(dist.distances.rows() == (int)seqs.size() && dist.distances.cols() == (int)seqs.size());

	PhyTree *tree = buildNJTree(seqs_order,dist,topo);
	
	if(cmdlineopts.wlsrefine_flag) {
		tree = LeastSquares::refineTree(tree,seqs_order,dist);
	}

	tree = midpointRoot(tree);

	return tree;
}

#endif /* TREENJ_H_ */
