#ifndef CLEANEDGRAPH_H_
#define CLEANEDGRAPH_H_

#include "Graph.h"
#include "main.h"

template <class ALPHABET>
class CleanedGraph: public Graph<ALPHABET> {
private:
	typedef typename Graph<ALPHABET>::PredIterator PredIterator;
	typedef typename Graph<ALPHABET>::AdjacencyMatrix AdjacencyMatrix;
	typedef typename Graph<ALPHABET>::TRMatrix TRMatrix;
	typedef typename Model<ALPHABET>::Profile Profile;
	typedef typename Model<ALPHABET>::Freqs Freqs;
	typedef typename Model<ALPHABET>::Subst Subst;

public:
	CleanedGraph(const Graph<ALPHABET> &original);
	virtual ~CleanedGraph();
	index_t getMapping(index_t i) const { return this->outmapping[i]; }
	void uncleanMapping(std::vector<index_t> &mapping) const;

private:
	std::vector<index_t> outmapping;
};



/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

template <class ALPHABET>
CleanedGraph<ALPHABET>::CleanedGraph(const Graph<ALPHABET> &original)
	: Graph<ALPHABET>(original)
{
	std::vector<bool> marked_fw(original.size());
	std::vector<bool> marked_bw(original.size());
	std::vector<index_t> mapping(original.size());

	dp_score_t repeatExt = cmdlineopts.repeatext_prob == 0 ? INFINITY : 0;

	/* keep nodes which are on a path from start to end without max_cost edge */
	for(index_t i=0; i<original.size(); ++i) {
		marked_fw[i] = false;
		marked_bw[i] = false;
		mapping[i] = (index_t)-1;
	}
	marked_fw[0] = true;
	marked_bw[original.size()-1] = true;
	mapping[0] = 0;
	index_t newDim = 1;

	/* mark nodes reachable from end */
	for(index_t i=original.size(); i>0; --i) {
		index_t to = i-1;
		if(!marked_bw[to]) continue;

		for (PredIterator from = original.getPreds(to,0,repeatExt); from; ++from) {
			dp_score_t c = from.value();
			if(c != (dp_score_t)INFINITY) {
				marked_bw[*from] = true;
			}
		}
	}

	/* mark nodes reachable from start */
	for(index_t i=1; i<original.size(); ++i) {
		index_t to = i;

		for (PredIterator from = original.getPreds(to,0,repeatExt); from; ++from) {
			dp_score_t c = from.value();
			if(c != (dp_score_t)INFINITY && marked_fw[*from]) {
				marked_fw[to] = true;
				if(marked_bw[to]) {
					mapping[to] = newDim;
					++newDim;
				}
				break;
			}
		}
	}

	assert(marked_bw[0] && marked_fw[original.size()-1]);

	/* copy over remaining edges */
	std::vector<Eigen::Triplet<dp_score_t> > edges;
	edges.reserve(this->edges.nonZeros());

	for (index_t to = 0; to < original.size(); ++to) {
		for (typename AdjacencyMatrix::InnerIterator from(this->edges,to); from; ++from) {
			index_t y = mapping[to];
			index_t x = mapping[from.index()];
			if(x != (index_t)-1 && y != (index_t)-1 && from.value() < 0) {
				edges.push_back(Eigen::Triplet<dp_score_t>(y,x,from.value()));
			}
		}
	}

	this->edges.setZero();
	this->edges.resize(newDim,newDim);
	this->edges.setFromTriplets(edges.begin(),edges.end());
	this->edges.makeCompressed();
	edges.clear();

	std::vector<Eigen::Triplet<index_t> > tredges;
	tredges.reserve(this->repeats.nonZeros());

	for (index_t to = 0; to < original.size(); ++to) {
		for (typename TRMatrix::InnerIterator from(this->repeats,to); from; ++from) {
			index_t y = mapping[to];
			index_t x = mapping[from.index()];
			if(x != (index_t)-1 && y != (index_t)-1 && from.value() > 0) {
				tredges.push_back(Eigen::Triplet<index_t>(y,x,from.value()));
			}
		}
	}

	this->repeats.setZero();
	this->repeats.resize(newDim,newDim);
	this->repeats.setFromTriplets(tredges.begin(),tredges.end());
	this->repeats.makeCompressed();
	tredges.clear();

	Profile newSites((int)ALPHABET::DIM,newDim);

	for(index_t i=0; i<original.size(); ++i) {
		if(mapping[i] != (index_t)-1) {
			newSites.col(mapping[i]) = original[i];
		}
	}

	this->sites = newSites;

	this->outmapping.resize(this->size());
	for(index_t i=0; i < original.size(); ++i) {
		if(mapping[i] != (index_t)-1) {
			this->outmapping[mapping[i]] = i;
		}
	}
}

template <class ALPHABET>
CleanedGraph<ALPHABET>::~CleanedGraph() {
}

template <class ALPHABET>
void CleanedGraph<ALPHABET>::uncleanMapping(std::vector<index_t> &mapping) const
{
	for(std::vector<index_t>::iterator it = mapping.begin(); it != mapping.end(); ++it) {
		if(*it != (index_t)-1) {
			*it = this->getMapping(*it);
		}
	}
}

#endif /* CLEANEDGRAPH_H_ */
