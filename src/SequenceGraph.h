#ifndef SEQUENCEGRAPH_H_
#define SEQUENCEGRAPH_H_

#include "Graph.h"
#include "Alphabet.h"
#include "Model.h"
#include "Repeat.h"
#include "CSProfile.h"

template <class ALPHABET>
class SequenceGraph: public Graph<ALPHABET> {
private:
	typedef typename Model<ALPHABET>::Freqs Freqs;
	typedef typename Model<ALPHABET>::Profile Profile;
	typedef typename Graph<ALPHABET>::AdjacencyMatrix AdjacencyMatrix;
	typedef typename Graph<ALPHABET>::TRMatrix TRMatrix;

protected:
	static Freqs char2vec(ALPHABET c);
	static std::vector<Freqs, Eigen::aligned_allocator<Freqs> > nodes2profile(const std::vector<ALPHABET> &nodes, bool startend=false);
	static std::vector<Freqs, Eigen::aligned_allocator<Freqs> > nodes2profile(const sequence_t<ALPHABET> &seq, bool startend=false);

public:
	SequenceGraph() {};
	virtual ~SequenceGraph() {};

	SequenceGraph(const sequence_t<ALPHABET> &seq);
	SequenceGraph(const std::vector<ALPHABET> &nodes);
	SequenceGraph(const sequence_t<AA> &seq, const CSProfile &csprofile, const Model<AA> &model);
	void addNodes(index_t before, const std::vector<ALPHABET> &nodes);
	void addNodesNoEdges(index_t before, const std::vector<ALPHABET> &nodes);
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
typename SequenceGraph<ALPHABET>::Freqs SequenceGraph<ALPHABET>::char2vec(ALPHABET c) {
	Freqs ret;
	if(c.isValid()) {
		ret = Freqs::Unit(c.value());
	} else {
		score_t s = 1.0/ALPHABET::DIM;
		ret = Freqs::Constant(s);
	}
	return ret;
}

template <class ALPHABET>
std::vector<typename SequenceGraph<ALPHABET>::Freqs, Eigen::aligned_allocator<typename SequenceGraph<ALPHABET>::Freqs> > SequenceGraph<ALPHABET>::nodes2profile(const std::vector<ALPHABET> &nodes, bool startend) {
	std::vector<Freqs, Eigen::aligned_allocator<Freqs> > newNodes;
	if(startend) {
		newNodes.reserve(nodes.size() + 2);
	} else {
		newNodes.reserve(nodes.size());
	}

	if(startend) {
		newNodes.push_back(Freqs::Zero());
	}
	for(index_t i = 0; i < nodes.size(); ++i) {
		newNodes.push_back(char2vec(nodes[i]));
	}
	if(startend) {
		newNodes.push_back(Freqs::Zero());
	}

	return newNodes;
}

template <class ALPHABET>
std::vector<typename SequenceGraph<ALPHABET>::Freqs, Eigen::aligned_allocator<typename SequenceGraph<ALPHABET>::Freqs> > SequenceGraph<ALPHABET>::nodes2profile(const sequence_t<ALPHABET> &seq, bool startend) {
	std::vector<Freqs, Eigen::aligned_allocator<Freqs> > newNodes;
	if(startend) {
		newNodes.reserve(seq.size() + 2);
	} else {
		newNodes.reserve(seq.size());
	}

	if(startend) {
		newNodes.push_back(Freqs::Zero());
	}
	for(index_t i = 0; i < seq.size(); ++i) {
		newNodes.push_back(char2vec(seq.at(i)));
	}
	if(startend) {
		newNodes.push_back(Freqs::Zero());
	}

	return newNodes;
}


template <class ALPHABET>
SequenceGraph<ALPHABET>::SequenceGraph(const sequence_t<ALPHABET> &seq)
	: Graph<ALPHABET>(nodes2profile(seq,true))
{ }

template <class ALPHABET>
SequenceGraph<ALPHABET>::SequenceGraph(const std::vector<ALPHABET> &nodes)
	: Graph<ALPHABET>(nodes2profile(nodes,true))
{ }

template <class ALPHABET>
SequenceGraph<ALPHABET>::SequenceGraph(const sequence_t<AA> &seq, const CSProfile &csprofile, const Model<AA> &model) {
	index_t dim = seq.length()+2;
	this->edges = AdjacencyMatrix(dim,dim);
	this->edges.setZero();
	this->repeats = TRMatrix(dim,dim);
	this->repeats.setZero();
	this->sites = csprofile.createProfile(seq,model);

	this->fillInitialEdges();
}

template <class ALPHABET>
void SequenceGraph<ALPHABET>::addNodes(index_t before, const std::vector<ALPHABET> &nodes) {
	this->addNodes(nodes2profile(nodes));
}

template <class ALPHABET>
void SequenceGraph<ALPHABET>::addNodesNoEdges(index_t before, const std::vector<ALPHABET> &nodes) {
	this->addNodesNoEdges(nodes2profile(nodes));
}

#endif /* SEQUENCEGRAPH_H_ */
