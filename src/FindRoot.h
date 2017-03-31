#ifndef FINDROOT_H_
#define FINDROOT_H_

#include "PhyTree.h"
#include "ProgressiveAlignment.h"
#include "RepeatDetectionTReks.h"
#include "CSProfile.h"
#include "Alphabet.h"
#include "debug.h"
#include "main.h"


namespace FindRoot {
	template <class ALPHABET>
	struct edge;

	template <class ALPHABET>
	struct node;

	template <class ALPHABET>
	struct node {
		edge<ALPHABET> *edges[3];
		ProgressiveAlignmentResult<ALPHABET>* cached_alignments[3];
		std::string name;

		node() {
			this->edges[0] = this->edges[1] = this->edges[2] = NULL;
			this->cached_alignments[0] = this->cached_alignments[1] = this->cached_alignments[2] = NULL;
			this->name = "";
		}

		~node() {
			if(this->cached_alignments[0]) delete this->cached_alignments[0];
			if(this->cached_alignments[1]) delete this->cached_alignments[1];
			if(this->cached_alignments[2]) delete this->cached_alignments[2];
		}

		node *getAdjacentNode(index_t index);

		bool isLeaf() const {
			return this->edges[1] == NULL;
		}

		ProgressiveAlignmentResult<ALPHABET> *getAlignment(edge<ALPHABET> *e, const ModelFactory<ALPHABET> &model_factory);
	};

	template <class ALPHABET>
	struct edge {
		ProgressiveAlignmentResult<ALPHABET>* alignment;
		node<ALPHABET> *nodes[2];
		double length;
		double support;

		edge() {
			this->alignment = NULL;
			this->length = -1;
			this->support = -1;
		}

		~edge() {
			if(this->alignment) delete this->alignment;
		}

		bool visited() const { return this->alignment != NULL; }

		ProgressiveAlignmentResult<ALPHABET> *getAlignment(const ModelFactory<ALPHABET> &model_factory);
	};
}

/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include "GapParsimony.h"
#include "SequenceGraph.h"

template <class ALPHABET>
ProgressiveAlignmentResult<ALPHABET> progressive_alignment_find_root(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const PhyTree &tree, const std::map<std::string, std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<ALPHABET> &model_factory)  throw (ProgressiveAlignmentException);

namespace FindRoot {
	template <class ALPHABET>
	static node<ALPHABET> *getOtherNode(const edge<ALPHABET>* e, const node<ALPHABET>* n);

	template <class ALPHABET>
	static double getEdgeLength(const edge<ALPHABET>* e); // FIXME sorry for this ugly C++ hack

	template <class ALPHABET>
	ProgressiveAlignmentResult<ALPHABET> *node<ALPHABET>::getAlignment(edge<ALPHABET> *e, const ModelFactory<ALPHABET> &model_factory) {
		int index = -1;
		if(e == this->edges[0]) {
			index = 0;
		} else if (e == this->edges[1]) {
			index = 1;
		} else if (e == this->edges[2]) {
			index = 2;
		} else {
			error("invalid edge pointer");
			return NULL;
		}

		// is result already cached?
		if(this->cached_alignments[index] == NULL) {
			assert(!this->isLeaf() && "leaf nodes should at this point already have a result attached");

			this->cached_alignments[index] = new ProgressiveAlignmentResult<ALPHABET>();

			int index1 = -1;
			int index2 = -1;
			switch(index) {
			case 0:
				index1 = 1; index2 = 2; break;
			case 1:
				index1 = 0; index2 = 2; break;
			case 2: default:
				index1 = 0; index2 = 1; break;
			}

			edge<ALPHABET> *e1 = this->edges[index1];
			edge<ALPHABET> *e2 = this->edges[index2];

			ProgressiveAlignmentResult<ALPHABET> *r1 = this->getAdjacentNode(index1)->getAlignment(e1, model_factory);
			ProgressiveAlignmentResult<ALPHABET> *r2 = this->getAdjacentNode(index2)->getAlignment(e2, model_factory);

			*(this->cached_alignments[index]) = align_progressive_results(*r1,*r2,getEdgeLength(e1),getEdgeLength(e2),getEdgeSupport(e1),getEdgeSupport(e2),model_factory);
		}

		return this->cached_alignments[index];
	}

	template <class ALPHABET>
	ProgressiveAlignmentResult<ALPHABET> *edge<ALPHABET>::getAlignment(const ModelFactory<ALPHABET> &model_factory) {
		if(this->alignment) return this->alignment;

		this->alignment = new ProgressiveAlignmentResult<ALPHABET>();

		ProgressiveAlignmentResult<ALPHABET> *r1 = this->nodes[0]->getAlignment(this, model_factory);
		ProgressiveAlignmentResult<ALPHABET> *r2 = this->nodes[1]->getAlignment(this, model_factory);

		*(this->alignment) = align_progressive_results(*r1,*r2,this->length/2,this->length/2,this->support,this->support,model_factory);

		return this->alignment;
	}

	template <class ALPHABET>
	node<ALPHABET> *node<ALPHABET>::getAdjacentNode(index_t index) {
		return getOtherNode(this->edges[index], this);
	}

	template <class ALPHABET>
	static node<ALPHABET> *getOtherNode(const edge<ALPHABET> *e, const node<ALPHABET> *n) {
		return reinterpret_cast<node<ALPHABET> *>(reinterpret_cast<long>(e->nodes[1]) ^ reinterpret_cast<long>(e->nodes[0]) ^ reinterpret_cast<long>(n));
	}

	template <class ALPHABET>
	static double getEdgeLength(const edge<ALPHABET> *e) {
		return e->length;
	}

	template <class ALPHABET>
	static double getEdgeSupport(const edge<ALPHABET> *e) {
		return e->support;
	}

	template <class ALPHABET>
	static void tree2graph(FindRoot::node<ALPHABET> *current, const PhyTree& tree, std::vector<FindRoot::node<ALPHABET>*> &nodes, std::vector<FindRoot::edge<ALPHABET>*> &edges, const std::map<std::string,sequence_t<ALPHABET> > &sequences, const std::map<std::string, std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<ALPHABET> &model_factory) {
		(void)csprofile;

		current->name = tree.getName();

		if(tree.isLeaf()) {
			typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator it = sequences.find(current->name);
			typename std::map<std::string, std::vector<repeat_t> >::const_iterator it2 = repeats.find(current->name);

			if(it != sequences.end()) {
				current->cached_alignments[0] = new ProgressiveAlignmentResult<ALPHABET>();
				current->cached_alignments[0]->is_csprofile = false;
				current->cached_alignments[0]->graph = SequenceGraph<ALPHABET>(it->second);

				if(it2 != repeats.end()) {
					for(std::vector<repeat_t>::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3) {
						std::vector<int> tr_hom(current->cached_alignments[0]->graph.size(),-1);
						std::copy(it3->tr_hom.begin(), it3->tr_hom.end(), tr_hom.begin()+it3->start+1);
						current->cached_alignments[0]->tr_homologies.push_back(tr_hom);
						current->cached_alignments[0]->tr_source.push_back(current->name);
					}
					current->cached_alignments[0]->graph.addRepeats(current->cached_alignments[0]->tr_homologies);
				}

				current->cached_alignments[0]->aligned_sequences[current->name] = it->second;
				current->cached_alignments[0]->score = 0;
			} else {
				error("unknown sequence name: %s",current->name.c_str());
			}
		} else {
			assert(tree.n_children() == 2);

			FindRoot::edge<ALPHABET> *e0 = new FindRoot::edge<ALPHABET>();
			edges.push_back(e0);
			e0->length = tree[0].getBranchLength();
			e0->support = tree[0].getBranchSupport();
			current->edges[1] = e0;

			FindRoot::node<ALPHABET> *n0 = new FindRoot::node<ALPHABET>();
			nodes.push_back(n0);
			n0->edges[0] = e0;
			e0->nodes[0] = current;
			e0->nodes[1] = n0;

			tree2graph(n0,tree[0],nodes,edges,sequences,repeats,csprofile,model_factory);


			FindRoot::edge<ALPHABET> *e1 = new FindRoot::edge<ALPHABET>();
			edges.push_back(e1);
			e1->length = tree[1].getBranchLength();
			e1->support = tree[1].getBranchSupport();
			current->edges[2] = e1;

			FindRoot::node<ALPHABET> *n1 = new FindRoot::node<ALPHABET>();
			nodes.push_back(n1);
			n1->edges[0] = e1;
			e1->nodes[0] = current;
			e1->nodes[1] = n1;

			tree2graph(n1,tree[1],nodes,edges,sequences,repeats,csprofile,model_factory);
		}
	}
}

template <class ALPHABET>
ProgressiveAlignmentResult<ALPHABET> progressive_alignment_find_root(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const PhyTree &tree, const std::map<std::string, std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<ALPHABET> &model_factory)  throw (ProgressiveAlignmentException)
{
	// transform tree into graph
	std::vector<FindRoot::node<ALPHABET>*> nodes;
	std::vector<FindRoot::edge<ALPHABET>*> edges;

	if (tree.n_children() == 2) {
		FindRoot::edge<ALPHABET> *e0 = new FindRoot::edge<ALPHABET>();
		edges.push_back(e0);
		e0->length = tree[0].getBranchLength() + tree[1].getBranchLength();
		e0->support = std::max(tree[0].getBranchSupport(),tree[1].getBranchSupport());
		FindRoot::node<ALPHABET> *n0 = new FindRoot::node<ALPHABET>();
		nodes.push_back(n0);
		n0->edges[0] = e0;
		e0->nodes[0] = n0;
		FindRoot::node<ALPHABET> *n1 = new FindRoot::node<ALPHABET>();
		nodes.push_back(n1);
		n1->edges[0] = e0;
		e0->nodes[1] = n1;

		FindRoot::tree2graph<ALPHABET>(n0,tree[0],nodes,edges,sequences,repeats,csprofile,model_factory);
		FindRoot::tree2graph<ALPHABET>(n1,tree[1],nodes,edges,sequences,repeats,csprofile,model_factory);
	} else if (tree.n_children() == 3) {
		FindRoot::node<ALPHABET> *n0 = new FindRoot::node<ALPHABET>();
		nodes.push_back(n0);

		for(int i=0; i<3; ++i) {
			FindRoot::edge<ALPHABET> *ei = new FindRoot::edge<ALPHABET>();
			edges.push_back(ei);
			ei->length = tree[i].getBranchLength();
			ei->support = tree[i].getBranchSupport();
			FindRoot::node<ALPHABET> *ni = new FindRoot::node<ALPHABET>();
			nodes.push_back(ni);
			ni->edges[0] = ei;
			ei->nodes[0] = n0;
			ei->nodes[1] = ni;
			n0->edges[i] = ei;
			FindRoot::tree2graph<ALPHABET>(ni,tree[i],nodes,edges,sequences,repeats,csprofile,model_factory);
		}
	} else {
		error("multifurcations not allowed");
	}

	ProgressiveAlignmentResult<ALPHABET> *best_result = edges[0]->getAlignment(model_factory);
	score_t best_score = GapParsimony::scoreAlignment<ALPHABET>(*best_result,edges[0]);

	if(cmdlineopts.reroot_flag == 1) {
		for(typename std::vector<FindRoot::edge<ALPHABET>*>::iterator it = ++edges.begin(); it != edges.end(); ++it) {
			ProgressiveAlignmentResult<ALPHABET> *result = (*it)->getAlignment(model_factory);
			score_t score = GapParsimony::scoreAlignment<ALPHABET>(*result,*it);
			if(score < best_score) {
				best_result = result;
				best_score = score;
			}
		}
	} else {
		/* heuristic try all adjacent edges as roots (leaving out the already marked) */
		FindRoot::node<ALPHABET> *best_node = NULL;
		FindRoot::edge<ALPHABET> *best_edge = edges[0];

		do {
			FindRoot::edge<ALPHABET> *old_edge = best_edge;
			FindRoot::node<ALPHABET> *old_node = best_node;

			for(index_t i=0; i<2; ++i) {
				FindRoot::node<ALPHABET> *n = old_edge->nodes[i];
				if(n == old_node) continue;

				for(index_t j=0; j<3; ++j) {
					FindRoot::edge<ALPHABET> *e = n->edges[j];
					if(e == old_edge || e == NULL) continue;

					ProgressiveAlignmentResult<ALPHABET> *result = e->getAlignment(model_factory);
					score_t score = GapParsimony::scoreAlignment<ALPHABET>(*result,e);
					if(score < best_score) {
						best_result = result;
						best_edge = e;
						best_score = score;
						best_node = n;
					}
				}
			}

			if(best_edge == old_edge) break;
		} while (1);
	}

	std::cerr << "best gap parsimony score: " << best_score << std::endl;

	ProgressiveAlignmentResult<ALPHABET> result = *best_result;

	// TODO discard cached alignments if running out of memory

	for(typename std::vector<FindRoot::node<ALPHABET>*>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
		delete *it;
	}

	for(typename std::vector<FindRoot::edge<ALPHABET>*>::iterator it = edges.begin(); it != edges.end(); ++it) {
		delete *it;
	}

	return result;
}

#endif /* FINDROOT_H_ */
