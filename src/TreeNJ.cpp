#include "main.h"
#include "debug.h"
#include "TreeNJ.h"

#include <fstream>
#include <vector>
#include <map>
#include <queue>

/* This file implements the BioNJ algorithm described in:
 *
 *   O Gascuel, BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data.
 *   Mol Biol Evol (1997) 14(7): 685-695
 */

#include <Eigen/Dense>

#define MIN_DIST 1e-4
#define MIN_VAR 1e-5
#define INVALID ((index_t)-1)

static double support(double d) {
	double s = 1.0-std::exp(-std::log(2.0)*d/cmdlineopts.edge_halflife);
	s = std::min(1.0,std::max(0.0,s));
	if(std::isnan(s)) {
		s = 0.0;
	}
	return s;
}

static void buidTopoPlanInit(std::map<const PhyTree*,index_t> &n_children_ready, std::map<const PhyTree*,index_t> &n_children_valid, const std::map<std::string,index_t> &orig_leaf_index, std::map<const PhyTree*,index_t> &my_leaf_index, std::vector<bool> &seq_in_tree, const PhyTree *topo) {
	if(!topo->isLeaf()) {
		n_children_ready[topo] = 0;
		n_children_valid[topo] = 0;
		for(std::vector<PhyTree*>::const_iterator child=topo->firstChild(); child != topo->lastChild(); ++child) {
			buidTopoPlanInit(n_children_ready, n_children_valid, orig_leaf_index, my_leaf_index, seq_in_tree, *child);
			n_children_valid[topo] += n_children_valid[*child];
			n_children_ready[topo] += (*child)->isLeaf();
		}
	} else {
		std::map<std::string,index_t>::const_iterator pos = orig_leaf_index.find(topo->getName());
		if(pos != orig_leaf_index.end()) {
			my_leaf_index[topo] = pos->second;
			n_children_valid[topo] = 1;
			seq_in_tree[pos->second] = true;
		} else {
			my_leaf_index[topo] = INVALID;
			n_children_valid[topo] = 0;
		}
	}
}

static std::queue<std::pair<index_t,index_t> > buildTopoPlan(const std::vector<std::string> &seqs_order, const PhyTree *topo) {
	std::map<const PhyTree*,index_t> n_children_ready;
	std::map<const PhyTree*,index_t> n_children_valid;
	std::map<std::string,index_t> orig_leaf_index;
	std::map<const PhyTree*,index_t> my_leaf_index;
	std::vector<bool> seq_in_tree(seqs_order.size());
	std::queue<std::pair<index_t,index_t> > topo_plan;
	std::queue<const PhyTree*> worklist;

	for(index_t i=0; i<seqs_order.size(); ++i) {
		orig_leaf_index[seqs_order[i]] = i;
		seq_in_tree[i] = false;
	}

	buidTopoPlanInit(n_children_ready,n_children_valid,orig_leaf_index,my_leaf_index,seq_in_tree,topo);

	for(index_t i=0; i<seq_in_tree.size(); ++i) {
		if(!seq_in_tree[i]) {
			error("sequence \"%s\"is missing in given topology",seqs_order[i].c_str());
		}
	}

	for(std::map<const PhyTree*,index_t>::const_iterator it=n_children_ready.begin(); it != n_children_ready.end(); ++it) {
		if(it->second == it->first->n_children()) {
			worklist.push(it->first);
		}
	}

	while(!worklist.empty()) {
		const PhyTree *node = worklist.front();
		worklist.pop();
		assert(node->n_children() == 2);

		const PhyTree &child1 = (*node)[0];
		const PhyTree &child2 = (*node)[1];

		index_t index1 = my_leaf_index.at(&child1);
		index_t index2 = my_leaf_index.at(&child2);

		my_leaf_index.erase(&child1);
		my_leaf_index.erase(&child2);

		assert(child1.getParent() == node && child2.getParent() == node);

		if(index1 == INVALID) {
			my_leaf_index[node] = index2;
		} else if(index2 == INVALID) {
			my_leaf_index[node] = index1;
		} else {
			if(index1 > index2) {
				index_t h = index1;
				index1 = index2;
				index2 = h;
			}

			assert(index1 != INVALID && index2 != INVALID);

			my_leaf_index[node] = index1;
			topo_plan.push(std::pair<index_t,index_t>(index1,index2));

			for(std::map<const PhyTree*,index_t>::iterator it=my_leaf_index.begin(); it != my_leaf_index.end(); ++it) {
				if(it->second > index2 && it->second != INVALID) {
					it->second -= 1;
				}
			}
		}

		const PhyTree *parent = node->getParent();
		if(parent) {
			n_children_ready[parent] += 1;
			if(n_children_ready[parent] == parent->n_children()) {
				worklist.push(parent);
			}
		}
	}

	return topo_plan;
}

PhyTree* buildNJTree(std::vector<std::string> seqs_order, DistanceMatrix dist, const PhyTree *topo)
{
	std::queue<std::pair<index_t,index_t> > topo_plan;
	if(topo != NULL) {
		topo_plan = buildTopoPlan(seqs_order,topo);
	}

	std::vector<PhyTree*> subtrees_order(seqs_order.size());
	index_t i = 0;
	for (std::vector<std::string>::const_iterator it = seqs_order.begin(); it != seqs_order.end(); ++it, ++i) {
		subtrees_order[i] = new PhyTree(*it);
	}

	for (index_t dim = seqs_order.size(); dim > 3; --dim) {
		int index1;
		int index2;
		double min = 0;

		dist.distances = dist.distances.array().max(MIN_DIST * Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(dim,dim).array());
		dist.variances = dist.variances.array().max(MIN_VAR * Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(dim,dim).array());

		dist.distances.diagonal().setZero();
		dist.variances.diagonal().setZero();

		Eigen::Matrix<double, 1, Eigen::Dynamic> sums = dist.distances.colwise().sum();

		if(topo_plan.empty()) {
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Q;
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(dim, dim) * sums.asDiagonal();

			Q = 0.5 * dist.distances - (0.5 / (dim - 2.0)) * (S + S.transpose());
			assert(Q.isApprox(Q.transpose(),1e-6) && "Q has to be symmetrical -> problem with libEigen (using buggy gcc-4.4?)!");
			Q.diagonal().setConstant(INFINITY);

			min = Q.minCoeff(&index2, &index1);
			if (index2 < index1) {
				int h = index2;
				index2 = index1;
				index1 = h;
			}

			assert(min < INFINITY);
			assert(index1 < index2);
		} else {
			index1 = topo_plan.front().first;
			index2 = topo_plan.front().second;
			topo_plan.pop();
		}
		(void) min;

		std::string name1 = seqs_order[index1];
		std::string name2 = seqs_order[index2];

		double dist1 = (dist.distances(index1, index2) + (sums(index1) - sums(index2)) / (dim - 2.0)) / 2.0;
		dist1 = std::min(std::max(dist1, MIN_DIST), dist.distances(index1, index2));
		double dist2 = std::max(dist.distances(index2, index1) - dist1, MIN_DIST);

		// reduce distance matrix
		DistanceMatrix reduced_dist(dim - 1);
		// NOTE: index2 is deleted and new node gets index1
		reduced_dist = dist.reduce(index2);

		// should never happen
		if (index1 == index2) {
#ifndef NDEBUG
			std::cerr << min << "," << index1 << "," << index2 << "," << (dim - 1) << std::endl;
#endif
			index1 = 0;
			index2 = 1;
		}

		// distances and variances to new node
#if 0
		// classical NJ
		double lambda = .5;
#else
		// BioNJ
		double lambda = .5 + (dist.variances.row(index2) - dist.variances.row(index1)).sum() / (2 * (dim - 2) * dist.variances(index1, index2));

		if (std::isnan(lambda))
			lambda = .5;
		else
			lambda = std::min(std::max(0.0, lambda), 1.0);
#endif

		Eigen::Matrix<double, 1, Eigen::Dynamic> new_dist(dim - 1);
		Eigen::Matrix<double, 1, Eigen::Dynamic> new_var(dim - 1);
		for (index_t i = 0; i < (index_t) index2; ++i) {
			new_dist(i) = lambda * (dist.distances(index1, i) - dist1) + (1.0
					- lambda) * (dist.distances(index2, i) - dist2);
			new_var(i) = lambda * dist.variances(index1, i)
					+ (1.0 - lambda) * dist.variances(index2, i)
					- lambda * (1.0 - lambda) * dist.variances(index1, index2);
		}

		new_dist(index1) = 0;
		new_var(index1) = 0;

		for (index_t i = index2; i < dim - 1; ++i) {
			new_dist(i) = lambda * (dist.distances(index1, i + 1) - dist1)
					+ (1.0 - lambda) * (dist.distances(index2, i + 1) - dist2);
			new_var(i) = lambda * dist.variances(index1, i + 1)
					+ (1.0 - lambda) * dist.variances(index2, i + 1)
					- lambda * (1.0 - lambda) * dist.variances(index1, index2);
		}

		reduced_dist.distances.block(index1, 0, 1, dim - 1) = new_dist;
		reduced_dist.distances.block(0, index1, dim - 1, 1)	= new_dist.transpose();

		reduced_dist.variances.block(index1, 0, 1, dim - 1) = new_var;
		reduced_dist.variances.block(0, index1, dim - 1, 1) = new_var.transpose();

		dist = reduced_dist;

		// adjust sequence list
		seqs_order.erase(seqs_order.begin() + index2);
		seqs_order[index1] = name1 + "," + name2;

		// create new node
		PhyTree *tree = new PhyTree(seqs_order[index1]);
		tree->addChild(subtrees_order[index1], dist1, support(dist1));
		tree->addChild(subtrees_order[index2], dist2, support(dist2));
		subtrees_order.erase(subtrees_order.begin() + index2);
		subtrees_order[index1] = tree;
	}

	PhyTree *tree = new PhyTree("root");
	if (seqs_order.size() == 2) {
		double d = dist.distances(0, 1) / 2;
		tree->addChild(subtrees_order[0], d, support(d));
		tree->addChild(subtrees_order[1], d, support(d));
	} else {
		assert(seqs_order.size() == 3);

		distance_t d0 = (dist.distances(0, 1) + dist.distances(0, 2) - dist.distances(1, 2)) / 2.0;
		d0 = std::min(std::max(d0, MIN_DIST), std::min(dist.distances(1, 0), dist.distances(2, 0)));
		distance_t d1 = std::max(dist.distances(1, 0) - d0, MIN_DIST);
		distance_t d2 = std::max(dist.distances(2, 0) - d0, MIN_DIST);

		PhyTree *tree2 = new PhyTree("root2");

		tree2->addChild(subtrees_order[0], d0, support(d0));
		tree2->addChild(subtrees_order[1], d1, support(d1));

		tree->addChild(subtrees_order[2], d2 / 2, support(d2));
		tree->addChild(tree2, d2 / 2, support(d2));
	}

	return tree;
}
