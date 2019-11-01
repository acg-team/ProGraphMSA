#include "main.h"
#include "LeastSquares.h"
#include "PhyTree.h"
#include "DistanceFactory.h"
#include "debug.h"
#include "NNLS.h"
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <iostream>

#define MAX_ITERS 20
#define MAX_ITERS5 5

static double support(double d) {
	double s = 1.0-std::exp(-std::log(2.0)*d/cmdlineopts.edge_halflife);
	s = std::min(1.0,std::max(0.0,s));
	if(std::isnan(s)) {
		s = 0.0;
	}
	return s;
}

namespace LeastSquares {

struct Node;

struct Edge {
	Node   *nodes[2];
	distance_t  length;
	double  support;

	Node& other(const Node &node) {
		return *(reinterpret_cast<Node*>(reinterpret_cast<long>(nodes[0]) ^ reinterpret_cast<long>(nodes[1]) ^ reinterpret_cast<long>(&node)));
	}
	Node& operator[](index_t i) {
		return *(nodes[i]);
	}
	const Node& other(const Node &node) const {
		return *(reinterpret_cast<const Node*>(reinterpret_cast<long>(nodes[0]) ^ reinterpret_cast<long>(nodes[1]) ^ reinterpret_cast<long>(&node)));
	}
	const Node& operator[](index_t i) const {
		return *(nodes[i]);
	}
};

struct Node {
	Edge *edges[3];
	index_t leaf;
	bool todo;

	bool isLeaf() const {
		return leaf != (index_t)-1;
	}
	Edge& operator[](index_t i) {
		return *(edges[i]);
	}
	const Edge& operator[](index_t i) const {
		return *(edges[i]);
	}
};

struct Graph {
	Graph(index_t leaves) {
		n_leaves = leaves;
		n_nodes = 2*n_leaves-2;
		n_edges = 2*n_leaves-3;
		labels = new std::string[n_leaves];
		nodes = new Node[n_nodes];		
		edges = new Edge[n_edges];
	}

	Graph(const PhyTree &tree, const std::vector<std::string>& leaves_order) {
		assert(tree.n_children() == 2);
		n_leaves = leaves_order.size();

		labels = new std::string[n_leaves];
		nodes = new Node[2*n_leaves-2];		
		edges = new Edge[2*n_leaves-3];

		for(index_t i=0; i < this->n_leaves; ++i) {
			this->labels[i] = leaves_order[i];
		}

		n_nodes = 0;
		n_edges = 1;
		this->edges[0].length =	tree[0].getBranchLength() + tree[1].getBranchLength();
		this->edges[0].nodes[0] = this->tree2graphR(tree[0],&this->edges[0]);
		this->edges[0].nodes[1] = this->tree2graphR(tree[1],&this->edges[0]);
		assert(n_nodes == 2*n_leaves - 2);
		assert(n_edges == 2*n_leaves - 3);
	}

	Node *nodes;
	Edge *edges;
	std::string *labels;
	index_t n_leaves;
	index_t n_nodes;
	index_t n_edges;

	const std::string& getLabel(const Node* node) const {
		static std::string empty("");
		if(node->leaf == (index_t)-1) {
			return empty;
		} else {
			return this->labels[node->leaf];
		}
	}

	std::map<const Node*,distance_t> subtreeDist(const Node *node, const Edge *from) const {
		std::map<const Node*,distance_t> dists;
		subtreeDistR(node,from,dists,0);
		return dists;
	}

	PhyTree* toTree() const {
		const Edge& e = this->edges[0];
		PhyTree *root = new PhyTree();
		root->addChild(this->toTreeR(&e[0],&e),e.length/2.0,e.support);
		root->addChild(this->toTreeR(&e[1],&e),e.length/2.0,e.support);
		return root;
	}

	void verify() const {
#ifdef DEBUG
		for(index_t i=0; i < this->n_edges; ++i) {
			const Edge& e = this->edges[i];
			const Node& n0 = e[0];
			const Node& n1 = e[1];
			assert(&n0 != &n1);
			assert(&n1 == &e.other(n0));
			assert(&n0[0] == &e || &n0[1] == &e || &n0[2] == &e);
			assert(&n1[0] == &e || &n1[1] == &e || &n1[2] == &e);
		}
#endif
	}

private:
	Node* tree2graphR(const PhyTree &tree, Edge* edge) {
		Node& node = this->nodes[this->n_nodes++];
		node.leaf = (index_t)-1;
		if(tree.isLeaf()) {
			node.edges[0] = edge;
			node.edges[1] = node.edges[2] = NULL;

			for(index_t i=0; i<this->n_leaves; ++i) {
				if(tree.getName() == this->labels[i]) {
					node.leaf = i;
					break;
				}
			}
			if(node.leaf == (index_t)-1) {
				error("unknown leaf name: %s",tree.getName().c_str());
			}
		} else {
			node.edges[0] = edge;

			node.edges[1] = &this->edges[this->n_edges++];
			node[1].length = tree[0].getBranchLength();
			node[1].nodes[0] = &node;
			node[1].nodes[1] = this->tree2graphR(tree[0],node.edges[1]);

			node.edges[2] = &this->edges[this->n_edges++];
			node[2].length = tree[1].getBranchLength();
			node[2].nodes[0] = &node;
			node[2].nodes[1] = this->tree2graphR(tree[1],node.edges[2]);
		}
		return &node;
	}

	static void subtreeDistR(const Node *node, const Edge *from, std::map<const Node*,distance_t> &dists, distance_t dist) {
		if(node->isLeaf()) {
			dists[node] = dist;
		} else {
			for(index_t i=0; i<3; ++i) {
				const Edge *e = &(*node)[i];
				if(e != from) {
					subtreeDistR(&e->other(*node),e,dists,dist+e->length);
				}
			}
		}
	}

	PhyTree* toTreeR(const Node *node, const Edge *from) const {
		PhyTree *tree;
		if(!node->isLeaf()) {
			tree = new PhyTree();
			for(index_t i=0; i<3; ++i) {
				const Edge *e = node->edges[i];
				if(e != from) {
					tree->addChild(this->toTreeR(&e->other(*node),e),e->length,e->support);
				}
			}
		} else {
			tree = new PhyTree(this->getLabel(node));
		}
		return tree;
	}
};

double computeFit(const Graph &g, const DistanceMatrix &dist) {
	double fit = 0.0;

	for(index_t i=0; i < g.n_nodes; ++i) {
		Node &n = g.nodes[i];
		if(!n.isLeaf()) continue;
		index_t i1 = n.leaf;

		std::map<const Node*,distance_t> leaves = g.subtreeDist(&n[0].other(n),&n[0]);
		for(std::map<const Node*,distance_t>::iterator it=leaves.begin(); it != leaves.end(); ++it) {
			index_t i2 = it->first->leaf;
			fit += (it->second + n[0].length - dist.distances(i1,i2))*dist.variances(i1,i2)*(it->second + n[0].length - dist.distances(i1,i2));
		}
	}
	return fit;
}

static double Opt4(const index_t leaf_mapping[4], Node* nodes[6], Edge* edges[5], const DistanceMatrix &dist, double &best_fit, bool apply) {
	static const double data[] = {
		1, 1, 0, 0, 0,
		1, 0, 1, 0, 1,
		1, 0, 0, 1, 1,
		0, 1, 1, 0, 1,
		0, 1, 0, 1, 1,
		0, 0, 1, 1, 0};

	Eigen::Matrix<double,6,1> dists;
	Eigen::Matrix<double,5,1> old_dists;
	Eigen::Matrix<double,5,1> new_dists;
	Eigen::Matrix<double,6,1> weights;
	Eigen::Matrix<double,6,5> A;

	dists[0] = dist.distances(leaf_mapping[0],leaf_mapping[1]);
	dists[1] = dist.distances(leaf_mapping[0],leaf_mapping[2]);
	dists[2] = dist.distances(leaf_mapping[0],leaf_mapping[3]);
	dists[3] = dist.distances(leaf_mapping[1],leaf_mapping[2]);
	dists[4] = dist.distances(leaf_mapping[1],leaf_mapping[3]);
	dists[5] = dist.distances(leaf_mapping[2],leaf_mapping[3]);

	weights[0] = dist.variances(leaf_mapping[0],leaf_mapping[1]);
	weights[1] = dist.variances(leaf_mapping[0],leaf_mapping[2]);
	weights[2] = dist.variances(leaf_mapping[0],leaf_mapping[3]);
	weights[3] = dist.variances(leaf_mapping[1],leaf_mapping[2]);
	weights[4] = dist.variances(leaf_mapping[1],leaf_mapping[3]);
	weights[5] = dist.variances(leaf_mapping[2],leaf_mapping[3]);

	A.noalias() = weights.asDiagonal() * Eigen::Map<const Eigen::Matrix<double,6,5,Eigen::RowMajor> >(data);

	new_dists = NNLS(A,dists);

	double fit = (A*new_dists - dists).squaredNorm();

	if(fit < best_fit && apply) {
		best_fit = fit;
		for(index_t i=0; i<4; ++i) {
			edges[leaf_mapping[i]]->length = new_dists(i);
			edges[leaf_mapping[i]]->nodes[0] = nodes[leaf_mapping[i]];
		}
		edges[4]->length = new_dists(4);
		edges[4]->nodes[0] = nodes[4];
		edges[4]->nodes[1] = nodes[5];

		nodes[4]->edges[0] = edges[4];
		nodes[4]->edges[1] = edges[leaf_mapping[0]];
		nodes[4]->edges[2] = edges[leaf_mapping[1]];

		nodes[5]->edges[0] = edges[4];
		nodes[5]->edges[1] = edges[leaf_mapping[2]];
		nodes[5]->edges[2] = edges[leaf_mapping[3]];

		edges[leaf_mapping[0]]->nodes[1] = nodes[4];
		edges[leaf_mapping[1]]->nodes[1] = nodes[4];
		edges[leaf_mapping[2]]->nodes[1] = nodes[5];
		edges[leaf_mapping[3]]->nodes[1] = nodes[5];
	}

	return fit;
}

static void OptimizeQuartet(Edge *e, Graph *g, const DistanceMatrix &all_dist, bool apply)
{
	double best_fit = INFINITY;
	Node *nodes[6];
	Edge *edges[5];

	if((*e)[0].isLeaf() || (*e)[1].isLeaf()) {
		e->support = support(e->length);
		return;
	}

	edges[4] = e;
	nodes[4] = &(*e)[0];
	nodes[5] = &(*e)[1];
	edges[0] = &(*nodes[4])[0] == e ? &(*nodes[4])[1] : &(*nodes[4])[0];
	edges[1] = &(*nodes[4])[2] == e ? &(*nodes[4])[1] : &(*nodes[4])[2];
	edges[2] = &(*nodes[5])[0] == e ? &(*nodes[5])[1] : &(*nodes[5])[0];
	edges[3] = &(*nodes[5])[2] == e ? &(*nodes[5])[1] : &(*nodes[5])[2];
	nodes[0] = &edges[0]->other(*nodes[4]);
	nodes[1] = &edges[1]->other(*nodes[4]);
	nodes[2] = &edges[2]->other(*nodes[5]);
	nodes[3] = &edges[3]->other(*nodes[5]);


	std::map<const Node*,distance_t> leaf_distances[4];
	for(index_t i=0; i<4; ++i) {
		leaf_distances[i] = g->subtreeDist(nodes[i],edges[i]);
	}
	DistanceMatrix dist(4);
	for(index_t i=0; i<4; ++i) {
		for(index_t j=i+1; j<4; ++j) {
			dist.distances(i,j) = 0;
			dist.variances(i,j) = 0;
			for(std::map<const Node*,distance_t>::iterator k = leaf_distances[i].begin(); k != leaf_distances[i].end(); ++k) {
				for(std::map<const Node*,distance_t>::iterator l = leaf_distances[j].begin(); l != leaf_distances[j].end(); ++l) {
					index_t kk = k->first->leaf;
					index_t ll = l->first->leaf;
					dist.distances(i,j) += all_dist.variances(kk,ll) * (all_dist.distances(kk,ll) - k->second - l->second);
					dist.variances(i,j) += all_dist.variances(kk,ll);
				}
			}
			dist.distances(j,i) = dist.distances(i,j);
			dist.variances(j,i) = dist.variances(i,j);
		}
	}
	dist.variances = dist.variances.array().sqrt();
	dist.distances.array() /= dist.variances.array();

	static index_t mapping1[4] = {0,1,2,3};
	double f1 = Opt4(mapping1,nodes,edges,dist,best_fit,apply);

	static index_t mapping2[4] = {0,2,1,3};
	double f2 = Opt4(mapping2,nodes,edges,dist,best_fit,apply);

	static index_t mapping3[4] = {0,3,1,2};
	double f3 = Opt4(mapping3,nodes,edges,dist,best_fit,apply);

	e->support = 1.0 / (1.0 + std::exp((f2-f1)/-2.0) + std::exp((f3-f1)/-2.0));
}

static void OptimizeQuartets(Graph *g, const DistanceMatrix &all_dist, bool apply)
{
	for(index_t i=0; i < g->n_edges; ++i) {
		Edge *e = &g->edges[i];
		OptimizeQuartet(e,g,all_dist,apply);
		g->verify();
	}
}

// Node 4 at branch to 0
static double Opt5v1(const index_t leaf_mapping[5], Node* nodes[8], Edge* edges[7], const DistanceMatrix &dist, double &best_fit, bool apply)
{
	static const double data[] = {
		1, 1, 0, 0, 0, 0, 0,
		1, 0, 1, 0, 0, 1, 0,
		1, 0, 0, 1, 0, 1, 0,
		1, 0, 0, 0, 1, 0, 0,
		0, 1, 1, 0, 0, 1, 0,
		0, 1, 0, 1, 0, 1, 0,
		0, 1, 0, 0, 1, 0, 1,
		0, 0, 1, 1, 0, 0, 0,
		0, 0, 1, 0, 1, 1, 1,
		0, 0, 0, 1, 1, 1, 1};

	Eigen::Matrix<double,10,1> dists;
	Eigen::Matrix<double,10,1> weights;
	Eigen::Matrix<double,7,1> new_dists;
	Eigen::Matrix<double,10,7> A;

	dists[0] = dist.distances(leaf_mapping[0],leaf_mapping[1]);
	dists[1] = dist.distances(leaf_mapping[0],leaf_mapping[2]);
	dists[2] = dist.distances(leaf_mapping[0],leaf_mapping[3]);
	dists[3] = dist.distances(leaf_mapping[0],leaf_mapping[4]);
	dists[4] = dist.distances(leaf_mapping[1],leaf_mapping[2]);
	dists[5] = dist.distances(leaf_mapping[1],leaf_mapping[3]);
	dists[6] = dist.distances(leaf_mapping[1],leaf_mapping[4]);
	dists[7] = dist.distances(leaf_mapping[2],leaf_mapping[3]);
	dists[8] = dist.distances(leaf_mapping[2],leaf_mapping[4]);
	dists[9] = dist.distances(leaf_mapping[3],leaf_mapping[4]);

	weights[0] = dist.variances(leaf_mapping[0],leaf_mapping[1]);
	weights[1] = dist.variances(leaf_mapping[0],leaf_mapping[2]);
	weights[2] = dist.variances(leaf_mapping[0],leaf_mapping[3]);
	weights[3] = dist.variances(leaf_mapping[0],leaf_mapping[4]);
	weights[4] = dist.variances(leaf_mapping[1],leaf_mapping[2]);
	weights[5] = dist.variances(leaf_mapping[1],leaf_mapping[3]);
	weights[6] = dist.variances(leaf_mapping[1],leaf_mapping[4]);
	weights[7] = dist.variances(leaf_mapping[2],leaf_mapping[3]);
	weights[8] = dist.variances(leaf_mapping[2],leaf_mapping[4]);
	weights[9] = dist.variances(leaf_mapping[3],leaf_mapping[4]);

	A.noalias() = weights.asDiagonal() *  Eigen::Map<const Eigen::Matrix<double,10,7,Eigen::RowMajor> >(data);

	new_dists = NNLS(A,dists);

	double fit = (A*new_dists - dists).squaredNorm();

	if(fit < best_fit && apply) {
		best_fit = fit;
		for(index_t i=0; i<5; ++i) {
			edges[leaf_mapping[i]]->length = new_dists(i);
			edges[leaf_mapping[i]]->nodes[0] = nodes[leaf_mapping[i]];
		}
		edges[5]->length = new_dists(5);
		edges[5]->nodes[0] = nodes[5];
		edges[5]->nodes[1] = nodes[6];

		edges[6]->length = new_dists(6);
		edges[6]->nodes[0] = nodes[5];
		edges[6]->nodes[1] = nodes[7];

		nodes[5]->edges[0] = edges[5];
		nodes[5]->edges[1] = edges[6];
		nodes[5]->edges[2] = edges[leaf_mapping[1]];

		nodes[6]->edges[0] = edges[5];
		nodes[6]->edges[1] = edges[leaf_mapping[2]];
		nodes[6]->edges[2] = edges[leaf_mapping[3]];

		nodes[7]->edges[0] = edges[6];
		nodes[7]->edges[1] = edges[leaf_mapping[0]];
		nodes[7]->edges[2] = edges[leaf_mapping[4]];

		edges[leaf_mapping[0]]->nodes[1] = nodes[7];
		edges[leaf_mapping[1]]->nodes[1] = nodes[5];
		edges[leaf_mapping[2]]->nodes[1] = nodes[6];
		edges[leaf_mapping[3]]->nodes[1] = nodes[6];
		edges[leaf_mapping[4]]->nodes[1] = nodes[7];

		nodes[5]->todo = true;
		nodes[6]->todo = true;
		nodes[7]->todo = true;
	}

	return fit;
}

// Node 4 in center position
static double Opt5v2(const index_t leaf_mapping[5], Node* nodes[8], Edge* edges[7], const DistanceMatrix &dist, double &best_fit, bool apply) {
	static const double data[] = {
		1, 1, 0, 0, 0, 0, 0,
		1, 0, 1, 0, 0, 1, 1,
		1, 0, 0, 1, 0, 1, 1,
		1, 0, 0, 0, 1, 1, 0,
		0, 1, 1, 0, 0, 1, 1,
		0, 1, 0, 1, 0, 1, 1,
		0, 1, 0, 0, 1, 1, 0,
		0, 0, 1, 1, 0, 0, 0,
		0, 0, 1, 0, 1, 0, 1,
		0, 0, 0, 1, 1, 0, 1};

	Eigen::Matrix<double,10,1> dists;
	Eigen::Matrix<double,10,1> weights;
	Eigen::Matrix<double,7,1> new_dists;
	Eigen::Matrix<double,10,7> A;

	dists[0] = dist.distances(leaf_mapping[0],leaf_mapping[1]);
	dists[1] = dist.distances(leaf_mapping[0],leaf_mapping[2]);
	dists[2] = dist.distances(leaf_mapping[0],leaf_mapping[3]);
	dists[3] = dist.distances(leaf_mapping[0],leaf_mapping[4]);
	dists[4] = dist.distances(leaf_mapping[1],leaf_mapping[2]);
	dists[5] = dist.distances(leaf_mapping[1],leaf_mapping[3]);
	dists[6] = dist.distances(leaf_mapping[1],leaf_mapping[4]);
	dists[7] = dist.distances(leaf_mapping[2],leaf_mapping[3]);
	dists[8] = dist.distances(leaf_mapping[2],leaf_mapping[4]);
	dists[9] = dist.distances(leaf_mapping[3],leaf_mapping[4]);

	weights[0] = dist.variances(leaf_mapping[0],leaf_mapping[1]);
	weights[1] = dist.variances(leaf_mapping[0],leaf_mapping[2]);
	weights[2] = dist.variances(leaf_mapping[0],leaf_mapping[3]);
	weights[3] = dist.variances(leaf_mapping[0],leaf_mapping[4]);
	weights[4] = dist.variances(leaf_mapping[1],leaf_mapping[2]);
	weights[5] = dist.variances(leaf_mapping[1],leaf_mapping[3]);
	weights[6] = dist.variances(leaf_mapping[1],leaf_mapping[4]);
	weights[7] = dist.variances(leaf_mapping[2],leaf_mapping[3]);
	weights[8] = dist.variances(leaf_mapping[2],leaf_mapping[4]);
	weights[9] = dist.variances(leaf_mapping[3],leaf_mapping[4]);

	A.noalias() = weights.asDiagonal() *  Eigen::Map<const Eigen::Matrix<double,10,7,Eigen::RowMajor> >(data);

	new_dists = NNLS(A,dists);

	double fit = (A*new_dists - dists).squaredNorm();

	if(fit < best_fit && apply) {
		best_fit = fit;
		for(index_t i=0; i<5; ++i) {
			edges[leaf_mapping[i]]->length = new_dists(i);
			edges[leaf_mapping[i]]->nodes[0] = nodes[leaf_mapping[i]];
		}
		edges[5]->length = new_dists(5);
		edges[5]->nodes[0] = nodes[5];
		edges[5]->nodes[1] = nodes[7];

		edges[6]->length = new_dists(6);
		edges[6]->nodes[0] = nodes[6];
		edges[6]->nodes[1] = nodes[7];

		nodes[5]->edges[0] = edges[5];
		nodes[5]->edges[1] = edges[leaf_mapping[0]];
		nodes[5]->edges[2] = edges[leaf_mapping[1]];

		nodes[6]->edges[0] = edges[6];
		nodes[6]->edges[1] = edges[leaf_mapping[2]];
		nodes[6]->edges[2] = edges[leaf_mapping[3]];

		nodes[7]->edges[0] = edges[5];
		nodes[7]->edges[1] = edges[6];
		nodes[7]->edges[2] = edges[leaf_mapping[4]];

		edges[leaf_mapping[0]]->nodes[1] = nodes[5];
		edges[leaf_mapping[1]]->nodes[1] = nodes[5];
		edges[leaf_mapping[2]]->nodes[1] = nodes[6];
		edges[leaf_mapping[3]]->nodes[1] = nodes[6];
		edges[leaf_mapping[4]]->nodes[1] = nodes[7];

		nodes[5]->todo = true;
		nodes[6]->todo = true;
		nodes[7]->todo = false;
	}

	return fit;
}

// pass center node and edge to node 4
static bool OptimizeQuintet(Node *n, Edge *e, Graph *g, const DistanceMatrix &all_dist, bool apply)
{
	double best_fit = INFINITY;
	Node *nodes[8];
	Edge *edges[7];

	edges[4] = e;
	nodes[7] = n;
	nodes[4] = &edges[4]->other(*n);

	if(nodes[7]->isLeaf()) return false;

	edges[5] = &(*nodes[7])[0] == e ? &(*nodes[7])[1] : &(*nodes[7])[0];
	edges[6] = &(*nodes[7])[2] == e ? &(*nodes[7])[1] : &(*nodes[7])[2];

	nodes[5] = &edges[5]->other(*nodes[7]);
	nodes[6] = &edges[6]->other(*nodes[7]);

	if(nodes[5]->isLeaf() || nodes[6]->isLeaf()) return false;

	edges[0] = &(*nodes[5])[0] == edges[5] ? &(*nodes[5])[1] : &(*nodes[5])[0];
	edges[1] = &(*nodes[5])[2] == edges[5] ? &(*nodes[5])[1] : &(*nodes[5])[2];
	edges[2] = &(*nodes[6])[0] == edges[6] ? &(*nodes[6])[1] : &(*nodes[6])[0];
	edges[3] = &(*nodes[6])[2] == edges[6] ? &(*nodes[6])[1] : &(*nodes[6])[2];

	nodes[0] = &edges[0]->other(*nodes[5]);
	nodes[1] = &edges[1]->other(*nodes[5]);
	nodes[2] = &edges[2]->other(*nodes[6]);
	nodes[3] = &edges[3]->other(*nodes[6]);


	std::map<const Node*,distance_t> leaf_distances[5];
	for(index_t i=0; i<5; ++i) {
		leaf_distances[i] = g->subtreeDist(nodes[i],edges[i]);
	}
	DistanceMatrix dist(5);
	for(index_t i=0; i<5; ++i) {
		for(index_t j=i+1; j<5; ++j) {
			dist.distances(i,j) = 0;
			dist.variances(i,j) = 0;
			for(std::map<const Node*,distance_t>::iterator k = leaf_distances[i].begin(); k != leaf_distances[i].end(); ++k) {
				for(std::map<const Node*,distance_t>::iterator l = leaf_distances[j].begin(); l != leaf_distances[j].end(); ++l) {
					index_t kk = k->first->leaf;
					index_t ll = l->first->leaf;
					dist.distances(i,j) += all_dist.variances(kk,ll) * (all_dist.distances(kk,ll) - k->second - l->second);
					dist.variances(i,j) += all_dist.variances(kk,ll);
				}
			}
			dist.distances(j,i) = dist.distances(i,j);
			dist.variances(j,i) = dist.variances(i,j);
		}
	}
	dist.variances = dist.variances.array().sqrt();
	dist.distances.array() /= dist.variances.array();

	static const index_t mapping11[5] = {0,1,2,3,4};
	double f1 = Opt5v2(mapping11,nodes,edges,dist,best_fit,apply);

	static const index_t mapping12[5] = {0,1,2,3,4};
	Opt5v1(mapping12,nodes,edges,dist,best_fit,apply);

	static const index_t mapping13[5] = {1,0,2,3,4};
	Opt5v1(mapping13,nodes,edges,dist,best_fit,apply);

	static const index_t mapping14[5] = {2,3,0,1,4};
	Opt5v1(mapping14,nodes,edges,dist,best_fit,apply);

	static const index_t mapping15[5] = {3,2,0,1,4};
	Opt5v1(mapping15,nodes,edges,dist,best_fit,apply);


	static const index_t mapping21[5] = {0,2,1,3,4};
	Opt5v2(mapping21,nodes,edges,dist,best_fit,apply);

	static const index_t mapping22[5] = {0,2,1,3,4};
	Opt5v1(mapping22,nodes,edges,dist,best_fit,apply);

	static const index_t mapping23[5] = {2,0,1,3,4};
	Opt5v1(mapping23,nodes,edges,dist,best_fit,apply);

	static const index_t mapping24[5] = {1,3,0,2,4};
	Opt5v1(mapping24,nodes,edges,dist,best_fit,apply);

	static const index_t mapping25[5] = {3,1,0,2,4};
	Opt5v1(mapping25,nodes,edges,dist,best_fit,apply);


	static const index_t mapping31[5] = {0,3,1,2,4};
	Opt5v2(mapping31,nodes,edges,dist,best_fit,apply);

	static const index_t mapping32[5] = {0,3,1,2,4};
	Opt5v1(mapping32,nodes,edges,dist,best_fit,apply);

	static const index_t mapping33[5] = {3,0,1,2,4};
	Opt5v1(mapping33,nodes,edges,dist,best_fit,apply);

	static const index_t mapping34[5] = {1,2,0,3,4};
	Opt5v1(mapping34,nodes,edges,dist,best_fit,apply);

	static const index_t mapping35[5] = {2,1,0,3,4};
	Opt5v1(mapping35,nodes,edges,dist,best_fit,apply);

	return best_fit < f1;
}

static void OptimizeQuintets(Graph *g, const DistanceMatrix &all_dist, bool apply)
{
	for(index_t i=0; i < g->n_nodes; ++i) {
		g->nodes[i].todo = true;
	}

	for(index_t k=0; k < MAX_ITERS5; ++k) {
		bool any = false;
		for(index_t i=0; i < g->n_nodes; ++i) {
			Node *n = &g->nodes[i];
			if(n->todo == false) continue;
			n->todo = false;
			if(n->isLeaf()) continue;

			for(index_t j=0; j<3; ++j) {
				Edge *e = n->edges[j];
				bool changed = OptimizeQuintet(n,e,g,all_dist,apply);
				g->verify();

				if(changed) {
					any = true;
					break;
				}
			}
		}
		if(!any) break;
	}
}


PhyTree* refineTree(PhyTree *tree, const std::vector<std::string> &leaf_order, const DistanceMatrix &dist) {
	Graph g(*tree,leaf_order);
	g.verify();

	DistanceMatrix weights = dist;
	weights.variances = weights.variances.array().inverse();

	double fit1 = computeFit(g,weights);
#ifdef DEBUG
	std::cerr << "Initial fit: " << fit1 << std::endl;
#endif

	OptimizeQuartets(&g,weights,true);
	double fit2 = computeFit(g,weights);
#ifdef DEBUG
	std::cerr << "After OptimizeQuartets: " << fit2 << std::endl;
#endif

	index_t i=0;
	do {
		fit1 = fit2;

		if(cmdlineopts.wlsrefine_flag > 1) {
			OptimizeQuintets(&g,weights,true);
			fit2 = computeFit(g,weights);
#ifdef DEBUG
			std::cerr << "After OptimizeQuintets: " << fit2 << std::endl;
#endif
		}

		OptimizeQuartets(&g,weights,true);
		fit2 = computeFit(g,weights);
#ifdef DEBUG
		std::cerr << "After OptimizeQuartets: " << fit2 << std::endl;
#endif

		++i;
	} while(fit2 < fit1 && i < MAX_ITERS);

#ifdef DEBUG
	std::cerr << "Final fit: " << fit2 << std::endl;
#endif

	// compute support only
	OptimizeQuartets(&g,weights,false);

	delete tree;
	tree = g.toTree();
	return tree;
}
}
