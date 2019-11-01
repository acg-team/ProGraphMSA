#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <iterator>
#include <cassert>
#include <iostream>
#include <cmath>
#include <utility>
#include <map>

#include "main.h"
#include "Model.h"
#include "Repeat.h"

#include <Eigen/SparseCore>

#define MAX_EDGE_COST ((dp_score_t)10000)

template <class ALPHABET>
class Graph {
private:
	typedef typename Model<ALPHABET>::Profile Profile;
	typedef typename Model<ALPHABET>::Freqs Freqs;
	typedef typename Model<ALPHABET>::Subst Subst;

protected:
	typedef Eigen::SparseMatrix<dp_score_t,Eigen::RowMajor> AdjacencyMatrix;
	typedef Eigen::SparseMatrix<index_t,Eigen::RowMajor> TRMatrix;

	Profile sites;
	AdjacencyMatrix edges;
	TRMatrix repeats;

	void fillInitialEdges() {
		Eigen::VectorXi cols = Eigen::VectorXi::Zero(this->size());
		for(index_t i=0; i<this->size(); ++i) {
			cols[i] = this->edges.row(i).nonZeros() + (i>0);
		}
		this->edges.reserve(cols);

		for(index_t i=0; i<this->size()-1; ++i) {
			this->edges.coeffRef(i+1,i) = -MAX_EDGE_COST;
		}
		this->edges.makeCompressed();
	}

	void getRepeatEdges(std::map<std::pair<index_t,index_t>,index_t> &repeat_map, const std::vector<int> &tr_homology, index_t offset = 0) {
		// XXX include the positions before first/after last i.e. deletion of first/last position?
		for(index_t from=0; from != tr_homology.size(); ++from) {
			if(tr_homology[from] < 0) continue;

			index_t n_units = 0;
			bool take_next = false;

			for(index_t to=from+1; to != tr_homology.size(); ++to) {
				if(tr_homology[to] < 0) continue;

				if(tr_homology[to] <= tr_homology[to-1]) {
					n_units += 1;
				}

				if(take_next) {
					std::pair<index_t,index_t> i(offset+to,offset+from);
					typename std::map<std::pair<index_t,index_t>,index_t>::iterator it = repeat_map.find(i);
					if(it != repeat_map.end()) {
						it->second = std::min(it->second,n_units);
					} else {
						repeat_map[i] = n_units;
					}
					take_next = false;
				}

				if(tr_homology[to] == tr_homology[from]) {
					take_next = true;
				}
			}
		}
	}

	void setEdgesFromMap(const std::map<std::pair<index_t,index_t>,dp_score_t> &edge_map) {
		std::vector<Eigen::Triplet<dp_score_t> > edges2;
		edges2.reserve(edge_map.size());
		for(std::map<std::pair<index_t,index_t>,dp_score_t>::const_iterator it = edge_map.begin(); it != edge_map.end(); ++it) {
			dp_score_t cost = std::min(it->second,MAX_EDGE_COST) - MAX_EDGE_COST;
			edges2.push_back(Eigen::Triplet<dp_score_t>(it->first.first,it->first.second,cost));
		}
		this->edges.setFromTriplets(edges2.begin(),edges2.end());
		this->edges.makeCompressed();
	}

	void setRepeatsFromMap(const std::map<std::pair<index_t,index_t>,index_t> &repeat_map) {
		std::vector<Eigen::Triplet<index_t> > repeats2;
		repeats2.reserve(repeat_map.size());
		for(std::map<std::pair<index_t,index_t>,index_t>::const_iterator it = repeat_map.begin(); it != repeat_map.end(); ++it) {
			repeats2.push_back(Eigen::Triplet<index_t>(it->first.first,it->first.second,it->second));
		}
		this->repeats.setFromTriplets(repeats2.begin(),repeats2.end());
		this->repeats.makeCompressed();
	}


public:
	Graph(const std::vector<Freqs, Eigen::aligned_allocator<Freqs> > &nodes) {
		index_t dim = nodes.size();

		assert(dim >= 2 && nodes[0].isZero(0) && nodes[dim-1].isZero(0));

		this->edges = AdjacencyMatrix(dim,dim);
		this->edges.setZero();
		this->sites = Profile::Zero(ALPHABET::DIM,dim);
		this->repeats = TRMatrix(dim,dim);
		this->repeats.setZero();

		for(index_t i=1; i<dim-1; ++i) {
			this->sites.col(i) = nodes[i];
		}

		fillInitialEdges();
	}

	Graph(const std::vector<Freqs, Eigen::aligned_allocator<Freqs> > &nodes, const std::map<std::pair<index_t,index_t>,dp_score_t> &edges, const std::map<std::pair<index_t,index_t>,index_t> &repeats) {
		index_t dim = nodes.size();

		assert(dim >= 2 && nodes[0].isZero(0) && nodes[dim-1].isZero(0));

		this->edges = AdjacencyMatrix(dim,dim);
		this->edges.setZero();
		this->sites = Profile::Zero(ALPHABET::DIM,dim);
		this->repeats = TRMatrix(dim,dim);
		this->repeats.setZero();

		for(index_t i=1; i<dim-1; ++i) {
			this->sites.col(i) = nodes[i];
		}

		this->setEdgesFromMap(edges);
		this->setRepeatsFromMap(repeats);
	}

	Graph() {
		index_t dim = 2;
		this->edges = AdjacencyMatrix(dim,dim);
		this->edges.setZero();
		this->sites = Profile::Zero(ALPHABET::DIM,2);
		this->repeats = TRMatrix(dim,dim);
		this->repeats.setZero();

		fillInitialEdges();
	}

	Graph(const Graph<ALPHABET> &other) {
		this->edges = other.edges;
		this->sites = other.sites;
		this->repeats = other.repeats;
	}

	Graph<ALPHABET>& operator=(const Graph<ALPHABET> &other) {
		this->edges = other.edges;
		this->sites = other.sites;
		this->repeats = other.repeats;

		return *this;
	}

	virtual ~Graph() {}

	index_t size() const {
		return this->sites.cols();
	}

	index_t size_edges() const {
		return this->edges.nonZeros();
	}

	index_t size_repeats() const {
		return this->repeats.nonZeros();
	}

	class PredIterator : public std::iterator<std::input_iterator_tag, index_t> {
	private:
		const Graph<ALPHABET> *graph;
		AdjacencyMatrix::InnerIterator it;
		TRMatrix::InnerIterator it2;
		dp_score_t repeatInit;
		dp_score_t repeatExt;
	public:
		PredIterator(const Graph<ALPHABET> &graph, index_t row, dp_score_t repeat_init, dp_score_t repeat_ext)
			: graph(&graph), it(AdjacencyMatrix::InnerIterator(graph.edges,row)), it2(TRMatrix::InnerIterator(graph.repeats,row)), repeatInit(repeat_init), repeatExt(repeat_ext)
		{}

		PredIterator(const PredIterator &other)
			: graph(other.graph), it(other.it), it2(other.it2), repeatInit(other.repeatInit), repeatExt(other.repeatExt)
		{}

		 PredIterator& operator++() {if(it) ++it; else ++it2; return *this;}

		 bool operator==(const PredIterator& rhs) const {
			 return graph==rhs.graph && it.col() == rhs.it.col() && it.row() == rhs.it.row()
					 && it2.col() == rhs.it2.col() && it2.row() == rhs.it2.row()
					 && repeatInit == rhs.repeatInit && repeatExt == rhs.repeatExt;
		 }

		 bool operator!=(const PredIterator& rhs) const {
			 return graph!=rhs.graph || it.col() != rhs.it.col() || it.row() != rhs.it.row()
					 || it2.col() != rhs.it2.col() || it2.row() != rhs.it2.row()
					 || repeatInit != rhs.repeatInit || repeatExt != rhs.repeatExt;
		 }

		 bool isRepeat() const {
			 return !it;
		 }

		 index_t repeatUnits() const {
			 assert(!it && it2);
			 return it2.value();
		 }

		 operator bool() const {
			 return (bool)it || (bool)it2;
		 }

		 dp_score_t value() const {
			 if(it) {
				dp_score_t c = it.value();
				if(c == 0) {
					return (dp_score_t)INFINITY;
				} else {
					return c + MAX_EDGE_COST;
				}
			 } else {
				index_t c = it2.value();
				if(c == 0) {
					return (dp_score_t)INFINITY;
				} else {
					return this->repeatInit + this->repeatExt * (c-1);
				}
			 }
		 }

		 index_t operator*() const {
			 if(it) {
				 return it.col();
			 } else {
				 return it2.col();
			 }
		 }
	};
	friend class PredIterator;

	void addNodes(index_t before, const std::vector<Freqs, Eigen::aligned_allocator<Freqs> > &nodes) {
		assert(before >= 1 && before <= this->size() - 1);

		index_t oldDim = this->size();
		index_t newDim = oldDim + nodes.size();

		std::vector<Eigen::Triplet<dp_score_t> > edges;
		edges.reserve(this->edges.nonZeros() + nodes.size() - 1);

		for (index_t i = 0; i < (index_t)this->edges.rows(); ++i) {
			for (AdjacencyMatrix::InnerIterator it(this->edges, i); it; ++it) {
				index_t y = i + (i >= before) * nodes.size();
				index_t x = it.index() + (it.index() >= before) * nodes.size();
				edges.push_back(Eigen::Triplet<dp_score_t>(y,x,it.value()));
			}
		}

		for(index_t i = before; i < before + nodes.size() - 1; ++i) {
			edges.push_back(Eigen::Triplet<dp_score_t>(i+1,i,-MAX_EDGE_COST));
		}

		this->edges.setZero();
		this->edges.resize(newDim,newDim);
		this->edges.setFromTriplets(edges.begin(),edges.end());
		this->edges.makeCompressed();
		edges.clear();


		std::vector<Eigen::Triplet<index_t> > tredges;
		tredges.reserve(this->repeats.nonZeros());

		for (index_t i = 0; i < (index_t)this->repeats.rows(); ++i) {
			for (TRMatrix::InnerIterator it(this->repeats, i); it; ++it) {
				index_t y = i + (i >= before) * nodes.size();
				index_t x = it.index() + (it.index() >= before) * nodes.size();
				tredges.push_back(Eigen::Triplet<index_t>(y,x,it.value()));
			}
		}

		this->repeats.setZero();
		this->repeats.resize(newDim,newDim);
		this->repeats.setFromTriplets(tredges.begin(),tredges.end());
		this->repeats.makeCompressed();
		tredges.clear();


		Profile newNodes((int)ALPHABET::DIM,nodes.size());

		for(index_t i = 0; i < nodes.size(); ++i) {
			newNodes.col(i) = nodes[i];
		}

		Profile newSites((int)ALPHABET::DIM,newDim);

		newSites << this->sites.block(0,0,this->sites.rows(),before),
				newNodes,
				this->sites.block(0,before,this->sites.rows(),this->sites.cols()-before);

		this->sites = newSites;
	}

	void addNodesNoEdges(index_t before, const std::vector<Freqs, Eigen::aligned_allocator<Freqs> > &nodes) {
		assert(before >= 1 && before <= this->size() - 1);

		index_t oldDim = this->size();
		index_t newDim = oldDim + nodes.size();

		std::vector<Eigen::Triplet<dp_score_t> > edges;
		edges.reserve(this->edges.nonZeros());

		for (index_t i = 0; i < (index_t)this->edges.rows(); ++i) {
			for (AdjacencyMatrix::InnerIterator it(this->edges, i); it; ++it) {
				index_t y = i + (i >= before) * nodes.size();
				index_t x = it.index() + ((index_t)it.index() >= before) * nodes.size();
				edges.push_back(Eigen::Triplet<dp_score_t>(y,x,it.value()));
			}
		}

		this->edges.setZero();
		this->edges.resize(newDim,newDim);
		this->edges.setFromTriplets(edges.begin(),edges.end());
		this->edges.makeCompressed();
		edges.clear();


		std::vector<Eigen::Triplet<index_t> > tredges;
		tredges.reserve(this->repeats.nonZeros());

		for (index_t i = 0; i < (index_t)this->repeats.rows(); ++i) {
			for (TRMatrix::InnerIterator it(this->repeats, i); it; ++it) {
				index_t y = i + (i >= before) * nodes.size();
				index_t x = it.index() + ((index_t)it.index() >= before) * nodes.size();
				tredges.push_back(Eigen::Triplet<index_t>(y,x,it.value()));
			}
		}

		this->repeats.setZero();
		this->repeats.resize(newDim,newDim);
		this->repeats.setFromTriplets(tredges.begin(),tredges.end());
		this->repeats.makeCompressed();
		tredges.clear();


		Profile newNodes((int)ALPHABET::DIM,nodes.size());

		for(index_t i = 0; i < nodes.size(); ++i) {
			newNodes.col(i) = nodes[i];
		}

		Profile newSites((int)ALPHABET::DIM,newDim);

		newSites << this->sites.block(0,0,this->sites.rows(),before),
				newNodes,
				this->sites.block(0,before,this->sites.rows(),this->sites.cols()-before);

		this->sites = newSites;
	}

	void reset() {
		this->sites.block(0,0,this->sites.rows(),this->sites.cols()-1).setOnes();
		this->sites.col(0).setZero();
		this->sites.col(this->sites.cols()-1).setZero();
	}

	void rmNodes(index_t first, index_t count = 1) {
		assert(first >= 1 && first+count <= this->size());

		index_t oldDim = this->size();
		index_t newDim = oldDim - count;

		std::vector<Eigen::Triplet<dp_score_t> > edges;
		edges.reserve(this->edges.nonZeros());

		for (index_t i = 0; i < (index_t)this->edges.rows(); ++i) {
			for (AdjacencyMatrix::InnerIterator it(this->edges, i); it; ++it) {
				if((i < first || i >= first+count) && ((index_t)it.index() < first || (index_t)it.index() >= first+count)) {
					index_t y = i - (i >= first) * count;
					index_t x = it.index() - ((index_t)it.index() >= first) * count;
					edges.push_back(Eigen::Triplet<dp_score_t>(y,x,it.value()));
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

		for (index_t i = 0; i < (index_t)this->repeats.rows(); ++i) {
			for (TRMatrix::InnerIterator it(this->repeats, i); it; ++it) {
				if((i < first || i >= first+count) && ((index_t)it.index() < first || (index_t)it.index() >= first+count)) {
					index_t y = i - (i >= first) * count;
					index_t x = it.index() - ((index_t)it.index() >= first) * count;
					tredges.push_back(Eigen::Triplet<index_t>(y,x,it.value()));
				}
			}
		}

		this->repeats.setZero();
		this->repeats.resize(newDim,newDim);
		this->repeats.setFromTriplets(tredges.begin(),tredges.end());
		this->repeats.makeCompressed();
		tredges.clear();


		Profile newSites((int)ALPHABET::DIM,newDim);

		newSites << this->sites.block(0,0,this->sites.rows(),first),
				this->sites.block(0,first+count,this->sites.rows(),this->sites.cols()-first-count);

		this->sites = newSites;
	}

	PredIterator getPreds(index_t node, dp_score_t repeat_init, dp_score_t repeat_ext) const {
		return PredIterator(*this, node, repeat_init, repeat_ext);
	}

	Freqs operator[](index_t i) const {
		return this->sites.col(i);
	}

	const Profile& getSites() const {
		return this->sites;
	}


#if 0
	void addRepeats(const std::vector<repeat_t> &repeats) {
		index_t dim = this->size();
		this->repeats = TRMatrix(dim,dim);
		this->repeats.setZero();
		std::vector<Eigen::Triplet<index_t> > tredges;

		for(std::vector<repeat_t>::const_iterator rep=repeats.begin(); rep != repeats.end(); ++rep) {
			this->get_repeat_edges(tredges,rep->tr_hom,rep->start);
		}

		this->setRepeatsFromMap(tredges);
	}
#endif


	void addRepeats(const std::vector<std::vector<int> > &tr_homologies) {
		index_t dim = this->size();
		this->repeats = TRMatrix(dim,dim);
		this->repeats.setZero();
		std::map<std::pair<index_t,index_t>,index_t> tredges;

		for(std::vector<std::vector<int> >::const_iterator tr_homology = tr_homologies.begin(); tr_homology != tr_homologies.end(); ++tr_homology) {
			this->getRepeatEdges(tredges,*tr_homology);
		}

		this->setRepeatsFromMap(tredges);
	}


	void printDot(const std::string &name="DAG") const {
		std::cout << "digraph \"" << name << "\" {" << std::endl;
		std::cout << "\trankdir=LR;" << std::endl;
		std::cout << "\tsplines=true;" << std::endl;

		std::cout << std::endl;

		std::cout << "\tn0 [label=\"START\"];" << std::endl;
		for(index_t i=1; i < this->size()-1; ++i) {
			typename Profile::Index x;
			this->sites.col(i).maxCoeff(&x);
			std::cout << "\tn" << i << " [label=\"" << i << ": " << ALPHABET((int)x).asString() << "\"];" << std::endl;
		}
		std::cout << "\tn" << this->size()-1 << " [label=\"END\"];" << std::endl;

		std::cout << std::endl;

		for (index_t i = 0; i < (index_t)this->edges.rows(); ++i) {
			for (AdjacencyMatrix::InnerIterator it(this->edges, i); it; ++it) {
				std::cout << "\tn" << it.index() << " -> n" << i
						<< " [label=\"" << (it.value()+MAX_EDGE_COST) << "\"];" << std::endl;
			}
			for (TRMatrix::InnerIterator it(this->repeats, i); it; ++it) {
				std::cout << "\tn" << it.index() << " -> n" << i
						<< " [label=\"" << it.value() << "\",color=\"blue\"];" << std::endl;
			}
		}

		std::cout << "}" << std::endl;
	}
};

#endif /* GRAPH_H_ */
