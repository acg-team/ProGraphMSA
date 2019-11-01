#ifndef PHYTREE_H_
#define PHYTREE_H_

#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <sstream>
#include <cmath>
#include "main.h"

class PhyTree {
private:
	std::vector<PhyTree*> children;
	PhyTree *parent;
	double branch_length;
	double branch_support;
	std::string name;

	void print_prefix(std::string prefix) const {
		std::cout << prefix << "(" << this->branch_length << ") " << this->name << std::endl;
		for(std::vector<PhyTree*>::const_iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->print_prefix(prefix+"  ");
		}
	}

	void fixDistancesR() {
		if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
			if(std::isnan(this->branch_length)) this->branch_length = cmdlineopts.max_dist;
			this->branch_length = std::min(std::max(cmdlineopts.min_dist,this->branch_length),cmdlineopts.max_dist);
		} else {
			if(std::isnan(this->branch_length)) this->branch_length = cmdlineopts.max_pdist;
			this->branch_length = std::min(std::max(cmdlineopts.min_pdist,this->branch_length),cmdlineopts.max_pdist);
		}
		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->fixDistancesR();
		}
	}

	std::string formatNewickR() const {
		if(this->isLeaf()) {
			return this->getName();
		} else {
			std::stringstream newick;
			newick << "(";
			std::vector<PhyTree*>::const_iterator i=this->children.begin();
			newick << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
			for(++i; i < this->children.end(); ++i) {
				newick << "," << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
			}
			newick << ")";
			return newick.str();
		}
	}

public:
	PhyTree(std::string name="") {
		this->parent = NULL;
		this->branch_length = 0;
		this->branch_support = 1;
		this->name = name;
	}

	~PhyTree() {
		assert(this->parent == NULL);
		for(std::vector<PhyTree*>::reverse_iterator i=this->children.rbegin(); i < this->children.rend(); ++i) {
			PhyTree *child = *i;
			child->parent = NULL;
			child->branch_length = 0;

			delete child;
		}
	}

	PhyTree* copy() {
		PhyTree* out = new PhyTree();
		out->branch_length = this->branch_length;
		out->branch_support = this->branch_support;
		out->name = this->name;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			out->addChild((*i)->copy(),(*i)->branch_length,(*i)->branch_support);
		}
		return out;
	}

	void addChild(PhyTree *child, double branch_length = 0, double branch_support = 1) {
		assert(child != this);
		assert(child->parent == NULL);

		this->children.push_back(child);
		child->parent = this;
		child->branch_length = branch_length;
		child->branch_support = branch_support;
	}

	index_t indexOf() {
		PhyTree *parent = this->parent;
		assert(parent != NULL);

		for(index_t i=0; i < parent->children.size(); ++i) {
			if(parent->children[i] == this) {
				return i;
			}
		}

		assert(false);
		return -1;
	}

	void pluck() {
		assert(this->parent != NULL);

		index_t index = this->indexOf();
		std::vector<PhyTree*>::iterator iter = this->parent->children.begin()+index;

		this->parent->children.erase(iter);
		this->parent = NULL;

		this->branch_length = 0;
		this->branch_support = 1;
	}

	PhyTree* pluckChild(index_t index) {
		std::vector<PhyTree*>::iterator iter = this->children.begin()+index;
		PhyTree *child = *iter;

		this->children.erase(iter);
		child->parent = NULL;
		child->branch_length = 0;
		this->branch_support = 1;

		return child;
	}

	void deleteChild(index_t index) {
		PhyTree *child = this->pluckChild(index);
		delete child;
	}

	void fixDistances() {
		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->fixDistancesR();
		}
	}

	index_t countLeaves() {
		if(this->isLeaf()) {
			return 1;
		} else {
			index_t leaves = 0;
			for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
				leaves += (*i)->countLeaves();
			}
			return leaves;
		}
	}

	double computeLength() {
		if(this->isLeaf()) {
			return 0;
		} else {
			double length = 0;
			for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
				length += (*i)->branch_length + (*i)->computeLength();
			}
			return length;
		}
	}

	std::string getName() const {
		return this->name;
	}

	PhyTree *getParent() {
		return this->parent;
	}

	const PhyTree *getParent() const {
		return this->parent;
	}

	double getBranchLength() const {
		return this->branch_length;
	}

	double getBranchSupport() const {
		return this->branch_support;
	}

	PhyTree& operator[](int i) {
		return *this->children[i];
	}

	const PhyTree& operator[](int i) const {
		return *this->children[i];
	}

	std::vector<PhyTree*>::iterator firstChild() {
		return this->children.begin();
	}

	std::vector<PhyTree*>::iterator lastChild() {
		return this->children.end();
	}

	std::vector<PhyTree*>::const_iterator firstChild() const {
		return this->children.begin();
	}

	std::vector<PhyTree*>::const_iterator lastChild() const {
		return this->children.end();
	}

	index_t n_children() const {
		return this->children.size();
	}

	bool isLeaf() const {
		return this->children.empty();
	}

	void print() const {
		this->print_prefix("");
	}

	std::string formatNewick() const {
		return this->formatNewickR() + ";";
	}
};

PhyTree* midpointRoot(PhyTree *root);

std::vector<std::string> get_tree_order_ancestral(const PhyTree *tree);
std::vector<std::string> get_tree_order(const PhyTree *tree);

#endif /* PHYTREE_H_ */
