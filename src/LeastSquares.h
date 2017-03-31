#ifndef LEASTSQUARES_H_
#define LEASTSQUARES_H_

#include "main.h"
#include "PhyTree.h"
#include "DistanceFactory.h"
#include <string>
#include <vector>

namespace LeastSquares {
PhyTree* refineTree(PhyTree *tree, const std::vector<std::string> &leaf_order, const DistanceMatrix &dist);
}

#endif
