#ifndef STOCKHOLM_H_
#define STOCKHOLM_H_

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "main.h"
#include "PhyTree.h"

void write_stockholm(const std::map<std::string,std::string > alignment, std::ostream &out);
void write_stockholm(const std::map<std::string,std::string > alignment, const std::vector<std::string> &order, const PhyTree &tree, std::ostream &out);
void write_stockholm(const std::map<std::string,std::string > alignment, const std::vector<std::string> &order, const PhyTree &tree, const std::vector<const PhyTree*> &trees, std::ostream &out);

#endif /* STOCKHOLM_H_ */
