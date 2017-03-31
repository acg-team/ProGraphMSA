#include "main.h"
#include "Stockholm.h"
#include "DateCompiled.h"

void write_stockholm(const std::map<std::string,std::string > alignment, std::ostream &out) {
	out << "# STOCKHOLM 1.0" << std::endl;
	out << "# created by ProGraphMSA " << LogoDate << std::endl;
	for(std::map<std::string,std::string>::const_iterator it=alignment.begin(); it != alignment.end(); ++it) {
		out << it->first << "\t" << it->second << std::endl;
	}
	out << "//" << std::endl;
}

void write_stockholm(const std::map<std::string,std::string > alignment, const std::vector<std::string> &order, const PhyTree &tree, std::ostream &out) {
	out << "# STOCKHOLM 1.0" << std::endl;
	out << "# created by ProGraphMSA " << LogoDate << std::endl;
	for(std::vector<std::string>::const_iterator it=order.begin(); it != order.end(); ++it) {
		out << *it << "\t" << alignment.at(*it) << std::endl;
	}
	out << "#=GF TN guide_tree" << std::endl;
	out << "#=GF NH\t" << tree.formatNewick() << std::endl;
	out << "//" << std::endl;
}

void write_stockholm(const std::map<std::string,std::string > alignment, const std::vector<std::string> &order, const PhyTree &tree, const std::vector<const PhyTree*> &trees, std::ostream &out) {
	out << "# STOCKHOLM 1.0" << std::endl;
	out << "# created by ProGraphMSA " << LogoDate << std::endl;
	for(std::vector<std::string>::const_iterator it=order.begin(); it != order.end(); ++it) {
		out << *it << "\t" << alignment.at(*it) << std::endl;
	}
	for(index_t i=0; i < trees.size(); ++i) {
		out << "#=GF TN guide_tree_iteration_" << i << std::endl;
		out << "#=GF NH\t" << trees[i]->formatNewick() << std::endl;
	}
	out << "#=GF TN guide_tree" << std::endl;
	out << "#=GF NH\t" << tree.formatNewick() << std::endl;
	out << "//" << std::endl;
}
