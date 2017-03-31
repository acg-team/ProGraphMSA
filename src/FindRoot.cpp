#include "FindRoot.h"
#include "GapParsimony.h"
#include "SequenceGraph.h"
#include "debug.h"
#include "main.h"

namespace FindRoot {

	template <>
	void tree2graph<AA>(FindRoot::node<AA> *current, const PhyTree& tree, std::vector<FindRoot::node<AA>*> &nodes, std::vector<FindRoot::edge<AA>*> &edges, const std::map<std::string,sequence_t<AA> > &sequences, const std::map<std::string, std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<AA> &model_factory) {
		current->name = tree.getName();

		if(tree.isLeaf()) {
			std::map<std::string,sequence_t<AA> >::const_iterator it = sequences.find(current->name);
			std::map<std::string, std::vector<repeat_t> >::const_iterator it2 = repeats.find(current->name);

			if(it != sequences.end()) {
				current->cached_alignments[0] = new ProgressiveAlignmentResult<AA>();
				current->cached_alignments[0]->is_csprofile = false;
				if(csprofile) {
					current->cached_alignments[0]->graph = SequenceGraph<AA>(it->second,*csprofile,model_factory.getModel(tree.getBranchLength()));
					current->cached_alignments[0]->is_csprofile = true;
				} else {
					current->cached_alignments[0]->graph = SequenceGraph<AA>(it->second);
				}

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

			FindRoot::edge<AA> *e0 = new FindRoot::edge<AA>();
			edges.push_back(e0);
			e0->length = tree[0].getBranchLength();
			current->edges[1] = e0;

			FindRoot::node<AA> *n0 = new FindRoot::node<AA>();
			nodes.push_back(n0);
			n0->edges[0] = e0;
			e0->nodes[0] = current;
			e0->nodes[1] = n0;

			tree2graph(n0,tree[0],nodes,edges,sequences,repeats,csprofile,model_factory);


			FindRoot::edge<AA> *e1 = new FindRoot::edge<AA>();
			edges.push_back(e1);
			e1->length = tree[1].getBranchLength();
			current->edges[2] = e1;

			FindRoot::node<AA> *n1 = new FindRoot::node<AA>();
			nodes.push_back(n1);
			n1->edges[0] = e1;
			e1->nodes[0] = current;
			e1->nodes[1] = n1;

			tree2graph(n1,tree[1],nodes,edges,sequences,repeats,csprofile,model_factory);
		}
	}
}
