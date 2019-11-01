#include <string>
#include <iostream>
#include <algorithm>
#include "debug.h"

#include "ProgressiveAlignment.h"
#include "ModelFactory.h"
#include "SequenceGraph.h"
#include "Alphabet.h"

template <>
ProgressiveAlignmentResult<AA> progressive_alignment<AA>(const std::map<std::string,sequence_t<AA> > &sequences, const PhyTree &tree, const std::map<std::string,std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<AA> &model_factory, std::map<const PhyTree*,ProgressiveAlignmentResult<AA> > &alignment_cache) throw (ProgressiveAlignmentException)
{
	ProgressiveAlignmentResult<AA> result;
	result.is_csprofile = false;

	if(tree.isLeaf()) {
		const std::string name = tree.getName();
		std::map<std::string,sequence_t<AA> >::const_iterator it = sequences.find(name);
		std::map<std::string,std::vector<repeat_t> >::const_iterator it2 = repeats.find(name);

		if(it != sequences.end()) {
			if(csprofile) {
				result.graph = SequenceGraph<AA>(it->second,*csprofile,model_factory.getModel(tree.getBranchLength()));
				result.is_csprofile = true;
			} else {
				result.graph = SequenceGraph<AA>(it->second);
			}
			result.aligned_sequences[name] = it->second;
			result.profiles[name] = result.graph.getSites().block(0,1,AA::DIM,result.graph.size()-2);

			if(it2 != repeats.end()) {
				for(std::vector<repeat_t>::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3) {
					std::vector<int> tr_hom(result.graph.size(),-1);
					std::copy(it3->tr_hom.begin(), it3->tr_hom.end(), tr_hom.begin()+it3->start+1);
					result.tr_homologies.push_back(tr_hom);
					result.tr_source.push_back(name);
				}
				result.graph.addRepeats(result.tr_homologies);
			}

			result.score = 0;
			result.n_tr_indels = 0;
		} else {
				error("unknown sequence name: %s",name.c_str());
		}
	} else {
		if(tree.n_children() != 2) error("only bifurcating trees allowed");

		ProgressiveAlignmentResult<AA> r1 = progressive_alignment<AA>(sequences,tree[0],repeats,csprofile,model_factory,alignment_cache);
		ProgressiveAlignmentResult<AA> r2 = progressive_alignment<AA>(sequences,tree[1],repeats,csprofile,model_factory,alignment_cache);

		double distance1 = tree[0].getBranchLength();
		double distance2 = tree[1].getBranchLength();

		double support1 = tree[0].getBranchSupport();
		double support2 = tree[1].getBranchSupport();

		result = align_progressive_results<AA>(r1,r2,distance1,distance2,support1,support2,model_factory);

		if(cmdlineopts.earlyref_flag) {
			result = earlyRefinement<AA>(result,tree,model_factory,alignment_cache);
		}
	}

	if(cmdlineopts.earlyref_flag) {
		alignment_cache[&tree] = result;
	}

	return result;
}
