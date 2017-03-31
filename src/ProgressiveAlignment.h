#ifndef PROGRESSIVEALIGNMENT_H_
#define PROGRESSIVEALIGNMENT_H_

template <class ALPHABET>
struct ProgressiveAlignmentResult;

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "Alphabet.h"
#include "GraphAlign.h"
#include "Graph.h"
#include "PhyTree.h"
#include "ModelFactory.h"
#include "debug.h"
#include "CSProfile.h"
#include "Repeat.h"
#include "SequenceGraph.h"

#undef SAMPLE_ANCESTRAL_SEQUENCES

template <class ALPHABET>
struct ProgressiveAlignmentResult {
	std::map<std::string,sequence_t<ALPHABET> > aligned_sequences;
	std::map<std::string,typename Model<ALPHABET>::Profile> profiles;
	std::vector<std::vector<int> > tr_homologies;
	std::vector<std::string> tr_source;
	Graph<ALPHABET> graph;
	score_t score;
	index_t n_tr_indels;
	bool is_csprofile;
};

class ProgressiveAlignmentException : public swps3_exception {
public:
	ProgressiveAlignmentException(std::string msg) : swps3_exception(msg) {}
};

template <class ALPHABET>
ProgressiveAlignmentResult<ALPHABET> align_progressive_results(const ProgressiveAlignmentResult<ALPHABET> &r1, const ProgressiveAlignmentResult<ALPHABET> &r2, double distance1, double distance2, double support1, double support2, const ModelFactory<ALPHABET> &model_factory);
template <class ALPHABET>
static ProgressiveAlignmentResult<ALPHABET> earlyRefinement(const ProgressiveAlignmentResult<ALPHABET> &old_result, const PhyTree &tree, const ModelFactory<ALPHABET> &model_factory, const std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > &alignment_cache);
template <class ALPHABET>
ProgressiveAlignmentResult<ALPHABET> progressive_alignment(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const PhyTree &tree, const std::map<std::string,std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<ALPHABET> &model_factory, std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > &alignment_cache) throw (ProgressiveAlignmentException);

/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

template <class ALPHABET>
ProgressiveAlignmentResult<ALPHABET> progressive_alignment(const std::map<std::string,sequence_t<ALPHABET> > &sequences, const PhyTree &tree, const std::map<std::string,std::vector<repeat_t> > &repeats, const CSProfile *csprofile, const ModelFactory<ALPHABET> &model_factory, std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > &alignment_cache) throw (ProgressiveAlignmentException)
{
	ProgressiveAlignmentResult<ALPHABET> result;
	result.is_csprofile = false;

	if(tree.isLeaf()) {
		const std::string name = tree.getName();
		typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator it = sequences.find(name);
		typename std::map<std::string,std::vector<repeat_t> >::const_iterator it2 = repeats.find(name);

		if(it != sequences.end()) {
			result.graph = SequenceGraph<ALPHABET>(it->second);
			result.aligned_sequences[name] = it->second;
			result.profiles[name] = result.graph.getSites().block(0,1,ALPHABET::DIM,result.graph.size()-2);
			result.score = 0;
			result.n_tr_indels = 0;

			if(it2 != repeats.end()) {
				for(std::vector<repeat_t>::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3) {
					std::vector<int> tr_hom(result.graph.size(),-1);
					std::copy(it3->tr_hom.begin(), it3->tr_hom.end(), tr_hom.begin()+it3->start+1);
					result.tr_homologies.push_back(tr_hom);
					result.tr_source.push_back(name);
				}
				result.graph.addRepeats(result.tr_homologies);
			}
		} else {
			error("unknown sequence name: %s",name.c_str());
		}
	} else {
		if(tree.n_children() != 2) error("only bifurcating trees allowed");

		ProgressiveAlignmentResult<ALPHABET> r1 = progressive_alignment<ALPHABET>(sequences,tree[0],repeats,csprofile,model_factory,alignment_cache);
		ProgressiveAlignmentResult<ALPHABET> r2 = progressive_alignment<ALPHABET>(sequences,tree[1],repeats,csprofile,model_factory,alignment_cache);

		double distance1 = tree[0].getBranchLength();
		double distance2 = tree[1].getBranchLength();

		double support1 = tree[0].getBranchSupport();
		double support2 = tree[1].getBranchSupport();

		result = align_progressive_results<ALPHABET>(r1,r2,distance1,distance2,support1,support2,model_factory);

		if(cmdlineopts.earlyref_flag) {
			result = earlyRefinement<ALPHABET>(result,tree,model_factory,alignment_cache);
		}
	}

	if(cmdlineopts.earlyref_flag) {
		alignment_cache[&tree] = result;
	}

	return result;
}

template <class ALPHABET>
static ProgressiveAlignmentResult<ALPHABET> earlyRefinement(const ProgressiveAlignmentResult<ALPHABET> &old_result, const PhyTree &tree, const ModelFactory<ALPHABET> &model_factory, const std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > &alignment_cache)
{
	ProgressiveAlignmentResult<ALPHABET> old_results[4];
	double distances[4];
	double gap_distances[4];
	index_t n_results = 0;

	assert(tree.n_children() == 2);

	if(tree[0].isLeaf() && tree[1].isLeaf()) {
		return old_result;
	} else {
		for(index_t i=0; i < tree.n_children(); ++i) {
			if(tree[i].isLeaf()) {
				old_results[n_results] = alignment_cache.at(&tree[i]);
				gap_distances[n_results] = tree[i].getBranchLength();
				if(old_results[n_results].is_csprofile) {
					distances[n_results] = 0;
				} else {
					distances[n_results] = tree[i].getBranchLength();
				}
				++n_results;
			} else {
				const PhyTree &parent=tree[i];
				assert(parent.n_children() == 2);

				for(index_t j=0; j < tree.n_children(); ++j) {
					old_results[n_results] = alignment_cache.at(&parent[j]);
					gap_distances[n_results] = parent.getBranchLength() + parent[j].getBranchLength();
					distances[n_results] = parent.getBranchLength();
					if(!old_results[n_results].is_csprofile) {
						distances[n_results] += parent[j].getBranchLength();
					}
					++n_results;
				}
			}
		}
		assert(n_results >= 2);

		ProgressiveAlignmentResult<ALPHABET> result;
		result.score = old_result.score;
		result.is_csprofile = false;
		result.n_tr_indels = old_result.n_tr_indels;

		std::vector<index_t> mappings[4];
		Graph<ALPHABET> anc_graph = old_result.graph;
		anc_graph.reset();
		std::vector<index_t> anc_mapping(anc_graph.size(),-2);
		for(index_t i=0; i < anc_graph.size(); ++i) {
			anc_mapping[i] = i;
		}

		for(index_t i = 0; i < n_results; ++i) {
			Model<ALPHABET> model = model_factory.getModel(distances[i],gap_distances[i]);
			AlignmentResult<ALPHABET> alignment_result = alignGraphs<ALPHABET>(old_result.graph,old_results[i].graph,model);
			for(index_t j=0; j < alignment_result.mapping1.size(); ++j) {
				if(alignment_result.mapping1[j] != (index_t)-1)
					alignment_result.mapping1[j] = anc_mapping[alignment_result.mapping1[j]];
			}

			AncestralResult<ALPHABET> ancestral_result = mergeGraphsIncremental<ALPHABET>(anc_graph,old_results[i].graph,alignment_result.mapping1,alignment_result.mapping2,model);
			anc_graph = ancestral_result.graph;
			mappings[i] = ancestral_result.mapping2;

			/* update anc_mapping: old_graph -> new anc_graph, new graph should be superset */
			std::vector<index_t> inv_mapping(anc_graph.size(),-2);
			for(index_t j=0; j < ancestral_result.mapping1.size(); ++j) {
				if(ancestral_result.mapping1[j] != (index_t)-1)
					inv_mapping[ancestral_result.mapping1[j]] = j;
			}
			for(index_t j=0; j < anc_mapping.size(); ++j) {
				anc_mapping[j] = inv_mapping[anc_mapping[j]];
			}

			for(index_t j = 0; j < i; ++j) {
				std::vector<index_t> new_mapping(anc_graph.size(),-2);
				for(index_t k = 0; k < anc_graph.size(); ++k) {
					new_mapping[k] = ancestral_result.mapping1[k];
					if(new_mapping[k] != (index_t)-1)
						new_mapping[k] = mappings[j][new_mapping[k]];
				}
				mappings[j] = new_mapping;
			}
		}

		/* remove unused nodes in anc_graph (and according mappings) */
		for(index_t i = 0; i < anc_graph.size(); ++i) {
			bool node_used = false;
			for(index_t j = 0; j < n_results; ++j) {
				if(mappings[j][i] != (index_t)-1) {
					node_used = true;
					break;
				}
			}

			if(!node_used) {
				for(index_t j = i+1; j < anc_graph.size(); ++j) {
					bool node_used2 = false;
					for(index_t k = 0; k < n_results; ++k) {
						if(mappings[k][j] != (index_t)-1) {
							node_used2 = true;
							break;
						}
					}
					if(node_used2) {
						anc_graph.rmNodes(i,j-i);
						for(index_t k = 0; k < n_results; ++k) {
							mappings[k].erase(mappings[k].begin()+i,mappings[k].begin()+j);
						}
						--i;
						break;
					}
				}
			}
		}

		result.graph = anc_graph;

		for(index_t i = 0; i < n_results; ++i) {
			extend_alignment(result,mappings[i],old_results[i].aligned_sequences);
			extend_tr_homologies(result,mappings[i],old_results[i].tr_homologies,old_results[i].tr_source);
		}

		result.graph.addRepeats(result.tr_homologies);

		return result;
	}
}

template <class ALPHABET>
static void extend_alignment(ProgressiveAlignmentResult<ALPHABET> &result, const std::vector<index_t> &mapping, const std::map<std::string,sequence_t<ALPHABET> > &aligned_sequences)
{
	for(typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator it=aligned_sequences.begin(); it != aligned_sequences.end();++it) {
		sequence_t<ALPHABET> extended(result.graph.size()-2,ALPHABET::X);
		const sequence_t<ALPHABET> &original = it->second;

		index_t k = 0;
		for(index_t j=1;j<result.graph.size()-1;++j) {
			if(mapping[j] != (index_t)-1) {
				assert(k == mapping[j]-1);
				extended[j-1] = original[k++];
			} else {
				extended[j-1] = ALPHABET::GAP;
			}
		}

		result.aligned_sequences[it->first] = extended;
	}
}

template <class ALPHABET>
static void extend_tr_homologies(ProgressiveAlignmentResult<ALPHABET> &result, const std::vector<index_t> &mapping, const std::vector<std::vector<int> > &tr_homologies, const std::vector<std::string> &tr_source)
{
	std::vector<std::string>::const_iterator source = tr_source.begin();
	for(std::vector<std::vector<int> >::const_iterator it=tr_homologies.begin(); it != tr_homologies.end();++it,++source) {
		std::vector<int> extended(result.graph.size()-2,-1);
		const std::vector<int> &original = *it;

		index_t k = 0;
		for(index_t j=1;j<result.graph.size()-1;++j) {
			if(mapping[j] != (index_t)-1) {
				assert(k == mapping[j]-1);
				extended[j-1] = original[k++];
			} else {
				extended[j-1] = -1;
			}
		}

		result.tr_homologies.push_back(extended);
		result.tr_source.push_back(*source);
	}
}

template <class ALPHABET>
static std::string create_ancestral_seq_name(const std::map<std::string,sequence_t<ALPHABET> > &aligned_seqs)
{
	std::vector<std::string> leaves;
	leaves.reserve(aligned_seqs.size());

	for(typename std::map<std::string,sequence_t<ALPHABET> >::const_iterator it=aligned_seqs.begin(); it != aligned_seqs.end();++it) {
		const std::string &name = it->first;
		if(name[0] != '(') {
			leaves.push_back(name);
		}
	}

	std::sort(leaves.begin(),leaves.end());

	std::stringstream ss;
	ss << "(";
	for(std::vector<std::string>::const_iterator it=leaves.begin(); it != leaves.end();++it) {
		const std::string &name = *it;
		if(it != leaves.begin()) {
			ss << ",";
		}
		ss << name;
	}
	ss << ")";
	return ss.str();
}

template <class ALPHABET>
static void prelim_ancestral_seq(ProgressiveAlignmentResult<ALPHABET> &result, const std::vector<bool> &is_matched, const Model<ALPHABET> &model)
{
	sequence_t<ALPHABET> extended(result.graph.size()-2,ALPHABET::X);
	std::string anc_name = create_ancestral_seq_name(result.aligned_sequences);
	index_t anc_length = 0;

	for(index_t i=1;i<result.graph.size()-1;++i) {
		if(is_matched[i]) {
			typename Model<ALPHABET>::Freqs col = result.graph[i];
			col.array() *= model.pi.array();

#ifdef SAMPLE_ANCESTRAL_SEQUENCES
			col /= col.sum();
			double rnd = drand48();
			int j=0;
			for(; j<ALPHABET::DIM-1; ++j) {
				rnd -= col(j);
				if(rnd <= 0) break;
			}
#else
			typename Model<ALPHABET>::Freqs::Index j;
			col.maxCoeff(&j);
#endif

			extended[i-1] = ALPHABET((int)j);
			anc_length += 1;
		} else {
			extended[i-1] = ALPHABET::GAP;
		}
	}

	result.aligned_sequences[anc_name] = extended;


	typename Model<ALPHABET>::Profile profile = Model<ALPHABET>::Profile::Zero(ALPHABET::DIM,anc_length);
	for(index_t i=1,j=0;i<result.graph.size()-1;++i) {
		if(is_matched[i]) {
			typename Model<ALPHABET>::Freqs col = result.graph[i];
			col.array() *= model.pi.array();
			col /= col.sum();
			profile.col(j++) = col;
		}
	}

	result.profiles[anc_name] = profile;
}


template <class ALPHABET>
static void final_ancestral_seq(ProgressiveAlignmentResult<ALPHABET> &result, const std::vector<index_t> &mapping, const std::vector<bool> &matched, const ProgressiveAlignmentResult<ALPHABET> &old_result, const Model<ALPHABET> &model)
{
	sequence_t<ALPHABET> extended(result.graph.size()-2,ALPHABET::X);
	std::string anc_name = create_ancestral_seq_name(old_result.aligned_sequences);
	index_t anc_length = 0;

	for(index_t i=1;i<result.graph.size()-1;++i) {
		if(matched[i] && mapping[i] != (index_t)-1) {
			typename Model<ALPHABET>::Freqs col = old_result.graph[mapping[i]];
			col.array() *= model.pi.array();

#ifdef SAMPLE_ANCESTRAL_SEQUENCES
			col /= col.sum();
			double rnd = drand48();
			int j=0;
			for(; j<ALPHABET::DIM-1; ++j) {
				rnd -= col(j);
				if(rnd <= 0) break;
			}
#else
			typename Model<ALPHABET>::Freqs::Index j;
			col.maxCoeff(&j);
#endif

			extended[i-1] = ALPHABET((int)j);
			anc_length += 1;
		} else {
			extended[i-1] = ALPHABET::GAP;
		}
	}

	result.aligned_sequences[anc_name] = extended;

	typename Model<ALPHABET>::Profile profile = Model<ALPHABET>::Profile::Zero(ALPHABET::DIM,anc_length);
	for(index_t i=1,j=0;i<result.graph.size()-1;++i) {
		if(matched[i] && mapping[i] != (index_t)-1) {
			typename Model<ALPHABET>::Freqs col = old_result.graph[mapping[i]];
			col.array() *= model.pi.array();
			col /= col.sum();
			profile.col(j++) = col;
		}
	}

	result.profiles[anc_name] = profile;
}

template <class ALPHABET>
ProgressiveAlignmentResult<ALPHABET> align_progressive_results(const ProgressiveAlignmentResult<ALPHABET> &r1, const ProgressiveAlignmentResult<ALPHABET> &r2, double distance1, double distance2, double support1, double support2, const ModelFactory<ALPHABET> &model_factory)
{
	ProgressiveAlignmentResult<ALPHABET> result;

	double gap_distance1 = distance1;
	double gap_distance2 = distance2;

	if(r1.is_csprofile) {
		distance1 = 0;
	}

	if(r2.is_csprofile) {
		distance2 = 0;
	}

	double gap_distance = gap_distance1 + gap_distance2;
	double distance = distance1 + distance2;

	Model<ALPHABET> model = model_factory.getModel(distance,gap_distance);
	Model<ALPHABET> model1 = model_factory.getModel(distance1,gap_distance1);
	Model<ALPHABET> model2 = model_factory.getModel(distance2,gap_distance2);

	CleanedGraph<ALPHABET> cg1(r1.graph);
	CleanedGraph<ALPHABET> cg2(r2.graph);

	AlignmentResult<ALPHABET> alignment_result = alignGraphs<ALPHABET>(cg1,cg2,model);
	result.score = alignment_result.score;
	result.is_csprofile = false;
	result.n_tr_indels = alignment_result.n_tr_indels + r1.n_tr_indels + r2.n_tr_indels;
	result.profiles.insert(r1.profiles.begin(), r1.profiles.end());
	result.profiles.insert(r2.profiles.begin(), r2.profiles.end());

	cg1.uncleanMapping(alignment_result.mapping1);
	cg2.uncleanMapping(alignment_result.mapping2);

	AncestralResult<ALPHABET> ancestral_result = mergeGraphs(r1.graph,r2.graph,alignment_result.mapping1,alignment_result.mapping2,model1,model2,support1,support2);
	result.graph = ancestral_result.graph;

	extend_alignment(result,ancestral_result.mapping1,r1.aligned_sequences);
	extend_alignment(result,ancestral_result.mapping2,r2.aligned_sequences);

	extend_tr_homologies(result,ancestral_result.mapping1,r1.tr_homologies,r1.tr_source);
	extend_tr_homologies(result,ancestral_result.mapping2,r2.tr_homologies,r2.tr_source);

	if(cmdlineopts.ancestral_flag) {
		if(r1.aligned_sequences.size() > 1) {
			final_ancestral_seq(result,ancestral_result.mapping1,ancestral_result.is_matched,r1,model1);
		}
		if(r2.aligned_sequences.size() > 1) {
			final_ancestral_seq(result,ancestral_result.mapping2,ancestral_result.is_matched,r2,model2);
		}
		prelim_ancestral_seq(result,ancestral_result.is_matched,model);
	}

	result.graph.addRepeats(result.tr_homologies);

	if(cmdlineopts.repeats_flag) {
		std::cerr << "TR indels at " << create_ancestral_seq_name(result.aligned_sequences)
			<< ": " <<  alignment_result.n_tr_indels << std::endl;
	}

	return result;
}

#endif /* PROGRESSIVEALIGNMENT_H_ */
