#include <tclap/CmdLine.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Fasta.h"
#include "Stockholm.h"
#include "ProgressiveAlignment.h"
#include "ModelFactoryWag.h"
#include "ModelFactoryDarwin.h"
#include "ModelFactoryEcm.h"
#include "ModelFactoryCustom.h"
#include "ModelFactoryPlusF.h"
#include "newick.h"
#include "profile.h"
#include "TreeNJ.h"
#include "FindRoot.h"
#include "DateCompiled.h"

#include "CSProfile.h"
#include "RepeatDetectionTReks.h"
#include "Alphabet.h"

#include "main.h"

// ugly global variable
cmdlineopts_t cmdlineopts;

template <class ALPHABET>
void doAlign(const std::map<std::string,std::string> &seqs, std::map<std::string,std::string> &out_aligned_seqs, std::vector<const PhyTree*> &out_all_trees);

int main(int argc, char** argv)
{
	//command line parsing
	try {

		TCLAP::CmdLine cmd(
				"ProGraphMSA, fast multiple sequence alignment", ' ', LogoDate);

		TCLAP::ValueArg<std::string> outputArg("o", "output", "Output file name",
				false, "/dev/null", "filename");
		cmd.add(outputArg);

		TCLAP::ValueArg<std::string> treeArg("t", "tree", "initial guide tree", false, "", "newick file");
		cmd.add(treeArg);

		TCLAP::ValueArg<std::string> topoArg("", "topology", "topology of initial guide tree (branch lengths will be estimated)", false, "", "newick file");
		cmd.add(topoArg);

#ifdef WITH_CODON
		TCLAP::SwitchArg codonArg("","codon","align DNA sequence based on a codon model");
		cmd.add(codonArg);
#endif

#ifdef WITH_DNA
		TCLAP::SwitchArg dnaArg("","dna","align DNA sequence");
		cmd.add(dnaArg);
#endif

		TCLAP::SwitchArg fastaArg("f","fasta","output fasta format (instead of stockholm)");
		cmd.add(fastaArg);

		TCLAP::ValueArg<double> indelArg("g", "indel_rate", "insertion/deletion rate", false, 0.0093359375, "rate");
		cmd.add(indelArg);

		TCLAP::ValueArg<double> gapextArg("e", "gap_ext", "gap extension probability", false, 0.6119140625, "probability");
		cmd.add(gapextArg);

		TCLAP::ValueArg<double> endindelArg("E", "end_indel_prob", "probability of mismatching sequence ends (set to -1 to disable this feature)", false, 0.12, "probability");
		cmd.add(endindelArg);

		TCLAP::ValueArg<double> edgehlArg("l", "edge_halflife", "edge half-life", false, 0.3, "distance");
		cmd.add(edgehlArg);

		TCLAP::ValueArg<double> altspliceArg("s", "altsplice_prob", "alternative splicing probability", false, 0.328125, "probability");
		cmd.add(altspliceArg);

		TCLAP::ValueArg<double> cutdistArg("x", "cutoff_dist", "cutoff value for pairwise distance estimation", false, 2.2, "distance");
		cmd.add(cutdistArg);

		TCLAP::ValueArg<double> mindistArg("d", "min_dist", "minimum distance for alignment", false, 0.05, "distance");
		cmd.add(mindistArg);

		TCLAP::ValueArg<double> maxdistArg("D", "max_dist", "maximum distance for alignment", false, 2.2, "distance");
		cmd.add(maxdistArg);

		TCLAP::ValueArg<double> minpdistArg("p", "min_pdist", "minimum p-distance (divergence) for alignment", false, 0.05, "distance");
		cmd.add(minpdistArg);

		TCLAP::ValueArg<double> maxpdistArg("P", "max_pdist", "maximum p-distance (divergence) for alignment", false, 0.8, "distance");
		cmd.add(maxpdistArg);

		TCLAP::SwitchArg noforcealignArg("A","no_force_align","do not force alignment of start/stop codons");
		cmd.add(noforcealignArg);

		TCLAP::ValueArg<double> reprateArg("", "repeat_indel_rate", "insertion/deletion rate for repeat units (per site)", false, 0.1, "rate");
		cmd.add(reprateArg);

		TCLAP::ValueArg<double> repextArg("", "repeat_indel_ext", "repeat indel extension probability", false, 0.3, "probability");
		cmd.add(repextArg);

		TCLAP::SwitchArg repalignArg("", "repalign", "re-align detected tandem repeat units");
		cmd.add(repalignArg);

		TCLAP::MultiSwitchArg repeatsArg("R", "repeats", "use T-Reks to identify tandem repeats");
		cmd.add(repeatsArg);

		TCLAP::ValueArg<std::string> readrepsArg("", "read_repeats", "read repeats from file", false, "", "T-Reks format output");
		cmd.add(readrepsArg);

		TCLAP::ValueArg<std::string> trdoutputArg("", "trd_output", "write TR detector output to file", false, "", "filename");
		cmd.add(trdoutputArg);

		TCLAP::ValueArg<std::string> customtrcmdArg("", "custom_tr_cmd", "custom command for detecting tandem-repeats", false, "", "command");
		cmd.add(customtrcmdArg);

		TCLAP::MultiSwitchArg rerootArg("r", "reroot", "reroot tree on all branches and minimize gap parsimony (specify twice for heuristic root search)");
		cmd.add(rerootArg);

		TCLAP::MultiSwitchArg wlsrefineArg("W","wls_refine","refine guide tree with weighted least-squares");
		cmd.add(wlsrefineArg);

		TCLAP::SwitchArg earlyrefArg("","early_refinement","perform early refinement");
		cmd.add(earlyrefArg);

		TCLAP::ValueArg<std::string> csArg("c", "cs_profile", "library of context-sensitive profiles", false, "", "file");
		cmd.add(csArg);

		TCLAP::SwitchArg darwinArg("w","darwin","use darwin's model of evolution (instead of WAG)");
		cmd.add(darwinArg);
		
		TCLAP::ValueArg<std::string> profileArg("", "profile_out", "output ancestral sequence profiles", false, "", "file");
		cmd.add(profileArg);

		TCLAP::ValueArg<std::string> cmodelArg("", "custom_model", "custom substitution model in qmat format", false, "", "file");
		cmd.add(cmodelArg);

		TCLAP::SwitchArg aafreqsArg("F", "estimate_aafreqs", "estimate equilibrium amino acid frequencies from input data");
		cmd.add(aafreqsArg);

		TCLAP::ValueArg<double> pcountArg("C", "aafreqs_pseudocount", "pseudo-count for estimating equilibrium amino acid frequencies", false, 1000.0, "count");
		cmd.add(pcountArg);

		TCLAP::SwitchArg nwdistArg("a","nwdist","estimate initial distance tree from NW alignments");
		cmd.add(nwdistArg);

		TCLAP::SwitchArg mldistArg("m","mldist","use ML distances");
		cmd.add(mldistArg);

		TCLAP::SwitchArg mldistgapArg("M","mldist_gap","use ML distances with gaps");
		cmd.add(mldistgapArg);

		TCLAP::SwitchArg inputorderArg("I","input_order","output sequences in input order (default: tree order)");
		cmd.add(inputorderArg);

		TCLAP::SwitchArg onlytreeArg("T","only_tree","only output final tree instead of MSA");
		cmd.add(onlytreeArg);

		TCLAP::ValueArg<int> iterArg("i","iterations","number of iterations re-estimating guide tree", false, 2, "iterations");
		cmd.add(iterArg);

		TCLAP::SwitchArg alltreesArg("","all_trees","output all intermediate guide trees");
		cmd.add(alltreesArg);

		TCLAP::SwitchArg ancestralArg("","ancestral_seqs","output all ancestral sequences");
		cmd.add(ancestralArg);

		TCLAP::UnlabeledValueArg<std::string> seqArg("sequences", "input sequences", true, "", "fasta file");
		cmd.add(seqArg);

		cmd.parse(argc, argv);

		cmdlineopts.output_file = outputArg.getValue();
		cmdlineopts.sequence_file = seqArg.getValue();
		cmdlineopts.tree_file = treeArg.getValue();
		cmdlineopts.topo_file = topoArg.getValue();
		cmdlineopts.cs_file = csArg.getValue();
		cmdlineopts.cmodel_file = cmodelArg.getValue();
		cmdlineopts.profile_file = profileArg.getValue();
		cmdlineopts.fasta_flag = fastaArg.getValue();
		cmdlineopts.indel_rate = indelArg.getValue();
		cmdlineopts.end_indel_prob = endindelArg.getValue();
		cmdlineopts.gapext_prob = gapextArg.getValue();
		cmdlineopts.edge_halflife = edgehlArg.getValue();
		cmdlineopts.altsplice_prob = altspliceArg.getValue();
		cmdlineopts.pseudo_count = pcountArg.getValue();
		cmdlineopts.darwin_flag = darwinArg.getValue();
		cmdlineopts.repeat_rate = reprateArg.getValue();
		cmdlineopts.repeatext_prob = repextArg.getValue();
		cmdlineopts.repalign_flag = repalignArg.getValue();
		cmdlineopts.repeats_flag = repeatsArg.getValue();
		cmdlineopts.readreps_file = readrepsArg.getValue();
		cmdlineopts.trdout_file = trdoutputArg.getValue();
		cmdlineopts.customtr_cmd = customtrcmdArg.getValue();
		cmdlineopts.reroot_flag = rerootArg.getValue();
		cmdlineopts.earlyref_flag = earlyrefArg.getValue();
		cmdlineopts.nwdist_flag = nwdistArg.getValue();
		cmdlineopts.mldist_flag = mldistArg.getValue();
		cmdlineopts.mldist_gap_flag = mldistgapArg.getValue();
		cmdlineopts.onlytree_flag = onlytreeArg.getValue();
		cmdlineopts.alltrees_flag = alltreesArg.getValue();
		cmdlineopts.ancestral_flag = ancestralArg.getValue();
		cmdlineopts.wlsrefine_flag = wlsrefineArg.getValue();
		cmdlineopts.aafreqs_flag = aafreqsArg.getValue();
		cmdlineopts.noforcealign_flag = noforcealignArg.getValue();
		cmdlineopts.inputorder_flag = inputorderArg.getValue();
		cmdlineopts.min_dist = mindistArg.getValue();
		cmdlineopts.max_dist = maxdistArg.getValue();
		cmdlineopts.min_pdist = minpdistArg.getValue();
		cmdlineopts.max_pdist = maxpdistArg.getValue();
		cmdlineopts.cutoff_dist = cutdistArg.getValue();
		cmdlineopts.iters = iterArg.getValue();

#ifdef WITH_CODON
		cmdlineopts.codon_flag = codonArg.isSet();
#else
		cmdlineopts.codon_flag = false;
#endif
#ifdef WITH_DNA
		cmdlineopts.dna_flag = dnaArg.isSet();
#else
		cmdlineopts.dna_flag = false;
#endif

#ifdef WITH_CODON
		/* scale default parameters for codon distances */
		if(codonArg.isSet()) {
			if(!indelArg.isSet()) {
				cmdlineopts.indel_rate /= 2.6;
			}
			if(!edgehlArg.isSet()) {
				cmdlineopts.edge_halflife *= 2.6;
			}
			if(!maxdistArg.isSet()) {
				cmdlineopts.max_dist = 5.0;
			}
			if(!cutdistArg.isSet()) {
				cmdlineopts.cutoff_dist = 5.0;
			}
		}
#endif

		/* do not iterate when guide tree provided */
		if(cmdlineopts.tree_file != "" && !iterArg.isSet()) {
			cmdlineopts.iters = 0;
		}

		std::vector<std::string> input_order;
		std::map<std::string,std::string> seqs = FastaLib(cmdlineopts.sequence_file.c_str()).readAll(input_order);

		std::ofstream custom_out;
		std::ostream *out;

		if(outputArg.isSet()) {
			custom_out.open(cmdlineopts.output_file.c_str());
			out = &custom_out;
		} else {
			out = &std::cout;
		}
		if(!*out) {
			error("error opening output file");
		}

		std::map<std::string,std::string> aligned_seqs;
		std::vector<const PhyTree *> all_trees;

		if(!cmdlineopts.codon_flag && !cmdlineopts.dna_flag) {
			doAlign<AA>(seqs,aligned_seqs,all_trees);
		}
#ifdef WITH_CODON
		else if(cmdlineopts.codon_flag) {
			doAlign<Codon>(seqs,aligned_seqs,all_trees);
		}
#endif
#ifdef WITH_DNA
		else if(cmdlineopts.dna_flag) {
			doAlign<DNA>(seqs,aligned_seqs,all_trees);
		}
#endif

		if(!cmdlineopts.onlytree_flag) {
			std::vector<std::string> order = input_order;
			if(!cmdlineopts.inputorder_flag) {
				order = get_tree_order(all_trees.back());
			}

			if(cmdlineopts.fasta_flag) {
				write_fasta(aligned_seqs,order,*out);
			} else {
				if(!cmdlineopts.alltrees_flag) {
					write_stockholm(aligned_seqs,order,*all_trees.back(),*out);
				} else {
					write_stockholm(aligned_seqs,order,*all_trees.back(),all_trees,*out);
				}
			}
		} else {
			if(cmdlineopts.alltrees_flag) {
				for(std::vector<const PhyTree*>::const_iterator it=all_trees.begin(); it < all_trees.end(); ++it) {
					*out << (*it)->formatNewick() << std::endl;
				}
			} else {
				*out << all_trees.back()->formatNewick() << std::endl;
			}
		}

		for(std::vector<const PhyTree*>::const_iterator it=all_trees.begin(); it < all_trees.end(); ++it) {
			delete *it;
		}
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "Command line error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

	return 0;
}

template <class ALPHABET>
void doAlign(const std::map<std::string,std::string> &seqs, std::map<std::string,std::string> &out_aligned_seqs, std::vector<const PhyTree*> &out_all_trees) {
	/* remove start and stop codons */
	bool any_startStripped = false;
	bool any_endStripped = false;
	std::map<std::string,bool> startStripped;
	std::map<std::string,bool> endStripped;

	typename std::map<std::string,sequence_t<ALPHABET> > seqs2;
	for (std::map<std::string, std::string>::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
		seqs2[it->first] = sequenceFromString<ALPHABET>(it->second);
		if(!cmdlineopts.noforcealign_flag) {
			sequence_t<ALPHABET> &seq = seqs2[it->first];
			if(ALPHABET::stripStart != ALPHABET::GAP && seq[0] == ALPHABET::stripStart) {
				seq = seq.substr(1);
				any_startStripped = true;
				startStripped[it->first] = true;
			} else {
				startStripped[it->first] = false;
			}

			if(ALPHABET::stripEnd != ALPHABET::GAP && seq[seq.length()-1] == ALPHABET::stripEnd) {
				seq = seq.substr(0,seq.length()-1);
				any_endStripped = true;
				endStripped[it->first] = true;
			} else {
				endStripped[it->first] = false;
			}
		}
	}


	/* load evolutionary model */
	ModelFactory<ALPHABET> *model_factory = ModelFactory<ALPHABET>::getDefault(seqs2);


	/* load context-specific profile library */
	CSProfile *csprofile = NULL;
	if(cmdlineopts.cs_file != "") {
		csprofile = new CSProfile(cmdlineopts.cs_file);
	}


	/* detect/read repeats */
	std::map< std::string, std::vector<repeat_t> > reps;
	if(cmdlineopts.readreps_file != "")
		reps = read_repeats(cmdlineopts.readreps_file,seqs2);
	else if(cmdlineopts.repeats_flag) {
		reps = detect_repeats(seqs2);

		if(cmdlineopts.repalign_flag) {
			reps = align_repeats(seqs2,reps,csprofile,model_factory);
		}
	}


	/* create/read initial guide tree */
	PhyTree* tree = NULL;
	PhyTree* topo = NULL;

	if(cmdlineopts.topo_file != "") {
		std::ifstream tree_str(cmdlineopts.topo_file.c_str());
		topo = newick_parser::parse_newick(&tree_str);
	}

	if(cmdlineopts.tree_file != "") {
		std::ifstream tree_str(cmdlineopts.tree_file.c_str());
		tree = newick_parser::parse_newick(&tree_str);
	} else {
		tree = TreeNJ<ALPHABET>(seqs2,false,model_factory,topo);
	}

	out_all_trees.push_back(tree->copy());

	ProgressiveAlignmentResult<ALPHABET> result;
	ProgressiveAlignmentResult<ALPHABET> old_result;


	/* further rounds of alignment followed by estimation of
	 * improved tree from induced pairwise distances */
	for(int i = 0; i < cmdlineopts.iters; ++i) {
		std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > alignment_cache;
		result = progressive_alignment<ALPHABET>(seqs2,*tree,reps,csprofile,*model_factory,alignment_cache);

		/* delete ancestral sequences */
		for(typename std::map<std::string,sequence_t<ALPHABET> >::iterator it=result.aligned_sequences.begin(); it != result.aligned_sequences.end();) {
			const std::string &name = it->first;
			if(name[0] == '(') {
				result.aligned_sequences.erase(it++);
			} else {
				++it;
			}
		}

		/* early exit if alignment iteration already converged */
		if(i > 0 && result.aligned_sequences == old_result.aligned_sequences)
			break;

		delete tree;

		tree = TreeNJ<ALPHABET>(result.aligned_sequences,true,model_factory,topo);

		out_all_trees.push_back(tree->copy());

		old_result = result;
	}


	/* compute final alignment */
	if(!cmdlineopts.onlytree_flag) {
		if(cmdlineopts.reroot_flag) {
			result = progressive_alignment_find_root<ALPHABET>(seqs2,*tree,reps,csprofile,*model_factory);
		} else {
			std::map<const PhyTree*,ProgressiveAlignmentResult<ALPHABET> > alignment_cache;
			result = progressive_alignment<ALPHABET>(seqs2,*tree,reps,csprofile,*model_factory,alignment_cache);
		}
	}

	delete tree;
	if(topo) delete topo;
	delete model_factory;
	if(csprofile) delete csprofile;

	if(cmdlineopts.repeats_flag) {
		std::cerr << "TR indels: " << result.n_tr_indels << std::endl;
	}

	if(cmdlineopts.profile_file != "") {
		std::ofstream profile_file;
		profile_file.open(cmdlineopts.profile_file.c_str());
		write_profile<ALPHABET>(result.profiles, profile_file);
		profile_file.close();
	}

	/* re-insert start/stop codons */
	for (typename std::map<std::string, sequence_t<ALPHABET> >::iterator it = result.aligned_sequences.begin(); it != result.aligned_sequences.end(); ++it) {
		sequence_t<ALPHABET> aseq = it->second;
		if(any_startStripped) {
			if(startStripped.find(it->first) != startStripped.end() && startStripped[it->first]) {
				aseq.insert(aseq.begin(),ALPHABET::X);
			} else {
				aseq.insert(aseq.begin(),ALPHABET::GAP);
			}
		}

		if(any_endStripped) {
			if(endStripped.find(it->first) != endStripped.end() && endStripped[it->first]) {
				aseq.insert(aseq.end(),ALPHABET::X);
			} else {
				aseq.insert(aseq.end(),ALPHABET::GAP);
			}
		}

		if(seqs.find(it->first) != seqs.end()) {
			out_aligned_seqs[it->first] = stringFromSequence<ALPHABET>(aseq,seqs.at(it->first));
		} else {
			out_aligned_seqs[it->first] = stringFromSequence<ALPHABET>(aseq);
		}
	}
}
