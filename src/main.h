/*
 * Copyright (c) 2007-2011 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef MAIN_H
#define MAIN_H

#include <string>

typedef double score_t;
typedef float dp_score_t;
typedef unsigned int index_t;
typedef double distance_t;


struct cmdlineopts_t {
	std::string output_file;
	std::string sequence_file;
	std::string tree_file;
	std::string topo_file;
	std::string cs_file;
	std::string cmodel_file;
	std::string readreps_file;
	std::string trdout_file;
	std::string profile_file;
	std::string customtr_cmd;
	int iters;
	int reroot_flag;
	int wlsrefine_flag;
	bool earlyref_flag;
	bool repeats_flag;
	bool repalign_flag;
	bool fasta_flag;
	bool noforcealign_flag;
	bool aafreqs_flag;
	bool darwin_flag;
	bool nwdist_flag;
	bool onlytree_flag;
	bool mldist_flag;
	bool mldist_gap_flag;
	bool alltrees_flag;
	bool ancestral_flag;
	bool codon_flag;
	bool dna_flag;
	bool inputorder_flag;
	double indel_rate;
	double end_indel_prob;
	double gapext_prob;
	double edge_halflife;
	double altsplice_prob;
	double pseudo_count;
	double cutoff_dist;
	double repeat_rate;
	double repeatext_prob;
	double max_dist;
	double min_dist;
	double max_pdist;
	double min_pdist;
};

extern cmdlineopts_t cmdlineopts;

#endif /* MAIN_H */
