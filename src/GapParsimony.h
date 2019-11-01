#include "FindRoot.h"
#include "ProgressiveAlignment.h"

#ifndef GAPPARSIMONY_H_
#define GAPPARSIMONY_H_

namespace GapParsimony {
	template <class ALPHABET>
	score_t scoreAlignment(const ProgressiveAlignmentResult<ALPHABET> &alignment, FindRoot::edge<ALPHABET> *root);
}

/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include <iostream>
#include <string>
#include "newick.h"

#include "main.h"

namespace GapParsimony {
	typedef unsigned long block_t;

	struct result {
		block_t* consensus;
		unsigned int score;
	};

	template <class ALPHABET>
	static result scoreSubtree(const ProgressiveAlignmentResult<ALPHABET> &alignment, FindRoot::node<ALPHABET> *subtree, FindRoot::edge<ALPHABET> *other)
	{
		result res;
		const index_t length = alignment.aligned_sequences.begin()->second.length();
		const index_t cperblock = sizeof(block_t)*4;
		const index_t nblocks = length/cperblock + (length%cperblock != 0);

		if(subtree->isLeaf()){
			res.consensus = new block_t[nblocks];
			const sequence_t<ALPHABET> &seq = alignment.aligned_sequences.at(subtree->name);
			assert(length == seq.length());

			for(index_t i=0; i<nblocks; ++i) {
				res.consensus[i] = 0UL;
			}

			/* initialize bitsets (lower bit set if character, high bit if gap) */
			for(index_t i=0; i<length; ++i) {
				index_t block = i/cperblock;
				index_t pos = i%cperblock;
				res.consensus[block] |= 1UL << (pos*2+(seq[i].isGap()));
			}

			/* fill remaining part with ones -> not counted as mismatches */
			for(index_t i=length%cperblock; i<cperblock; ++i) {
				res.consensus[nblocks-1] |= 3UL << (i*2);
			}


			res.score = 0;
		} else {
			FindRoot::edge<ALPHABET> *edge1 = (subtree->edges[0] == other) ? subtree->edges[1] : subtree->edges[0];
			FindRoot::edge<ALPHABET> *edge2 = (subtree->edges[2] == other) ? subtree->edges[1] : subtree->edges[2];
			FindRoot::node<ALPHABET> *node1 = (edge1->nodes[0] == subtree) ? edge1->nodes[1] : edge1->nodes[0];
			FindRoot::node<ALPHABET> *node2 = (edge2->nodes[0] == subtree) ? edge2->nodes[1] : edge2->nodes[0];

			result res1 = scoreSubtree<ALPHABET>(alignment,node1,edge1);
			result res2 = scoreSubtree<ALPHABET>(alignment,node2,edge2);

			res.score = res1.score + res2.score;
			res.consensus = res1.consensus;
			for(index_t i=0; i<nblocks; ++i) {
				res.consensus[i] &= res2.consensus[i];
				block_t tmp1 = ~res.consensus[i];
				tmp1 = tmp1 & (tmp1 << 1) & (~0UL/3UL*2UL); // == 0xAAAAAAAAAAAAAAAAUL

				res.score += __builtin_popcountl(tmp1);

				tmp1 |= tmp1 >> 1;
				res.consensus[i] |= tmp1;
			}

			delete [] res2.consensus;
		}

		return res;
	}

	template <class ALPHABET>
	score_t scoreAlignment(const ProgressiveAlignmentResult<ALPHABET> &alignment, typename FindRoot::edge<ALPHABET> *root)
	{
		const index_t length = alignment.aligned_sequences.begin()->second.length();
		const index_t cperblock = sizeof(block_t)*4;
		const index_t nblocks = length/cperblock + (length%cperblock != 0);

		GapParsimony::result result1 = GapParsimony::scoreSubtree<ALPHABET>(alignment,root->nodes[0],root);
		GapParsimony::result result2 = GapParsimony::scoreSubtree<ALPHABET>(alignment,root->nodes[1],root);

		unsigned int score = result1.score + result2.score;

		for(index_t i=0; i<nblocks; ++i) {
			block_t consensus = result1.consensus[i] & result2.consensus[i];
			block_t tmp1 = ~consensus;
			tmp1 = tmp1 & (tmp1 << 1) & (~0UL/3UL*2UL);

			score += __builtin_popcountl(tmp1);
		}

		delete [] result1.consensus;
		delete [] result2.consensus;

		return score;
	}
}

#endif /* GAPPARSIMONY_H_ */
