#ifndef GRAPHALIGN_H_
#define GRAPHALIGN_H_

#include "Graph.h"

template <class ALPHABET>
struct AlignmentResult {
	dp_score_t score;
	index_t n_tr_indels;
	std::vector<index_t> mapping1;
	std::vector<index_t> mapping2;
};

template <class ALPHABET>
struct AncestralResult {
	Graph<ALPHABET> graph;
	std::vector<index_t> mapping1;
	std::vector<index_t> mapping2;
	std::vector<bool> is_matched;
};

template <class ALPHABET>
AlignmentResult<ALPHABET> alignGraphs(const Graph<ALPHABET> &g1, const Graph<ALPHABET> &g2, const Model<ALPHABET> &model);

template <class ALPHABET>
AncestralResult<ALPHABET> mergeGraphs(const Graph<ALPHABET> &g1, const Graph<ALPHABET> &g2, const std::vector<index_t> &mapping1_in, const std::vector<index_t> &mapping2_in, const Model<ALPHABET> &model1, const Model<ALPHABET> &model2);

/*
 *  ___                 _                           _        _   _
 * |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
 *  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
 *  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
 * |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
 *               |_|
 */

#include "main.h"
#include "CleanedGraph.h"
#include "debug.h"
#include "Alphabet.h"
#include "ls_log.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

#ifdef USE_LS_LOG
   #define LOG(x) (std::log(x)/std::log(2))
#else
   #define LOG(x) (std::log(x))
#endif

#define EPSILON	1e-6

typedef Eigen::Matrix<dp_score_t, Eigen::Dynamic, Eigen::Dynamic> DynProgMatrix;

template <class ALPHABET>
static double averageAlignmentLengthRec(const Graph<ALPHABET> &g, index_t current, double *cache) {
	if (cache[current] == -1.0) {
		double sum = 0;
		index_t paths = 0;

		for (typename Graph<ALPHABET>::PredIterator xit = g.getPreds(current,INFINITY,INFINITY); xit; ++xit) {
			index_t xp = *xit;
			if (xit.value() == 0.0) {
				double res = averageAlignmentLengthRec<ALPHABET>(g,xp,cache);
				if(res >= 0.0) {
					sum += res + 1.0;
					++paths;
				}
			}
		}

		if (paths > 0) {
			cache[current] = sum/paths;
		} else {
			cache[current] = -2.0;
		}
	}
	return cache[current];
}

template <class ALPHABET>
static double averageAlignmentLength(const Graph<ALPHABET> &g) {
	if (g.size() == 0) return 0;

	double *cache = new double[g.size()];
	for(index_t i=1; i < g.size(); ++i) {
		cache[i] = -1.0;
	}
	cache[0] = 0;

	double ret = averageAlignmentLengthRec<ALPHABET>(g,g.size()-1,cache);

	delete [] cache;
	return ret;
}

template <class ALPHABET>
struct DynProgScores {
	DynProgScores(const Graph<ALPHABET> &g1, const Graph<ALPHABET> &g2, const Model<ALPHABET> &model) {
		const double l1 = averageAlignmentLength(g1);
		const double l2 = averageAlignmentLength(g2);
		const double exp_length = std::max(l1, l2) * std::exp(model.distance * cmdlineopts.indel_rate * (model.epsilon / (1.0 - model.epsilon) + 1.0));
		const double nu = 2.0 / (2.0 + l1 + l2);

		double ttau = 1.0 / (1.0 + exp_length);
		if (model.epsilon + ttau >= 1.0) {
			ttau = (1.0 - model.epsilon)/2.0;
		}
		const double tau = ttau;

		gap_init = LOG(model.delta * (1.0 - model.epsilon - tau) / (1.0 - nu));
		gap_extend = LOG(model.epsilon / (1.0 - nu));
		match_init = LOG((1.0 - 2.0 * model.delta) * (1.0 - tau) / (1.0 - nu) / (1.0 - nu));
		end_skip = LOG(tau);
		if(cmdlineopts.end_indel_prob >= 0 && cmdlineopts.end_indel_prob <= 1) {
			end_match = LOG(tau  * (1.0 - cmdlineopts.end_indel_prob) / (1.0 - 2.0 * model.delta) / (1.0 - tau)); // XXX incorrect  if one sequence has zero length
			end_gap = LOG(tau * cmdlineopts.end_indel_prob / 2.0 / (1.0 - model.epsilon - tau) / model.delta);
			start_gap = LOG(cmdlineopts.end_indel_prob / 2.0 * (1.0 - model.epsilon - tau) / (1.0 - cmdlineopts.end_indel_prob) / (1.0 - nu));
			start_init = LOG((1.0 - tau) * (1.0 - cmdlineopts.end_indel_prob));
		} else {
			end_match = LOG(tau / (1.0 - tau));
			end_gap = LOG(tau / (1.0 - model.epsilon - tau));
			start_gap = LOG(model.delta * (1.0 - model.epsilon - tau) / (1.0 - nu));
			start_init = LOG(1.0 - tau);
		}

		const double repeat_prob = 1.0 - std::exp(-model.distance*cmdlineopts.repeat_rate);
		repeat_init = -LOG(std::min<double>(1,repeat_prob/(1-repeat_prob) * (1 - cmdlineopts.repeatext_prob)));
		repeat_ext = -LOG(std::min<double>(1,std::max<double>(0,cmdlineopts.repeatext_prob)));
	}

	dp_score_t gap_init;
	dp_score_t gap_extend;
	dp_score_t match_init;
	dp_score_t end_match;
	dp_score_t end_gap;
	dp_score_t end_skip;
	dp_score_t start_gap;
	dp_score_t start_init;
	dp_score_t repeat_init;
	dp_score_t repeat_ext;
};

template <class ALPHABET>
static DynProgMatrix precomputeScores(const Graph<ALPHABET> &g1, const Graph<ALPHABET> &g2, const Model<ALPHABET> &model, const DynProgScores<ALPHABET> &scores) {
	Eigen::Matrix<dp_score_t,ALPHABET::DIM,Eigen::Dynamic> g1s = g1.getSites().template cast<dp_score_t>();
	Eigen::Matrix<dp_score_t,ALPHABET::DIM,Eigen::Dynamic> g2s = g2.getSites().template cast<dp_score_t>();
	Eigen::Matrix<dp_score_t,ALPHABET::DIM,ALPHABET::DIM> M = model.M.template cast<dp_score_t>();
	Eigen::Matrix<dp_score_t,ALPHABET::DIM,1> pi = model.pi.template cast<dp_score_t>();

	DynProgMatrix S;

	S = (g1s.transpose() * (M.transpose() * g2s)).array() / ((g1s.transpose() * pi) * (pi.transpose() * g2s)).array();

#ifndef USE_LS_LOG
	S = S.array().log() + scores.match_init;
#else
	ls_log_add(S.data(),scores.match_init,S.rows()*S.cols());
#endif

	return S;
}

template <class ALPHABET>
void markAlternativePath(index_t start, index_t end, const Graph<ALPHABET> &g, std::vector<index_t> &mapping,  std::vector<index_t> &other_mapping) {
	index_t len = end - start + 1;
	std::vector<dp_score_t> score(len,-INFINITY);
	std::vector<index_t> prev(len,(index_t)-1);

	score[0] = 0;
	for (index_t i=1; i < len; ++i) {
		index_t real_ix = i + start;
		for (typename Graph<ALPHABET>::PredIterator it = g.getPreds(real_ix,INFINITY,INFINITY); it; ++it) {
			if (*it >= start && *it <= end) {
				index_t i2 = *it - start;
				if (score[i] <= score[i2] - it.value()) {
					score[i] = score[i2] - it.value();
					prev[i] = i2;
				}
			}
		}
	}

	if (score[len-1] > -INFINITY) {
		index_t i  = prev[len-1];
		while (i != 0) {
			mapping.push_back(i + start);
			other_mapping.push_back(-1);
			i = prev[i];
		}
	} else {
#if DEBUG
		g.printDot();
		error("no alternative path found");
#endif
	}
}

template <class ALPHABET>
AlignmentResult<ALPHABET> alignGraphs(const Graph<ALPHABET> &g1, const Graph<ALPHABET> &g2, const Model<ALPHABET> &model) {
	const DynProgScores<ALPHABET> s(g1,g2,model);
	const DynProgMatrix S = precomputeScores(g1,g2,model,s);

	const dp_score_t minfty = (dp_score_t) -INFINITY;
	DynProgMatrix M = DynProgMatrix::Constant((int) g1.size(), (int) g2.size(), minfty);
	DynProgMatrix X = DynProgMatrix::Constant((int) g1.size(), (int) g2.size(), minfty);
	DynProgMatrix Y = DynProgMatrix::Constant((int) g1.size(), (int) g2.size(), minfty);
	DynProgMatrix W = DynProgMatrix::Constant((int) g1.size(), (int) g2.size(), minfty);

	/* init horizontal and vertical start scores depending on distance to start node */
	W(0, 0) = s.start_init; // XXX everything starting from start node needs to be multiplied by 1/(1-tau) (?)

	for (index_t y = 1; y < g1.size() - 1; ++y) {
		dp_score_t Sy = minfty;
		for (typename Graph<ALPHABET>::PredIterator yit = g1.getPreds(y,s.repeat_init,s.repeat_ext); yit; ++yit) {
			index_t yp = *yit;
			Sy = std::max(Sy, std::max(Y(yp, 0) + s.gap_extend, W(yp, 0)
					+ s.start_gap) - yit.value());
		}
		Y(y, 0) = Sy;
		W(y, 0) = Y(y, 0);
	}

	for (index_t x = 1; x < g2.size() - 1; ++x) {
		dp_score_t Sx = minfty;
		for (typename Graph<ALPHABET>::PredIterator xit = g2.getPreds(x,s.repeat_init,s.repeat_ext); xit; ++xit) {
			index_t xp = *xit;
			Sx = std::max(Sx, std::max(X(0, xp) + s.gap_extend, W(0, xp)
					+ s.start_gap) - xit.value());
		}
		X(0, x) = Sx;
		W(0, x) = X(0, x);
	}

	/* fill scoring matrix */

	for (index_t y = 1; y < g1.size() - 1; ++y) {
		for (index_t x = 1; x < g2.size() - 1; ++x) {
			for (typename Graph<ALPHABET>::PredIterator yit = g1.getPreds(y,s.repeat_init,s.repeat_ext); yit; ++yit) {
				for (typename Graph<ALPHABET>::PredIterator xit = g2.getPreds(x,s.repeat_init,s.repeat_ext); xit; ++xit) {
					index_t yp = *yit;
					index_t xp = *xit;

					dp_score_t Sm = W(yp, xp) + S(y,x)
							- yit.value() - xit.value();
					dp_score_t Sx = std::max(X(y, xp) + s.gap_extend, W(y, xp)
							+ s.gap_init) - xit.value();
					dp_score_t Sy = std::max(Y(yp, x) + s.gap_extend, W(yp, x)
							+ s.gap_init) - yit.value();
					dp_score_t Sw = std::max(Sm, std::max(Sx, Sy));

					M(y, x) = std::max((dp_score_t)M(y, x), Sm);
					X(y, x) = std::max((dp_score_t)X(y, x), Sx);
					Y(y, x) = std::max((dp_score_t)Y(y, x), Sy);
					W(y, x) = std::max((dp_score_t)W(y, x), Sw);
				}
			}
		}
	}

	/* match of end-node */

	dp_score_t Wend = minfty;

	for (typename Graph<ALPHABET>::PredIterator yit = g1.getPreds(g1.size() - 1,s.repeat_init,s.repeat_ext); yit; ++yit) {
		for (typename Graph<ALPHABET>::PredIterator xit = g2.getPreds(g2.size() - 1,s.repeat_init,s.repeat_ext); xit; ++xit) {
			index_t yp = *yit;
			index_t xp = *xit;

			if(xp == 0 && yp == 0) {
				Wend = std::max(s.end_skip - yit.value() - xit.value(), Wend);
			} else {
				//Wend = std::max(W(yp, xp) - yit.value() - xit.value(), Wend);
				Wend = std::max(X(yp, xp) + s.end_gap - yit.value() - xit.value(), Wend);
				Wend = std::max(Y(yp, xp) + s.end_gap - yit.value() - xit.value(), Wend);
				Wend = std::max(M(yp, xp) + s.end_match - yit.value() - xit.value(), Wend);
			}
		}
	}


	/* backtracking */

	AlignmentResult<ALPHABET> result;
	result.score = Wend;
	result.n_tr_indels = 0;

	enum State {
		State_m, State_x, State_y
	} current_state = State_m;
	dp_score_t current_score = minfty;
	index_t y = g1.size() - 1, x = g2.size() - 1;
	std::vector<index_t> mapping1;
	std::vector<index_t> mapping2;

	mapping1.reserve(g1.size() + g2.size());
	mapping2.reserve(g1.size() + g2.size());

#define mapping(y,x) 	{mapping1.push_back(y); mapping2.push_back(x);}

	mapping(g1.size()-1,g2.size()-1);

	// transitions to the end state
	bool tr_indel_x = false;
	bool tr_indel_y  = false;
	dp_score_t best_match = INFINITY;
	for (typename Graph<ALPHABET>::PredIterator yit = g1.getPreds(g1.size() - 1,s.repeat_init,s.repeat_ext); yit; ++yit) {
		for (typename Graph<ALPHABET>::PredIterator xit = g2.getPreds(g2.size() - 1,s.repeat_init,s.repeat_ext); xit; ++xit) {
			if (best_match > std::abs(Wend - (M(*yit,*xit) + s.end_match - yit.value() - xit.value()))) {
				best_match = std::abs(Wend - (M(*yit,*xit) + s.end_match - yit.value() - xit.value()));
				tr_indel_x = xit.isRepeat();
				tr_indel_y = yit.isRepeat();
				current_score = M(*yit, *xit);
				current_state = State_m;
				y = *yit;
				x = *xit;
			}
			if (best_match > std::abs(Wend - (Y(*yit,*xit) + s.end_gap - yit.value() - xit.value()))) {
				best_match = std::abs(Wend - (Y(*yit,*xit) + s.end_gap - yit.value() - xit.value()));
				tr_indel_x = xit.isRepeat();
				tr_indel_y = yit.isRepeat();
				current_score = Y(*yit, *xit);
				current_state = State_y;
				y = *yit;
				x = *xit;
			}
			if (best_match > std::abs(Wend - (X(*yit,*xit) + s.end_gap - yit.value() - xit.value()))) {
				best_match = std::abs(Wend - (X(*yit,*xit) + s.end_gap - yit.value() - xit.value()));
				tr_indel_x = xit.isRepeat();
				tr_indel_y = yit.isRepeat();
				current_score = X(*yit, *xit);
				current_state = State_x;
				y = *yit;
				x = *xit;
			}
			if(*xit == 0 && *yit == 0 && best_match > std::abs(Wend - (s.end_skip - yit.value() - xit.value()))) {
				best_match = std::abs(Wend - (s.end_skip - yit.value() - xit.value()));
				tr_indel_x = xit.isRepeat();
				tr_indel_y = yit.isRepeat();
				y = *yit;
				x = *xit;
			}
		}
	}
	result.n_tr_indels += tr_indel_x + tr_indel_y;

	if(tr_indel_y) {
		markAlternativePath(y,g1.size()-1,g1,mapping1,mapping2);
	}
	if(tr_indel_x) {
		markAlternativePath(x,g2.size()-1,g2,mapping2,mapping1);
	}

#ifdef DEBUG
	if(!(best_match/Wend < EPSILON)) error("backtracking failed (smallest relative difference %g)",best_match/Wend,best_match,Wend);
#endif

	if(x != 0 || y != 0) {
		if(current_state == State_m) {
			assert(x != 0 && y != 0);
			mapping(y,x);
		} else if(current_state == State_x) {
			mapping(-1,x);
		} else if(current_state == State_y) {
			mapping(y,-1);
		}
	}

	// not-end transitions
	dp_score_t next_score = INFINITY;
	State next_state = State_m;
	index_t next_x = -1, next_y = -1;
	while (x != 0 || y != 0) {
		best_match = INFINITY;

		// gap extensions and openings
		if (current_state == State_y) {
			for (typename Graph<ALPHABET>::PredIterator yit = g1.getPreds(y,s.repeat_init,s.repeat_ext); yit; ++yit) {
				index_t yp = *yit;

				if (best_match > std::abs(current_score - (Y(yp, x) + s.gap_extend - yit.value()))) {
					best_match = std::abs(current_score - (Y(yp, x) + s.gap_extend - yit.value()));
					tr_indel_x = false;
					tr_indel_y = yit.isRepeat();
					next_x = x;
					next_y = yp;
					next_score = Y(next_y, next_x);
					next_state = State_y;
				}

				if (best_match > std::abs(current_score - (W(yp, x) + s.gap_init - yit.value()))) {
					best_match = std::abs(current_score - (W(yp, x) + s.gap_init - yit.value()));
					tr_indel_x = false;
					tr_indel_y = yit.isRepeat();
					next_x = x;
					next_y = yp;

					if (next_x != 0 || next_y != 0) {
						if (W(next_y,next_x) == M(next_y,next_x)) {
							next_score = M(next_y, next_x);
							next_state = State_m;
						} else if (W(next_y,next_x) == Y(next_y,next_x)) {
							next_score = Y(next_y, next_x);
							next_state = State_y;
						} else if (W(next_y,next_x) == X(next_y,next_x)) {
							next_score = X(next_y, next_x);
							next_state = State_x;
						} else {
							error("backtracking failed");
						}
					}
				}
			}
		}

		if (current_state == State_x) {
			for (typename Graph<ALPHABET>::PredIterator xit = g2.getPreds(x,s.repeat_init,s.repeat_ext); xit; ++xit) {
				index_t xp = *xit;

				if (best_match > std::abs(current_score - (X(y, xp) + s.gap_extend - xit.value()))) {
					best_match = std::abs(current_score - (X(y, xp) + s.gap_extend - xit.value()));
					tr_indel_x = xit.isRepeat();
					tr_indel_y = false;
					next_x = xp;
					next_y = y;
					next_score = X(next_y, next_x);
					next_state = State_x;
				}

				if (best_match > std::abs(current_score - (W(y, xp) + s.gap_init - xit.value()))) {
					best_match = std::abs(current_score - (W(y, xp) + s.gap_init - xit.value()));
					tr_indel_x = xit.isRepeat();
					tr_indel_y = false;
					next_x = xp;
					next_y = y;

					if (next_x != 0 || next_y != 0) {
						if (W(next_y,next_x) == M(next_y,next_x)) {
							next_score = M(next_y, next_x);
							next_state = State_m;
						} else if (W(next_y,next_x) == Y(next_y,next_x)) {
							next_score = Y(next_y, next_x);
							next_state = State_y;
						} else if (W(next_y,next_x) == X(next_y,next_x)) {
							next_score = X(next_y, next_x);
							next_state = State_x;
						} else {
							error("backtracking failed");
						}
					}
				}
			}
		}

		// transitions from match state
		if (current_state == State_m) {
			for (typename Graph<ALPHABET>::PredIterator yit = g1.getPreds(y,s.repeat_init,s.repeat_ext); yit; ++yit) {
				for (typename Graph<ALPHABET>::PredIterator xit = g2.getPreds(x,s.repeat_init,s.repeat_ext); xit; ++xit) {
					index_t yp = *yit;
					index_t xp = *xit;

					if (best_match > std::abs(current_score - (W(yp,xp) + S(y,x)
									- yit.value() - xit.value()))) {
						best_match = std::abs(current_score - (W(yp,xp) + S(y,x) - yit.value() - xit.value()));
						tr_indel_x = xit.isRepeat();
						tr_indel_y = yit.isRepeat();
						next_y = yp;
						next_x = xp;

						if (next_x != 0 || next_y != 0) {
							if (W(next_y,next_x) == M(next_y,next_x)) {
								next_score = M(next_y, next_x);
								next_state = State_m;
							} else if (W(next_y,next_x) == Y(next_y,next_x)) {
								next_score = Y(next_y, next_x);
								next_state = State_y;
							} else if (W(next_y,next_x) == X(next_y,next_x)) {
								next_score = X(next_y, next_x);
								next_state = State_x;
							} else {
								error("backtracking failed");
							}
						}
					}
				}
			}
		}
		result.n_tr_indels += tr_indel_x + tr_indel_y;

		if(tr_indel_y) {
			markAlternativePath(next_y,y,g1,mapping1,mapping2);
		}
		if(tr_indel_x) {
			markAlternativePath(next_x,x,g2,mapping2,mapping1);
		}

#ifdef DEBUG
		if(!(best_match/current_score < EPSILON)) error("backtracking failed (smallest relative difference %g)",best_match/current_score,best_match,current_score);
#endif

		x = next_x;
		y = next_y;
		current_state = next_state;
		current_score = next_score;

		if (x != 0 || y != 0) {
			if(current_state == State_m) {
				assert(x != 0 && y != 0);
				mapping(y,x);
			} else if(current_state == State_x) {
				mapping(-1,x);
			} else if(current_state == State_y) {
				mapping(y,-1);
			}
		}
	}
	mapping(0,0);

	std::reverse(mapping1.begin(), mapping1.end());
	std::reverse(mapping2.begin(), mapping2.end());

	/* free DynProg matrices */
	W.resize(0,0);
	Y.resize(0,0);
	X.resize(0,0);
	M.resize(0,0);
	const_cast<DynProgMatrix*>(&S)->resize(0,0);

	result.mapping1 = mapping1;
	result.mapping2 = mapping2;

	return result;
}

#undef EPSILON
#undef mapping

template <class M, class V>
inline void updateEdge(M &map, index_t from, index_t to, V cost) {
	std::pair<index_t,index_t> i(to,from);
	typename M::iterator it = map.find(i);
	if(it != map.end()) {
		it->second = std::min(it->second,cost);
	} else {
		map[i] = cost;
	}
}

template <class ALPHABET>
AncestralResult<ALPHABET> mergeGraphs(const Graph<ALPHABET> &g1, const Graph<ALPHABET> &g2, const std::vector<index_t> &mapping1, const std::vector<index_t> &mapping2, const Model<ALPHABET> &model1, const Model<ALPHABET> &model2, double support1, double support2)
{
	std::vector<typename Model<ALPHABET>::Freqs,Eigen::aligned_allocator<typename Model<ALPHABET>::Freqs> > nodes;
	nodes.reserve(g1.size() + g2.size());

	std::map<std::pair<index_t,index_t>, dp_score_t> edges;
	std::map<std::pair<index_t,index_t>, index_t> repeats;

	assert(mapping1.size() == mapping2.size());

	AncestralResult<ALPHABET> result;
	result.mapping1.reserve(g1.size() + g2.size());
	result.mapping2.reserve(g1.size() + g2.size());
	result.is_matched.reserve(g1.size() + g2.size());


	/* unify graphs */

	for (index_t i1 = 0, i2 = 0, j = 0; j < mapping1.size(); ++j) {
		index_t k1 = mapping1[j];
		index_t k2 = mapping2[j];

		if (k1 != (index_t) -1) {
			if (i1 != k1) {
				assert(i1 <= k1);
				for (; i1 != k1; ++i1) {
					typename Model<ALPHABET>::Freqs p = model1.P * g1[i1];
					nodes.push_back(p.norm() == 0 ? p : p.normalized());
					result.mapping1.push_back(i1);
					result.mapping2.push_back(-1);
					result.is_matched.push_back(false);
				}
			}
			++i1;
		}

		if (k2 != (index_t) -1) {
			if (i2 != k2) {
				assert(i2 <= k2);
				for (; i2 != k2; ++i2) {
					typename Model<ALPHABET>::Freqs p = model1.P * g2[i2];
					nodes.push_back(p.norm() == 0 ? p : p.normalized());
					result.mapping1.push_back(-1);
					result.mapping2.push_back(i2);
					result.is_matched.push_back(false);
				}
			}
			++i2;
		}

		if(k1 != (index_t)-1 && k2 != (index_t)-1) {
			typename Model<ALPHABET>::Freqs p = (model1.P * g1[k1]).array() * (model2.P * g2[k2]).array();
			nodes.push_back(p.norm() == 0 ? p : p.normalized());
			result.mapping1.push_back(k1);
			result.mapping2.push_back(k2);
		} else if(k1 != (index_t)-1) {
			typename Model<ALPHABET>::Freqs p = model1.P * g1[k1];
			nodes.push_back(p.norm() == 0 ? p : p.normalized());
			result.mapping1.push_back(k1);
			result.mapping2.push_back(-1);
		} else if(k2 != (index_t)-1) {
			typename Model<ALPHABET>::Freqs p = model2.P * g2[k2];
			nodes.push_back(p.norm() == 0 ? p : p.normalized());
			result.mapping1.push_back(-1);
			result.mapping2.push_back(k2);
		} else {
			assert(false && "error in mapping");
		}
		result.is_matched.push_back(true);
	}

	assert(nodes.size() == result.mapping1.size() && nodes.size() == result.mapping2.size() && nodes.size() == result.is_matched.size());

	/* homologous path and allow skipping newly inserted gaps (insertion) */

	index_t last_xy = 0;
	index_t last_x = 0;
	index_t last_y = 0;
	index_t last_mapped = 0;

	for (index_t i = 1; i < nodes.size(); ++i) {
		if(!result.is_matched[i]) continue;

		updateEdge(edges,last_mapped,i,(dp_score_t)0);
		last_mapped = i;

		if(result.mapping1[i] != (index_t) -1 && result.mapping2[i] != (index_t) -1) {
			if(last_xy != i-1) {
				updateEdge(edges,last_xy,i,(dp_score_t)0);
			}
			last_xy = i;
		}

		if(result.mapping1[i] != (index_t) -1) {
			if(last_y != i-1) {
				updateEdge(edges,last_y,i,(dp_score_t)0);
			}
			last_y = i;
		}

		if(result.mapping2[i] != (index_t) -1) {
			if(last_x != i-1) {
				updateEdge(edges,last_x,i,(dp_score_t)0);
			}
			last_x = i;
		}
	}

	/* construct inverse mapping for efficiency */

	std::vector<index_t> inv_mapping1(g1.size(), 0);
	std::vector<index_t> inv_mapping2(g2.size(), 0);

	for (index_t i = 0; i < result.mapping1.size(); ++i) {
		index_t m = result.mapping1[i];
		if (m != (index_t) -1)
			inv_mapping1[m] = i;
	}
	for (index_t i = 0; i < result.mapping2.size(); ++i) {
		index_t m = result.mapping2[i];
		if (m != (index_t) -1)
			inv_mapping2[m] = i;
	}

	/* compute penalties for unused edges */

	double unused_prob1 = cmdlineopts.altsplice_prob + (1.0-cmdlineopts.altsplice_prob)*(1.0-support1);
	dp_score_t unused_penalty1 = -LOG(unused_prob1);

	double unused_prob2 = cmdlineopts.altsplice_prob + (1.0-cmdlineopts.altsplice_prob)*(1.0-support2);
	dp_score_t unused_penalty2 = -LOG(unused_prob2);

	/* add missing edges */

	for (index_t to = 0; to < g1.size(); ++to) {
		for (typename Graph<ALPHABET>::PredIterator from = g1.getPreds(to,0,0); from; ++from) {
			index_t y = inv_mapping1[*from];
			index_t x = inv_mapping1[to];

			// add additional penalty if this is an edge leading in or out of the matched subgraph
			if(!from.isRepeat()) {
				if(result.is_matched[*from] && result.is_matched[to]) {
					updateEdge(edges, y, x, from.value() + unused_penalty1);
				} else if(result.is_matched[*from] || result.is_matched[to]) {
					updateEdge(edges, y, x, from.value() + unused_penalty1/2);
				} else {
					updateEdge(edges, y, x, from.value());
				}
			} else {
				updateEdge(repeats, y, x, from.repeatUnits());
			}
		}
	}
	for (index_t to = 0; to < g2.size(); ++to) {
		for (typename Graph<ALPHABET>::PredIterator from = g2.getPreds(to,0,0); from; ++from) {
			index_t y = inv_mapping2[*from];
			index_t x = inv_mapping2[to];

			// add additional penalty if this is an edge leading in or out of the matched subgraph
			if(!from.isRepeat()) {
				if(result.is_matched[*from] && result.is_matched[to]) {
					updateEdge(edges, y, x, from.value() + unused_penalty2);
				} else if(result.is_matched[*from] || result.is_matched[to]) {
					updateEdge(edges, y, x, from.value() + unused_penalty2/2);
				} else {
					updateEdge(edges, y, x, from.value());
				}
			} else {
				updateEdge(repeats, y, x, from.repeatUnits());
			}
		}
	}

	result.graph = Graph<ALPHABET>(nodes,edges,repeats);

	return result;
}

template <class ALPHABET>
AncestralResult<ALPHABET> mergeGraphsIncremental(const Graph<ALPHABET> &anc_graph, const Graph<ALPHABET> &graph, const std::vector<index_t> &anc_mapping, const std::vector<index_t> &mapping, const Model<ALPHABET> &model)
{
	std::vector<typename Model<ALPHABET>::Freqs,Eigen::aligned_allocator<typename Model<ALPHABET>::Freqs> > nodes;
	nodes.reserve(anc_graph.size() + graph.size());

	std::map<std::pair<index_t,index_t>, dp_score_t> edges;
	std::map<std::pair<index_t,index_t>, index_t> repeats;

	assert(anc_mapping.size() == mapping.size());

	AncestralResult<ALPHABET> result;
	result.mapping1.reserve(anc_graph.size() + graph.size());
	result.mapping2.reserve(anc_graph.size() + graph.size());
	result.is_matched.reserve(anc_graph.size() + graph.size());


	/* unify graphs */

	for (index_t i1 = 0, i2 = 0, j = 0; j < anc_mapping.size(); ++j) {
		index_t k1 = anc_mapping[j];
		index_t k2 = mapping[j];

		if (k1 != (index_t) -1) {
			if (i1 != k1) {
				assert(i1 <= k1);
				for (; i1 != k1; ++i1) {
					typename Model<ALPHABET>::Freqs p = anc_graph[i1];
					nodes.push_back(p.norm() == 0 ? p : p.normalized());
					result.mapping1.push_back(i1);
					result.mapping2.push_back(-1);
					result.is_matched.push_back(false);
				}
			}
			++i1;
		}

		if (k2 != (index_t) -1) {
			if (i2 != k2) {
				assert(i2 <= k2);
				for (; i2 != k2; ++i2) {
					typename Model<ALPHABET>::Freqs p = model.P * graph[i2];
					nodes.push_back(p.norm() == 0 ? p : p.normalized());
					result.mapping1.push_back(-1);
					result.mapping2.push_back(i2);
					result.is_matched.push_back(false);
				}
			}
			++i2;
		}

		if(k1 != (index_t)-1 && k2 != (index_t)-1) {
			typename Model<ALPHABET>::Freqs p = anc_graph[k1].array() * (model.P * graph[k2]).array();
			nodes.push_back(p.norm() == 0 ? p : p.normalized());
			result.mapping1.push_back(k1);
			result.mapping2.push_back(k2);
		} else if(k1 != (index_t)-1) {
			typename Model<ALPHABET>::Freqs p = anc_graph[k1];
			nodes.push_back(p.norm() == 0 ? p : p.normalized());
			result.mapping1.push_back(k1);
			result.mapping2.push_back(-1);
		} else if(k2 != (index_t)-1) {
			typename Model<ALPHABET>::Freqs p = model.P * graph[k2];
			nodes.push_back(p.norm() == 0 ? p : p.normalized());
			result.mapping1.push_back(-1);
			result.mapping2.push_back(k2);
		} else {
			assert(false && "error in mapping");
		}
		result.is_matched.push_back(true);
	}

	assert(nodes.size() == result.mapping1.size() && nodes.size() == result.mapping2.size() && nodes.size() == result.is_matched.size());

	/* homologous path and allow skipping newly inserted gaps (insertion) */

	index_t last_xy = 0;
	index_t last_x = 0;
	index_t last_y = 0;
	index_t last_mapped = 0;

	for (index_t i = 1; i < nodes.size(); ++i) {
		if(!result.is_matched[i]) continue;

		updateEdge(edges,last_mapped,i,(dp_score_t)0);
		last_mapped = i;

		if(result.mapping1[i] != (index_t) -1 && result.mapping2[i] != (index_t) -1) {
			if(last_xy != i-1) {
				updateEdge(edges,last_xy,i,(dp_score_t)0);
			}
			last_xy = i;
		}

		if(result.mapping1[i] != (index_t) -1) {
			if(last_y != i-1) {
				updateEdge(edges,last_y,i,(dp_score_t)0);
			}
			last_y = i;
		}

		if(result.mapping2[i] != (index_t) -1) {
			if(last_x != i-1) {
				updateEdge(edges,last_x,i,(dp_score_t)0);
			}
			last_x = i;
		}
	}

	/* construct inverse mapping for efficiency */

	std::vector<index_t> inv_mapping1(anc_graph.size(), 0);
	std::vector<index_t> inv_mapping2(graph.size(), 0);

	for (index_t i = 0; i < result.mapping1.size(); ++i) {
		index_t m = result.mapping1[i];
		if (m != (index_t) -1)
			inv_mapping1[m] = i;
	}
	for (index_t i = 0; i < result.mapping2.size(); ++i) {
		index_t m = result.mapping2[i];
		if (m != (index_t) -1)
			inv_mapping2[m] = i;
	}

	/* add missing edges */

	for (index_t to = 0; to < anc_graph.size(); ++to) {
		for (typename Graph<ALPHABET>::PredIterator from = anc_graph.getPreds(to,0,0); from; ++from) {
			index_t y = inv_mapping1[*from];
			index_t x = inv_mapping1[to];

			if(!from.isRepeat()) {
				updateEdge(edges, y, x, from.value());
			} else {
				updateEdge(repeats, y, x, from.repeatUnits());
			}
		}
	}
	for (index_t to = 0; to < graph.size(); ++to) {
		for (typename Graph<ALPHABET>::PredIterator from = graph.getPreds(to,0,0); from; ++from) {
			index_t y = inv_mapping2[*from];
			index_t x = inv_mapping2[to];

			if(!from.isRepeat()) {
				updateEdge(edges, y, x, from.value());
			} else {
				updateEdge(repeats, y, x, from.repeatUnits());
			}
		}
	}

	result.graph = Graph<ALPHABET>(nodes,edges,repeats);

	return result;
}

#undef LOG

#endif /* GRAPHALIGN_H_ */
