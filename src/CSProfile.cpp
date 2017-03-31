#include "CSProfile.h"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>

#include "main.h"
#include "Alphabet.h"
#include "Model.h"

/* This file implements the context-specific profile algorithm from:
 *
 *   Biegert, A. and Soding, J. (2009) Sequence context-specific profiles for homology searching.
 *   Proc Natl Acad Sci USA, 106 (10), 3770-3775
 */

static const double w_center = .26236426446749105203; /* log(1.3) */
static const double beta = -.10536051565782630122; /* log(.9) */
static const double log_2 = .69314718055994530941;
static const double rlog_2 = 1.0/.69314718055994530941;
static const double lambda = .34657359027997265470;
static const double rlambda = 1.0/.34657359027997265470;

#define MATRIX_DIM (AA::DIM)
typedef Model<AA>::Freqs AAfreq;
typedef Model<AA>::Subst AAsubst;

CSProfile::CSProfile(const std::string &filename) {
	this->nprof = -1;
	this->ncols = -1;

	std::ifstream file(filename.c_str());
	std::string line;

	if(!std::getline(file,line) || line.find("ProfileLibrary") != 0)
		throw csprofile_exception("error opening profile library");

	// parse header
	while(std::getline(file,line)) {
		if(line.at(0) == '#' || line.length() == 0) {
			continue;
		} else if(line.find("NPROF") == 0) {
			std::string ignore;
			std::istringstream is(line);
			is >> ignore >> this->nprof;

			if(!is || this->nprof <= 0) throw csprofile_exception("parse error: " + line);
		} else if(line.find("NCOLS") == 0) {
			std::string ignore;
			std::istringstream is(line);
			is >> ignore >> this->ncols;

			if(!is || this->ncols <= 0) throw csprofile_exception("parse error: " + line);
		} else if(line.find("ITERS") == 0) {
			// ignore
		} else if(line.find("LOG") == 0) {
			// ignore
		} else if(line.find("ContextProfile") == 0) {
			break;
		} else {
			throw csprofile_exception("parse error: " + line);
		}
	}

	if(this->nprof <= 0 || this->ncols <= 0) {
		throw csprofile_exception("missing information in header");
	}

	this->profiles.resize(this->nprof,Eigen::Matrix<cs_score_t,Eigen::Dynamic,Eigen::Dynamic>::Zero(this->ncols,MATRIX_DIM));
	this->lprofiles.resize(this->nprof,Eigen::Matrix<cs_score_t,Eigen::Dynamic,Eigen::Dynamic>::Zero(this->ncols,MATRIX_DIM+1));
	this->priors.resize(this->nprof,0);

	this->weights.resize(this->ncols);
	const int center = this->ncols / 2;
	for (int j = -center; j <= center; ++j) {
		this->weights(center + j) = std::exp(w_center + beta * std::abs(j));
	}

	// parse profiles
	do {
		if(line.at(0) == '#' || line.length() == 0) {
			continue;
		} else if(line.find("ContextProfile") == 0) {
			int index = -1;
			cs_score_t prior = -1;
			Eigen::Matrix<cs_score_t,Eigen::Dynamic,Eigen::Dynamic> profile = Eigen::Matrix<cs_score_t,Eigen::Dynamic,Eigen::Dynamic>::Zero(this->ncols,MATRIX_DIM);

			while(std::getline(file,line)) {
				if(line.at(0) == '#' || line.length() == 0 || line.find("ITERS") == 0) {
					continue;
				} else if(line.find("INDEX") == 0) {
					std::string ignore;
					std::istringstream is(line);
					is >> ignore >> index;

					if(!is || index < 0 || index >= this->nprof) throw csprofile_exception("parse error: " + line);
				} else if(line.find("PRIOR") == 0) {
					std::string ignore;
					std::istringstream is(line);
					is >> ignore >> prior;

					if(!is || prior <= 0) throw csprofile_exception("parse error: " + line);
				} else if(line.find("NCOLS") == 0) {
					std::string ignore;
					std::istringstream is(line);
					int pncols;

					is >> ignore >> pncols;

					if(!is || pncols != this->ncols) throw csprofile_exception("parse error: " + line);
				} else if(line.find("ALPH") == 0) {
					std::string ignore;
					std::istringstream is(line);

					int palph;

					is >> ignore >> palph;

					if(!is || palph != MATRIX_DIM) throw csprofile_exception("parse error: " + line);
				} else if(line.find("LOG") == 0) {
					// ignore
				} else if(std::isspace(line[0])) {
					std::istringstream scols(line);
					std::vector<int> cols(MATRIX_DIM);

					for(index_t i=0; i<MATRIX_DIM; ++i) {
						char s = 0;
						scols >> s;
						cols[i] = AA(s).value();
						if(!scols || cols[i] < 0 || cols[i] >= MATRIX_DIM) throw csprofile_exception("parse error in column names");
					}

					while(std::getline(file,line)) {
						if(line == "//") {
							goto finish_profile;
						} else {
							std::istringstream row(line);
							int col = 0;
							row >> col;

							if(!row || col <= 0 || col > this->ncols) throw csprofile_exception("parse error: invalid column number");

							for(index_t i=0; i<MATRIX_DIM; ++i) {
								row >> profile(col-1,cols[i]);

								if(!row || profile(col-1,cols[i]) < 0) throw csprofile_exception("parse error in profile");
							}
						}
					}
				} else {
					throw csprofile_exception("parse error: " + line);
				}
			}
			finish_profile:
			if(index < 0 || index >= this->nprof) throw csprofile_exception("parse error: invalid index");
			profile = (profile * (-log_2/1000.0)).array().exp();
			this->lprofiles[index].block(0,0,this->ncols,MATRIX_DIM) = profile.array().log();
			this->lprofiles[index].block(0,0,this->ncols,MATRIX_DIM) -=  profile.rowwise().sum().array().log().matrix().asDiagonal() * Eigen::Matrix<cs_score_t,Eigen::Dynamic,Eigen::Dynamic>::Ones(this->ncols,MATRIX_DIM);
			this->lprofiles[index].col(MATRIX_DIM).setZero();
			this->profiles[index] = this->lprofiles[index].block(0,0,this->ncols,MATRIX_DIM).array().exp();
			this->lprofiles[index] = this->weights.asDiagonal() * this->lprofiles[index];

			if(prior <= 0) throw csprofile_exception("parse error: invalid prior");
			this->priors[index] = std::log(prior);
		} else {
			throw csprofile_exception("parse error: " + line);
		}
	} while(std::getline(file,line));
}

CSProfile::~CSProfile() {
}

Model<AA>::Profile CSProfile::createProfile(const sequence_t<AA> &seq, const Model<AA> &model) const {
	Eigen::Matrix<cs_score_t,Eigen::Dynamic,MATRIX_DIM> profile = Eigen::Matrix<cs_score_t,Eigen::Dynamic,MATRIX_DIM>::Zero(seq.length() + 2, MATRIX_DIM);

	const cs_score_t tau = model.divergence / 0.8;
	const int center = this->ncols / 2;

	std::vector<index_t> tseq(seq.length());
	for(index_t i=0; i < seq.length(); ++i) {
		if(!seq[i].isValid()) {
			tseq[i] = MATRIX_DIM;
		} else {
			tseq[i] = seq[i].value();
		}
	}

	for (index_t k = 0; (int)k < this->nprof; ++k) {
		for (index_t i = 0; i < seq.length(); ++i) {
			/* compute pk = log(P(p_k,X_i)) */
			cs_score_t pk = this->priors[k];
			for (int j = -center; j <= center; ++j) {
				if ((int)i + j >= 0 && i + j < seq.length()) {
					int cj = tseq[i + j];
					pk += this->lprofiles[k](j + center, cj);
				}
			}

			profile.row(i + 1) += this->profiles[k].row(center) * std::exp(pk);
		}
	}

	for (index_t i = 0; i < seq.length(); ++i) {
		int c = tseq[i];

		if(profile.row(i + 1).sum() <= 0) {
			score_t s = 1.0/20;
			profile.row(i + 1) = (model.P * AAfreq::Constant(s)).cast<cs_score_t>().transpose();
		} else if (c < 0 || c > 19) {
			profile.row(i + 1) *= 1.0 / profile.row(i + 1).sum();
			profile.row(i + 1).array() *= ((1.0/20.0) * model.pi.array().inverse().matrix().transpose()).cast<cs_score_t>().array();
		} else {
			profile.row(i + 1) *= tau / profile.row(i + 1).sum();
			profile(i + 1, c) += 1.0 - tau;
			if(profile(i + 1, c) <= 0.0) {
				profile(i + 1, c) = 1e-3;
			}
			profile.row(i + 1).array() *= ((1.0/20.0) * model.pi.array().inverse().matrix().transpose()).cast<cs_score_t>().array();
		}
	}

	return profile.cast<score_t>().transpose();
}
