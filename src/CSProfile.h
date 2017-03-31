#ifndef CSPROFILE_H_
#define CSPROFILE_H_

#include <string>
#include <vector>
#include <Eigen/Core>

#include "Model.h"
#include "debug.h"
#include "Alphabet.h"

class CSProfile {
public:
	CSProfile(const std::string &filename);
	~CSProfile();

	class csprofile_exception: public swps3_exception {
	public:
		csprofile_exception(std::string str) : swps3_exception(str) {};
	};
	Model<AA>::Profile createProfile(const sequence_t<AA> &seq, const Model<AA> &model) const;

protected:
	int nprof;
	int ncols;
	typedef double cs_score_t;

	std::vector<Eigen::Matrix<cs_score_t,Eigen::Dynamic,(AA::DIM)> > profiles;
	std::vector<Eigen::Matrix<cs_score_t,Eigen::Dynamic,(AA::DIM)+1> > lprofiles;
	std::vector<cs_score_t> priors;
	Eigen::Matrix<cs_score_t,Eigen::Dynamic,1> weights;
};

#endif /* CSPROFILE_H_ */
