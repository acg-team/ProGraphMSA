#include <Eigen/Dense>
#include <bitset>

#define TOL 1e-6
#define MAX_ITER 100

template <class Derived, class OtherDerived>
Eigen::Matrix<double,Derived::ColsAtCompileTime,1> NNLS(const Eigen::MatrixBase<Derived> &Z, const Eigen::MatrixBase<OtherDerived> &x) {
	Eigen::Matrix<double,Derived::ColsAtCompileTime,1>  d(Z.cols());
	Eigen::Matrix<bool,Derived::ColsAtCompileTime,1>    P(Z.cols());
	Eigen::Matrix<double,Derived::ColsAtCompileTime,1>  w(Z.cols());
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Zp(Z.rows(),Z.cols());
	Eigen::Matrix<double,Eigen::Dynamic,1>              dp(Z.cols());
	Eigen::Matrix<double,Eigen::Dynamic,1>              sp(Z.cols());
	Eigen::Matrix<double,Eigen::Dynamic,1>              alpha(Z.cols());
	Eigen::Matrix<index_t,Eigen::Dynamic,1>             mapping(Z.cols());

	Zp = Z;
	d = Zp.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(x);
	if(d.minCoeff() >= 0) {
		return d;
	}
	
	P.setZero();
	d.setZero();

	w = Z.transpose()*(x - Z*d);
	w.array() *= 1.0 - P.template cast<double>().array();
	typename Eigen::Matrix<double,Derived::RowsAtCompileTime,1>::Index iw;
	index_t iiw = 0;
	index_t n_iter = 0;

	while(!P.all() && w.maxCoeff(&iw) > TOL) {
		P[iw] = true;

		if(n_iter++ > MAX_ITER) {
#ifdef DEBUG
			std::cerr << "NNLS max iterations reached!" << std::endl;
			std::cerr << "Z =" << std::endl;
			std::cerr << Z << std::endl;
			std::cerr << "x =" << std::endl;
			std::cerr << x.transpose() << std::endl;
			std::cerr << "d =" << std::endl;
			std::cerr << d.transpose() << std::endl;
#endif
			return d;
		}

		do {
			index_t P_size = P.count();

			Zp.resize(Z.rows(),P_size);
			sp.resize(P_size);
			dp.resize(P_size);
			mapping.resize(P_size);

			index_t k=0;
			for(index_t i=0; i < (index_t)Z.cols(); ++i) {
				if(P[i]) {
					Zp.col(k) = Z.col(i);
					dp[k] = d[i];
					mapping[k] = i;
					if(i == iw) {
						iiw = k;
					}
					++k;
				}
			}
			sp = Zp.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(x);

			if(sp.minCoeff() > 0) {
				for(index_t i=0; i < (index_t)sp.rows(); ++i) {
					d[mapping[i]] = sp[i];
				}
				w = Z.transpose()*(x - Z*d);
				w.array() *= 1.0 - P.template cast<double>().array();

				break;
			} else if(sp[iiw] <= 0) {
				if(sp[iiw] <= 0) {
					w[iw] = 0;
				}

				break;
			}

			alpha = dp.array()/(dp.array()-sp.array());
			//alpha.array() *= (sp.array() <= 0).template cast<double>();
			for(index_t i=0; i < (index_t)sp.rows(); ++i) {
				if(sp[i] > 0) {
					alpha[i] = INFINITY;
				}
			}

			Eigen::Matrix<double,Eigen::Dynamic,1>::Index ia;
			double a = alpha.minCoeff(&ia);
			assert(a >= 0 && a <= 1);
			dp = dp + a*(sp-dp);
			for(index_t i=0; i < (index_t)dp.rows(); ++i) {
				if(dp[i] <= 0 || i == ia) {
					P[mapping[i]] = false;
					d[mapping[i]] = 0;
				} else {
					d[mapping[i]] = dp[i];
				}
			}
		} while(true);
	}

	return d;
}
