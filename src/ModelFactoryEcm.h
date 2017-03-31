#ifndef MODELFACTORYECM_H_
#define MODELFACTORYECM_H_

#ifdef WITH_CODON

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryEcm : public ModelFactory<Codon> {
public:
	ModelFactoryEcm();
	virtual ~ModelFactoryEcm() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#endif

#endif /* MODELFACTORYECM_H_ */
