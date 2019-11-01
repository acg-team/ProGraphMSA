#ifndef MODELFACTORYWAG_H_
#define MODELFACTORYWAG_H_

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryWag : public ModelFactory<AA> {
public:
	ModelFactoryWag();
	virtual ~ModelFactoryWag() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MODELFACTORYWAG_H_ */
