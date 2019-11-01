#ifndef MODELFACTORYDARWIN_H_
#define MODELFACTORYDARWIN_H_

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryDarwin : public ModelFactory<AA> {
public:
	ModelFactoryDarwin();
	virtual double getEpsilon(distance_t distance) const;
	virtual double getDelta(distance_t distance) const;
	virtual ~ModelFactoryDarwin() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MODELFACTORYDARWIN_H_ */
