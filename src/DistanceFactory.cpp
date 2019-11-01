#include "main.h"
#include "DistanceFactory.h"
#include "ModelFactory.h"
#include "DistanceFactoryAlign.h"
#include "DistanceFactoryAngle.h"
#include "DistanceFactoryPrealigned.h"

template <>
DistanceFactory<AA> *DistanceFactory<AA>::getDefault(const ModelFactory<AA> *model_factory, bool prealigned) {
	DistanceFactory<AA> *dist_factory = NULL;
	if(!prealigned) {
		if(cmdlineopts.nwdist_flag) {
			dist_factory = new DistanceFactoryAlign<AA>(model_factory);
		} else {
			dist_factory = new DistanceFactoryAngle<AA,2>();
		}
	} else {
		dist_factory = new DistanceFactoryPrealigned<AA>(model_factory);
	}
	return dist_factory;
}

#ifdef WITH_DNA
template <>
DistanceFactory<DNA> *DistanceFactory<DNA>::getDefault(const ModelFactory<DNA> *model_factory, bool prealigned) {
	DistanceFactory<DNA> *dist_factory = NULL;
	if(!prealigned) {
		if(cmdlineopts.nwdist_flag) {
			dist_factory = new DistanceFactoryAlign<DNA>(model_factory);
		} else {
			dist_factory = new DistanceFactoryAngle<DNA,6>();
		}
	} else {
		dist_factory = new DistanceFactoryPrealigned<DNA>(model_factory);
	}
	return dist_factory;
}
#endif

#ifdef WITH_CODON
template <>
DistanceFactory<Codon> *DistanceFactory<Codon>::getDefault(const ModelFactory<Codon> *model_factory, bool prealigned) {
	DistanceFactory<Codon> *dist_factory = NULL;
	if(!prealigned) {
		if(cmdlineopts.nwdist_flag) {
			dist_factory = new DistanceFactoryAlign<Codon>(model_factory);
		} else {
			dist_factory = new DistanceFactoryAngle<Codon,2>();
		}
	} else {
		dist_factory = new DistanceFactoryPrealigned<Codon>(model_factory);
	}
	return dist_factory;
}
#endif
