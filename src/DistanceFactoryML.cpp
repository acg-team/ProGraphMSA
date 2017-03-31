#include "DistanceFactoryML.h"

template <>
const double DistanceFactoryML<AA>::DIST_MAX = 2.2;

template <>
const double DistanceFactoryML<AA>::VAR_MAX = 1e3;

template <>
const double DistanceFactoryML<AA>::VAR_MIN = 1e-5;

#ifdef WITH_CODON
template <>
const double DistanceFactoryML<Codon>::DIST_MAX = 5.2;

template <>
const double DistanceFactoryML<Codon>::VAR_MAX = 5e3;

template <>
const double DistanceFactoryML<Codon>::VAR_MIN = 1e-5;
#endif

#ifdef WITH_DNA
template <>
const double DistanceFactoryML<DNA>::DIST_MAX = 2.2;

template <>
const double DistanceFactoryML<DNA>::VAR_MAX = 1e3;

template <>
const double DistanceFactoryML<DNA>::VAR_MIN = 1e-5;
#endif
