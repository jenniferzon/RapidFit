#ifndef DP_MASS_SHAPE
#define DP_MASS_SHAPE

#include "TComplex.h"

class DPMassShape
{

  public:

    DPMassShape() {};
    ~DPMassShape() {};

    virtual TComplex massShape(double m) = 0;

    virtual void setResonanceParameters(double mass, double sigma) = 0;

};

#endif