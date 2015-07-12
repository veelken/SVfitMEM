#ifndef TauAnalysis_SVfitMEM_SVfitIntegratorBase_h
#define TauAnalysis_SVfitMEM_SVfitIntegratorBase_h

/** \class SVfitIntegratorBase
 *
 * Base class for Markov Chain and VEGAS integration classes.
 *
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include <iostream>

namespace svFitMEM
{
  class SVfitIntegratorBase
  {
   public:
    SVfitIntegratorBase() {}
    ~SVfitIntegratorBase() {}

    /// compute integral of function g 
    /// the points xl and xh represent the lower left and upper right corner of a Hypercube in d-dimensional integration space 
    typedef double (*gPtr)(double*, size_t, void*);
    virtual void integrate(gPtr, const double*, const double*, unsigned, double&, double&) = 0;

    virtual void print(std::ostream&) const {}
  };
}

#endif

