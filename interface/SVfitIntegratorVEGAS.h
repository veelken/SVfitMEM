#ifndef TauAnalysis_SVfitMEM_SVfitIntegratorVEGAS_h
#define TauAnalysis_SVfitMEM_SVfitIntegratorVEGAS_h

/** \class SVfitIntegratorVEGAS
 *
 * Interface to VEGAS integration algorithm.
 *
 * The VEGAS algorithm is documented in:
 *  [1] "A New Algorithm for Adaptive Multidimensional Integration",
 *      G.P. Lepage, J. Comput. Phys. 27 (1978) 192.
 *
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorBase.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <vector>
#include <string>
#include <iostream>

namespace svFitMEM
{
  class SVfitIntegratorVEGAS : public SVfitIntegratorBase
  {
   public:
    SVfitIntegratorVEGAS(unsigned, unsigned, double, unsigned);
    ~SVfitIntegratorVEGAS();

    void integrate(SVfitIntegratorBase::gPtr, const double*, const double*, unsigned, double&, double&);

    void print(std::ostream&) const {}

   protected:
    void setIntegrand(SVfitIntegratorBase::gPtr, const double*, const double*, unsigned);

    SVfitIntegratorBase::gPtr integrand_;

    gsl_monte_function* vegasIntegrand_;
    gsl_monte_vegas_state* vegasWorkspace_;
    mutable gsl_rng* vegasRnd_;
    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;
    double maxChi2_;
    unsigned maxIntEvalIter_;
    double precision_;
    unsigned numDimensions_;

    /// lower and upper boundary of integration region
    mutable double* xl_;
    mutable double* xu_;
  };
}

#endif

