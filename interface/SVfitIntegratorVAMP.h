#ifndef TauAnalysis_SVfitMEM_SVfitIntegratorVAMP_h
#define TauAnalysis_SVfitMEM_SVfitIntegratorVAMP_h

/** \class SVfitIntegratorVAMP
 *
 * Interface to "Vegas AMPlified: Anisotropy, Multi-channel sampling and Parallelization" (VAMP) integration algorithm.
 *
 * The VAMP algorithm is documented in:
 *  [1] "Vegas revisited: Adaptive Monte Carlo integration beyond factorization",
 *      T. Ohl, J. Comput. Phys. 120 (1999) 13.
 *  [2] https://whizard.hepforge.org/vamp.pdf
 *
 * \author Christian Veelken, NICPB Tallinn
 *
 */

#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorBase.h"

#include <vector>
#include <string>
#include <iostream>

namespace svFitMEM
{
  class SVfitIntegratorVAMP : public SVfitIntegratorBase
  {
   public:
    SVfitIntegratorVAMP(unsigned, unsigned);
    ~SVfitIntegratorVAMP();

    void integrate(SVfitIntegratorBase::gPtr_Fortran, const double*, const double*, unsigned, double&, double&);

    void print(std::ostream&) const {}

   protected:
    void setIntegrand(SVfitIntegratorBase::gPtr_Fortran, const double*, const double*, unsigned);

    SVfitIntegratorBase::gPtr_Fortran integrand_;

    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;
    unsigned numDimensions_;

    /// lower and upper boundary of integration region
    mutable double* xl_;
    mutable double* xu_;
  };
}

#endif

