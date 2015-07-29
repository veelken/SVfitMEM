#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVAMP.h"

#include <assert.h>

using namespace svFitMEM;

extern "C" 
{
  void vamp_integrate_(SVfitIntegratorBase::gPtr_Fortran, const double*, const double*, int*, int*, int*, double*, double*);
}

SVfitIntegratorVAMP::SVfitIntegratorVAMP(unsigned numCallsGridOpt, unsigned numCallsIntEval)
  : integrand_(0),
    numCallsGridOpt_(numCallsGridOpt),
    numCallsIntEval_(numCallsIntEval)
{}

SVfitIntegratorVAMP::~SVfitIntegratorVAMP()
{}

void SVfitIntegratorVAMP::setIntegrand(SVfitIntegratorBase::gPtr_Fortran g, const double* xl, const double* xu, unsigned d)
{
  numDimensions_ = d;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
   for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xl_[iDimension] = xl[iDimension];
    xu_[iDimension] = xu[iDimension];
  }

  integrand_ = g;
}

void SVfitIntegratorVAMP::integrate(SVfitIntegratorBase::gPtr_Fortran g, const double* xl, const double* xu, unsigned d, double& integral, double& integralErr)
{
  setIntegrand(g, xl, xu, d);
  
  if ( !integrand_ ) {
    std::cerr << "<SVfitIntegratorVAMP>:"
	      << "No integrand function has been set yet --> ABORTING !!\n";
    assert(0);
  }

  int numDimensions_int = numDimensions_;
  int numCallsGridOpt_int = numCallsGridOpt_;
  int numCallsIntEval_int = numCallsIntEval_;
  vamp_integrate_(integrand_, xl, xu, &numDimensions_int, &numCallsGridOpt_int, &numCallsIntEval_int, &integral, &integralErr);
  
  delete [] xl_;
  delete [] xu_;
}
