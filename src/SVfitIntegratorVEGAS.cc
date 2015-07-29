#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVEGAS.h"

#include <assert.h>

using namespace svFitMEM;

SVfitIntegratorVEGAS::SVfitIntegratorVEGAS(unsigned numCallsGridOpt, unsigned numCallsIntEval, double maxChi2, unsigned maxIntEvalIter)
  : integrand_(0),
    numCallsGridOpt_(numCallsGridOpt),
    numCallsIntEval_(numCallsIntEval),
    maxChi2_(maxChi2),
    maxIntEvalIter_(maxIntEvalIter),
    precision_(1.e-5)
{}

SVfitIntegratorVEGAS::~SVfitIntegratorVEGAS()
{}

void SVfitIntegratorVEGAS::setIntegrand(SVfitIntegratorBase::gPtr_C g, const double* xl, const double* xu, unsigned d)
{
  numDimensions_ = d;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
   for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xl_[iDimension] = xl[iDimension];
    xu_[iDimension] = xu[iDimension];
  }

  integrand_ = g;

  vegasIntegrand_ = new gsl_monte_function;
  vegasIntegrand_->f = integrand_;
  vegasIntegrand_->dim = numDimensions_;
  vegasIntegrand_->params = new double[1];
  vegasWorkspace_ = gsl_monte_vegas_alloc(numDimensions_);
  gsl_rng_env_setup();
  vegasRnd_ = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(vegasRnd_, 12345); 

  gsl_monte_vegas_init(vegasWorkspace_);
  vegasWorkspace_->stage = 0;
  double integral = 0.;
  double integralErr = 0.;
  gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsGridOpt_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &integral, &integralErr);
  vegasWorkspace_->stage = 1;
}

void SVfitIntegratorVEGAS::integrate(SVfitIntegratorBase::gPtr_C g, const double* xl, const double* xu, unsigned d, double& integral, double& integralErr)
{
  setIntegrand(g, xl, xu, d);
  
  if ( !integrand_ ) {
    std::cerr << "<SVfitIntegratorVEGAS>:"
	      << "No integrand function has been set yet --> ABORTING !!\n";
    assert(0);
  }

  integral = 0.;
  integralErr = 0.;

  // CV: repeat integration in case chi2 of estimated integral/uncertainty values
  //     indicates that result of integration cannot be trusted
  //    (up to maxIntEvalIter times in total)
  unsigned iteration = 0;
  double chi2 = -1.;
  do {
    gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsIntEval_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &integral, &integralErr);
    vegasWorkspace_->stage = 3;
    ++iteration;
    chi2 = vegasWorkspace_->chisq;
  } while ( (chi2 > maxChi2_ || integralErr > (integral*precision_)) && iteration < maxIntEvalIter_ );	
  
  delete [] xl_;
  delete [] xu_;

  delete [] (double*)vegasIntegrand_->params;
  delete vegasIntegrand_;
  gsl_monte_vegas_free(vegasWorkspace_);
  gsl_rng_free(vegasRnd_);
}
