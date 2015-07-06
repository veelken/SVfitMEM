#include "TauAnalysis/SVfitMEM/interface/HttXsectionStableHiggs.h"

#include <TGraphErrors.h>
#include <TH1.h>

#include <algorithm>

using namespace svFitMEM;

namespace 
{
  double g(double* x, size_t dim, void* param)
  {    
    return HttXsectionIntegrandStableHiggs::gHttXsectionIntegrandStableHiggs->Eval(x);
  }
}

HttXsectionStableHiggs::HttXsectionStableHiggs(const std::string& madgraphFileName, double sqrtS, double mH, const std::string& pdfFileName, bool applyNWA, int verbosity) 
  : integrand_(0),
    sqrtS_(sqrtS),
    mH_(mH),
    applyNWA_(applyNWA),
    vegasIntegrand_(0),
    vegasWorkspace_(0),
    vegasRnd_(0),
    numCallsGridOpt_(10000),
    numCallsIntEval_(40000),
    //numCallsGridOpt_(2000),
    //numCallsIntEval_(8000),
    maxChi2_(2.),
    //maxIntEvalIter_(20),
    maxIntEvalIter_(5),
    precision_(1.e-2),
    xl_(0),
    xu_(0),
    verbosity_(verbosity)
{ 
  integrand_ = new HttXsectionIntegrandStableHiggs(madgraphFileName, sqrtS_, mH_, pdfFileName, applyNWA_, verbosity_);

  numDimensions_ = ( applyNWA_ ) ? 1 : 2;

  clock_ = new TBenchmark();
}

HttXsectionStableHiggs::~HttXsectionStableHiggs() 
{
  delete integrand_;

  delete clock_;
}

void
HttXsectionStableHiggs::integrate()
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<HttXsectionStableHiggs::integrate>:" << std::endl;
  }
  if ( verbosity_ >= 0 ) {
    clock_->Start("<HttXsectionStableHiggs::integrate>");
  }
  
  HttXsectionIntegrandStableHiggs::gHttXsectionIntegrandStableHiggs = integrand_;

  vegasIntegrand_ = new gsl_monte_function;
  vegasIntegrand_->f = &g;
  vegasIntegrand_->dim = numDimensions_;
  vegasIntegrand_->params = new double[1];
  vegasWorkspace_ = gsl_monte_vegas_alloc(numDimensions_);
  gsl_rng_env_setup();
  vegasRnd_ = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(vegasRnd_, 12345); 
  
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
  xl_[0] = 0.;
  xu_[0] = 1.;
  if ( !applyNWA_ ) {
    xl_[1] = -0.5*TMath::Pi(); // tk, as defined by Eq. (8) in arXiv:1010.2263v3
    xu_[1] = +0.5*TMath::Pi();
  }

  gsl_monte_vegas_init(vegasWorkspace_);
  xSection_ = 0.;
  xSectionErr_ = 0.;
  vegasWorkspace_->stage = 0;
  gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsGridOpt_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &xSection_, &xSectionErr_);
  xSection_ = 0.;
  xSectionErr_ = 0.;
  vegasWorkspace_->stage = 1;

  // CV: repeat integration in case chi2 of estimated integral/uncertainty values
  //     indicates that result of integration cannot be trusted
  //    (up to maxIntEvalIter times in total)
  unsigned iteration = 0;
  double chi2 = -1.;
  do {
    gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsIntEval_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &xSection_, &xSectionErr_);
    vegasWorkspace_->stage = 3;
    ++iteration;
    chi2 = vegasWorkspace_->chisq;
  } while ( (chi2 > maxChi2_ || xSectionErr_ > (xSection_*precision_)) && iteration < maxIntEvalIter_ );	
  
  if ( verbosity_ >= 1 ) {
    std::cout << "--> cross-section = " << xSection_ << " +/- " << xSectionErr_ << std::endl;
  }

  delete [] xl_;
  delete [] xu_;

  delete [] (double*)vegasIntegrand_->params;
  delete vegasIntegrand_;
  gsl_monte_vegas_free(vegasWorkspace_);
  gsl_rng_free(vegasRnd_);

  if ( verbosity_ >= 0 ) {
    clock_->Show("<HttXsectionStableHiggs::integrate>");
  }
}
