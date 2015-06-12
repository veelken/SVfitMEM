#ifndef TauAnalysis_SVfitMEM_HttXsectionStableTaus_h
#define TauAnalysis_SVfitMEM_HttXsectionStableTaus_h

#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandStableTaus.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <TBenchmark.h>
#include <TMath.h>

class HttXsectionStableTaus
{
 public:
  HttXsectionStableTaus(const std::string&, double, double, const std::string&, int = 0);
  ~HttXsectionStableTaus();

  /// number of function calls for VEGAS integration (default is 10000)
  void setMaxObjFunctionCalls(unsigned value) 
  { 
    numCallsGridOpt_ = TMath::Nint(0.20*value);
    numCallsIntEval_ = TMath::Nint(0.80*value);
  }

  /// run integration 
  void integrate();

  /// return cross-section
  double xSection() const { return xSection_; }
  /// return uncertainty on cross-section
  double xSectionErr() const { return xSectionErr_; }  

  double GammaH() const { return integrand_->GammaH(); }

 protected:

  svFitMEM::HttXsectionIntegrandStableTaus* integrand_;
  double sqrtS_;
  double mH_;

  /// auxiliary variables for VEGAS integration
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
  double* xl_;
  double* xu_;
  
  /// cross-section
  double xSection_;
  /// uncertainty on cross-section
  double xSectionErr_;
  
  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;

  /// verbosity level
  int verbosity_;
};

#endif
