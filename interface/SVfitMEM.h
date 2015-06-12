#ifndef TauAnalysis_SVfitMEM_SVfitMEM_h
#define TauAnalysis_SVfitMEM_SVfitMEM_h

#include "TauAnalysis/SVfitMEM/interface/SVfitIntegrand.h"
#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <TBenchmark.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TMath.h>

class SVfitMEM
{
 public:
  SVfitMEM(const std::string&, double, const std::string&, int = 0);
  ~SVfitMEM();

  /// take resolution on energy of hadronic tau decays into account
  void shiftVisPt(bool value, TFile* inputFile);

  /// number of function calls for VEGAS integration (default is 10000)
  void setMaxObjFunctionCalls(unsigned value) 
  { 
    numCallsGridOpt_ = TMath::Nint(0.20*value);
    numCallsIntEval_ = TMath::Nint(0.80*value);
  }

  /// run integration 
  void integrate(const std::vector<svFitMEM::MeasuredTauLepton>&, double, double, const TMatrixD&, const std::string& = "");

  /// return mass of the di-tau system 
  double mass() const { return mass_; }
  /// return uncertainty on the mass of the di-tau system
  double massErr() const { return massErr_; }
  /// return maximum of likelihood function
  double Lmax() const { return Lmax_; }  

 protected:

  svFitMEM::SVfitIntegrand* integrand_;
  double sqrtS_;

  std::vector<svFitMEM::MeasuredTauLepton> measuredTauLeptons_;

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
  
  /// mass of di-tau system 
  double mass_;
  /// uncertainty on the mass of di-tau system 
  double massErr_;
  /// maximum of likelihood function
  double Lmax_;

  /// resolution on Pt and mass of hadronic taus
  bool shiftVisPt_;  
  const TH1* lutVisPtResDM0_;
  const TH1* lutVisPtResDM1_;
  const TH1* lutVisPtResDM10_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;

  /// verbosity level
  int verbosity_;
};

#endif
