#ifndef TauAnalysis_SVfitMEM_HttXsectionWithTauDecays_h
#define TauAnalysis_SVfitMEM_HttXsectionWithTauDecays_h

#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandWithTauDecays.h"
#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <TBenchmark.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TMath.h>

class HttXsectionWithTauDecays
{
 public:
  HttXsectionWithTauDecays(const std::string&, double, double, const std::string&, int = 0);
  ~HttXsectionWithTauDecays();

  /// set Higgs -> tautau decay branching fraction
  void setBR(double br) { integrand_->setBR(br); }

  /// take resolution on energy of hadronic tau decays into account
  void shiftVisPt(bool value, TFile* inputFile);

  /// number of function calls for VEGAS integration (default is 10000)
  void setMaxObjFunctionCalls(unsigned value) 
  { 
    numCallsGridOpt_ = TMath::Nint(0.20*value);
    numCallsIntEval_ = TMath::Nint(0.80*value);
  }

  /// enable/disable acceptance cuts
  void enableAcceptanceCuts(double (*acceptance)(const svFitMEM::LorentzVector&, const svFitMEM::LorentzVector&, double, double))
  {
    integrand_->enableAcceptanceCuts(acceptance); 
  }
  void disableAcceptanceCuts()
  {
    integrand_->disableAcceptanceCuts(); 
  }
    
  /// run integration 
  void integrate(int, int, double, int, int, double, const TMatrixD&);

  /// return cross-section
  double xSection() const { return xSection_; }
  /// return uncertainty on cross-section
  double xSectionErr() const { return xSectionErr_; }  

 protected:

  svFitMEM::HttXsectionIntegrandWithTauDecays* integrand_;
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
