#ifndef TauAnalysis_SVfitMEM_HttXsectionWithTauDecays2_h
#define TauAnalysis_SVfitMEM_HttXsectionWithTauDecays2_h

#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandWithTauDecays2.h"
#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorBase.h"

#include <TBenchmark.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TMath.h>

class HttXsectionWithTauDecays2
{
 public:
  HttXsectionWithTauDecays2(double, double, const std::string&, int = HttXsectionIntegrandWithTauDecays2::kLiterature, const std::string& = "", int = 0); 
  ~HttXsectionWithTauDecays2();

  /// set Higgs -> tautau decay branching fraction
  void setBR(double br) { integrand_->setBR(br); }

  /// take resolution on energy of hadronic tau decays into account
  void shiftVisPt(bool value, TFile* inputFile);

  /// number of function calls for Markov Chain and VEGAS integration (default is 100000)
  void setMaxObjFunctionCalls(unsigned maxObjFunctionCalls) 
  { 
    maxObjFunctionCalls_ = maxObjFunctionCalls;
  }

  /// set integration algorithm (either Markov Chain integration, VEGAS or VAMP algorithm)
  void setIntMode(int intMode)
  {
    intMode_ = intMode;
  }

  /// enable/disable acceptance cuts
  void enableAcceptanceCuts(const acceptanceBaseType& acceptance)
  {
    integrand_->enableAcceptanceCuts(acceptance);
    applyAcceptanceCuts_ = true;
  }
  void disableAcceptanceCuts()
  {
    integrand_->disableAcceptanceCuts(); 
    applyAcceptanceCuts_ = false;
  }

  void enableZeroPtConstraint() { applyZeroPtConstraint_ = true; }
  void disableZeroPtConstraint() { applyZeroPtConstraint_ = false; }
  void setMax_v_r_(double max_v_r) { max_v_r_ = max_v_r; }
    
  /// run integration 
  void integrate(int, int, double, int, int, double, const TMatrixD&);

  /// return cross-section
  double xSection() const { return xSection_; }
  /// return uncertainty on cross-section
  double xSectionErr() const { return xSectionErr_; }  

  enum { kMarkovChain, kVEGAS, kVAMP }; 

 protected:

  bool applyMEtTF_;
  bool applyAcceptanceCuts_;
  bool applyZeroPtConstraint_;
  double max_v_r_;

  HttXsectionIntegrandWithTauDecays2* integrand_;
  double sqrtS_;
  double mH_;
  double mH2_;

  /// interface to integration algorithm (either Markov Chain integration, VEGAS or VAMP algorithm)
  int intMode_;
  svFitMEM::SVfitIntegratorBase* intAlgo_;
  unsigned maxObjFunctionCalls_;

  /// dimension of integration region
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
