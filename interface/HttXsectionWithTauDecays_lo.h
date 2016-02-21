#ifndef TauAnalysis_SVfitMEM_HttXsectionWithTauDecays_lo_h
#define TauAnalysis_SVfitMEM_HttXsectionWithTauDecays_lo_h

#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandWithTauDecays_lo.h"
#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorBase.h"

#include <TBenchmark.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TMath.h>

class HttXsectionWithTauDecays_lo
{
 public:
  HttXsectionWithTauDecays_lo(double, double, const std::string&, int = HttXsectionIntegrandWithTauDecays_lo::kLiterature, const std::string& = "", int = 0); 
  ~HttXsectionWithTauDecays_lo();

  /// set Higgs -> tautau decay branching fraction
  void setBR(double br) { integrand_->setBR(br); }

  /// set transfer functions for pT of hadronic tau decays
  void setHadTauTF(const HadTauTFBase* hadTauTF) 
  { 
    integrand_->setHadTauTF(hadTauTF);
  }
  /// enable/disable use of transfer functions for hadronic tau decays
  void enableHadTauTF() 
  { 
    integrand_->enableHadTauTF();
    useHadTauTF_ = true; 
  }
  void disableHadTauTF() 
  { 
    integrand_->disableHadTauTF();
    useHadTauTF_ = false; 
  }

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

  HttXsectionIntegrandWithTauDecays_lo* integrand_;
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
  
  /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
  bool useHadTauTF_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;

  /// verbosity level
  int verbosity_;
};

#endif
