#ifndef TauAnalysis_SVfitMEM_SVfitMEM_lo_h
#define TauAnalysis_SVfitMEM_SVfitMEM_lo_h

#include "TauAnalysis/SVfitMEM/interface/SVfitIntegrand_lo.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorBase.h"
#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"

#include <TBenchmark.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include <TMath.h>

class SVfitMEM_lo
{
 public:
  SVfitMEM_lo(double, const std::string&, int = svFitMEM::SVfitIntegrand_lo::kLiterature, const std::string& = "", int = 0); 
  ~SVfitMEM_lo();

  /// take cross-section*signal acceptance/efficiency into account
  void setCrossSection(const TGraphErrors*);  
  void setCrossSection_and_Acc(const TGraphErrors*, const TGraphErrors*, double);

  /// add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
  void addLogM(bool value, double power = 1.) 
  { 
    addLogM_ = value; 
    addLogM_power_ = power; 
  }

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

  /// set correlation between hadronic tau pT and MET
  void setRhoHadTau(double rhoHadTau) 
  { 
    integrand_->setRhoHadTau(rhoHadTau);
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
  
  /// run integration 
  void integrate(const std::vector<svFitMEM::MeasuredTauLepton>&, double, double, const TMatrixD&, const std::string& = "");

  /// return mass of the di-tau system 
  double mass() const { return mass_; }
  /// return uncertainty on the mass of the di-tau system
  double massErr() const { return massErr_; }
  /// return maximum of likelihood function
  double Lmax() const { return Lmax_; }  
  /// return flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution() { return isValidSolution_; }

  enum { kMarkovChain, kVEGAS, kVAMP }; 

 protected:

  svFitMEM::SVfitIntegrand_lo* integrand_;
  double sqrtS_;

  std::vector<svFitMEM::MeasuredTauLepton> measuredTauLeptons_;

  /// interface to integration algorithm (either Markov Chain integration, VEGAS or VAMP)
  int intMode_;
  svFitMEM::SVfitIntegratorBase* intAlgo_;
  unsigned maxObjFunctionCalls_;
  double precision_;

  /// dimension of integration region
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
  /// flag indicating if algorithm succeeded to find valid solution
  bool isValidSolution_;

  /// cross-section*signal acceptance/efficiency as function of mass
  const TGraphErrors* graph_xSection_;
  const TGraphErrors* graph_Acc_;
  double minAcc_;

  /// flag to enable/disable addition of log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution 
  bool addLogM_; 
  double addLogM_power_; 

  /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
  bool useHadTauTF_;

  /// clock for measuring run-time of algorithm
  TBenchmark* clock_;

  /// verbosity level
  int verbosity_;
};

#endif
