#ifndef TauAnalysis_SVfitMEM_HttXsectionIntegrandWithTauDecays_lo_h
#define TauAnalysis_SVfitMEM_HttXsectionIntegrandWithTauDecays_lo_h

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lo_mg5.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lo_lit.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"

#include "LHAPDF/LHAPDF.h"

#include <TMatrixD.h>

class acceptanceBaseType
{
 public:
  virtual double operator()(const svFitMEM::LorentzVector&, const svFitMEM::LorentzVector&, double, double) const = 0;
};

class HttXsectionIntegrandWithTauDecays_lo
{
 public:
  HttXsectionIntegrandWithTauDecays_lo(double, double, const std::string&, int, const std::string&, int);
  ~HttXsectionIntegrandWithTauDecays_lo();
  
  /// set Higgs -> tautau decay branching fraction
  void setBR(double br) { me_lit_.setBR(br); }

  /// enable/disable acceptance cuts
  void enableAcceptanceCuts(const acceptanceBaseType& acceptance)
  {
    //std::cout << "<HttXsectionIntegrandWithTauDecays::enableAcceptanceCuts>: acceptance = " << &acceptance << std::endl;
    acceptance_ = &acceptance;      
  }
  void disableAcceptanceCuts()
  {
    //std::cout << "<HttXsectionIntegrandWithTauDecays::disableAcceptanceCuts>:" << std::endl;
    acceptance_ = 0;
  }
  
  void setIdxLeg1_X(int idx) { idxLeg1_X_ = idx; updateNumDimensions(); }
  void setIdxLeg1_phi(int idx) { idxLeg1_phi_ = idx; updateNumDimensions(); }
  void setIdxLeg1VisPtShift(int idx) { idxLeg1VisPtShift_ = idx; updateNumDimensions(); }
  void setIdxLeg1_mNuNu(int idx) { idxLeg1_mNuNu_ = idx; updateNumDimensions(); }
  void setIdxLeg2_t(int idx) { idxLeg2_t_ = idx; updateNumDimensions(); }
  void setIdxLeg2_phi(int idx) { idxLeg2_phi_ = idx; updateNumDimensions(); }
  void setIdxLeg2VisPtShift(int idx) { idxLeg2VisPtShift_ = idx; updateNumDimensions(); }
  void setIdxLeg2_mNuNu(int idx) { idxLeg2_mNuNu_ = idx; updateNumDimensions(); }
  
  /// apply/do not apply transfer function for missing transverse energy
  void setApplyMEtTF(bool applyMEtTF) { applyMEtTF_ = applyMEtTF; }  

  /// set transfer functions for pT of hadronic tau decays
  void setHadTauTF(const HadTauTFBase* hadTauTF) 
  { 
    hadTauTF1_ = hadTauTF->Clone("leg1"); 
    hadTauTF2_ = hadTauTF->Clone("leg2");
  }
  /// enable/disable use of transfer functions for hadronic tau decays
  void enableHadTauTF() 
  { 
    if ( !(hadTauTF1_ && hadTauTF2_) ) {
      std::cerr << "No tau pT transfer functions defined, call 'setHadTauTF' function first !!" << std::endl;
      assert(0);
    }      
    useHadTauTF_ = true; 
  }
  void disableHadTauTF() 
  { 
    useHadTauTF_ = false; 
  }
    
  /// set momenta of visible tau decay products and of reconstructed missing transverse energy
  void setInputs(int, double, int, double, const TMatrixD&);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

  /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS integration)
  static const HttXsectionIntegrandWithTauDecays_lo* gHttXsectionIntegrandWithTauDecays;
  static int gNumInstances;

  enum { kMadgraph, kLiterature };

 protected:
  void updateNumDimensions();
  double compProb(const double*, const svFitMEM::LorentzVector&, double, const svFitMEM::LorentzVector&, double, double) const;

  /// measured tau leptons
  bool leg1isLep_;
  double leg1Mass_;
  double leg1Mass2_;
  mutable svFitMEM::Vector eZ1_;
  mutable svFitMEM::Vector eY1_;
  mutable svFitMEM::Vector eX1_; 
  mutable double leg1eX_x_;
  mutable double leg1eX_y_;
  mutable double leg1eX_z_;
  mutable double leg1eY_x_;
  mutable double leg1eY_y_;
  mutable double leg1eY_z_;
  mutable double leg1eZ_x_;
  mutable double leg1eZ_y_;
  mutable double leg1eZ_z_;
  bool leg2isLep_;
  double leg2Mass_;
  double leg2Mass2_;
  mutable svFitMEM::Vector eZ2_;
  mutable svFitMEM::Vector eY2_;
  mutable svFitMEM::Vector eX2_; 
  mutable double leg2eX_x_;
  mutable double leg2eX_y_;
  mutable double leg2eX_z_;
  mutable double leg2eY_x_;
  mutable double leg2eY_y_;
  mutable double leg2eY_z_;
  mutable double leg2eZ_x_;
  mutable double leg2eZ_y_;
  mutable double leg2eZ_z_;
  
  int mode_;

  bool applyMEtTF_;

  double mH_;
  double mH2_;
  mutable double GammaH_;
  mutable double q2_;
  mutable double GammaH_times_mH_;
  mutable double GammaH2_times_mH2_;
  
  double sqrtS_;
  double s_;
  double invSqrtS_;
  svFitMEM::Vector beamAxis_;
  
  /// inverse of MET covariance matrix
  TMatrixD invCovMET_;
  double invCovMETxx_;
  double invCovMETxy_;
  double invCovMETyx_;
  double invCovMETyy_;
  double const_MET_;
  
  /// account for resolution on pT of hadronic tau decays via appropriate transfer functions
  const HadTauTFBase* hadTauTF1_;
  const HadTauTFBase* hadTauTF2_;
  bool useHadTauTF_;
  
  int idxLeg1_X_;
  int idxLeg1_phi_;
  int idxLeg1VisPtShift_;
  int idxLeg1_mNuNu_;
  int idxLeg2_t_;
  int idxLeg2_phi_;
  int idxLeg2VisPtShift_;
  int idxLeg2_mNuNu_;
  int numDimensions_;

  LHAPDF::PDF* pdf_;
  bool pdfIsInitialized_;

  mutable me_ggH_lo_mg5 me_madgraph_;
  bool me_madgraph_isInitialized_;
  mutable me_ggH_lo_lit me_lit_;
  double* madgraphGluon1P4_;
  double* madgraphGluon2P4_;
  double* madgraphTau1P4_;
  double* madgraphTau2P4_;
  mutable vector<double*> madgraphMomenta_;
  
  /// acceptance function
  const acceptanceBaseType* acceptance_;

  /// data-members to record minimum and maximum of integrand
  mutable double minIntegrand_;
  mutable double maxIntegrand_;
  mutable bool isFirstCall_;
  mutable TH1* histogramLogIntegrand_;
  std::string outputFileName_;

  /// verbosity level
  int verbosity_;
};

#endif
