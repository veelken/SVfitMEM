#ifndef TauAnalysis_SVfitMEM_HttXsectionIntegrandWithTauDecays_h
#define TauAnalysis_SVfitMEM_HttXsectionIntegrandWithTauDecays_h

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_mg5.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lit.h"

#include <TMatrixD.h>

class HttXsectionIntegrandWithTauDecays 
{
 public:
  HttXsectionIntegrandWithTauDecays(double, double, const std::string&, int, const std::string&, int);
  ~HttXsectionIntegrandWithTauDecays();
  
  /// set Higgs -> tautau decay branching fraction
  void setBR(double br) { me_lit_.setBR(br); }

  /// enable/disable acceptance cuts
  void enableAcceptanceCuts(double (*acceptance)(const svFitMEM::LorentzVector&, const svFitMEM::LorentzVector&, double, double))
  {
    //std::cout << "<HttXsectionIntegrandWithTauDecays::enableAcceptanceCuts>: acceptance = " << &acceptance << std::endl;
    acceptance_ = acceptance;      
  }
  void disableAcceptanceCuts()
  {
    //std::cout << "<HttXsectionIntegrandWithTauDecays::disableAcceptanceCuts>:" << std::endl;
    acceptance_ = 0;
  }
  
  void setIdxLeg1_X(int idx) { idxLeg1_X_ = idx; }
  void setIdxLeg1_phi(int idx) { idxLeg1_phi_ = idx; }
  void setIdxLeg1VisPtShift(int idx) { idxLeg1VisPtShift_ = idx; }
  void setIdxLeg1_mNuNu(int idx) { idxLeg1_mNuNu_ = idx; }
  void setIdxLeg2_t(int idx) { idxLeg2_t_ = idx; }
  void setIdxLeg2_phi(int idx) { idxLeg2_phi_ = idx; }
  void setIdxLeg2VisPtShift(int idx) { idxLeg2VisPtShift_ = idx; }
  void setIdxLeg2_mNuNu(int idx) { idxLeg2_mNuNu_ = idx; }
  
  /// apply/do not apply transfer function for missing transverse energy
  void setApplyMEtTF(bool applyMEtTF) { applyMEtTF_ = applyMEtTF; }  

  /// take resolution on energy and mass of hadronic tau decays into account
  void shiftVisPt(bool, const TH1*, const TH1*);
    
  /// set momenta of visible tau decay products and of reconstructed missing transverse energy
  void setInputs(int, double, int, double, const TMatrixD&);

  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;

  /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS integration)
  static const HttXsectionIntegrandWithTauDecays* gHttXsectionIntegrandWithTauDecays;

  enum { kMadgraph, kLiterature };

 protected:
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
  
  bool shiftVisPt_;
  const TH1* leg1lutVisPtRes_;
  const TH1* leg2lutVisPtRes_;
  
  int idxLeg1_X_;
  int idxLeg1_phi_;
  int idxLeg1VisPtShift_;
  int idxLeg1_mNuNu_;
  int idxLeg2_t_;
  int idxLeg2_phi_;
  int idxLeg2VisPtShift_;
  int idxLeg2_mNuNu_;
  
  static bool pdfIsInitialized_;

  mutable me_ggH_mg5 me_madgraph_;
  bool me_madgraph_isInitialized_;
  mutable me_ggH_lit me_lit_;
  double* madgraphGluon1P4_;
  double* madgraphGluon2P4_;
  double* madgraphTau1P4_;
  double* madgraphTau2P4_;
  mutable vector<double*> madgraphMomenta_;
  
  /// acceptance function
  double (*acceptance_)(const svFitMEM::LorentzVector&, const svFitMEM::LorentzVector&, double, double);

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
