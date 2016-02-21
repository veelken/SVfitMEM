#ifndef TauAnalysis_SVfitMEM_me_ggH_nlo_lit_h
#define TauAnalysis_SVfitMEM_me_ggH_nlo_lit_h

/**
   \class me_ggH_nlo_lit me_ggH_nlo_lit.h "TauAnalysis/SVfitStandalone/interface/me_ggH_nlo_lit.h"
   \brief Compute leading order gg -> Higgs+jet matrix element. 
      The formulas are taken from 
        arXiv:hep-ph/0201114
        "Next-to-leading order QCD corrections to differential distributions of Higgs boson production in hadron-hadron collisions", V. Ravindran, J. Smith and W.L. Van Neerven
*/

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"

#include "LHAPDF/LHAPDF.h"

class me_ggH_nlo_lit
{
 public:
  /// Constructor and destructor
  me_ggH_nlo_lit(const std::string&, bool = false, bool = true);
  ~me_ggH_nlo_lit();

  void setS(double);
  void setHiggsMass(double);
  void setHiggsWidth(double);
  void setBR(double);
  
  /// Get and set momenta for matrix element evaluation
  void setMomenta(std::vector<const svFitMEM::LorentzVector*>&);

  /// Get matrix element.
  double getMatrixElement_gg() const;
  double getMatrixElement_gg_woHtoTauTauDecay() const;
  double getMatrixElement_qq() const;
  double getMatrixElement_qq_woHtoTauTauDecay() const;
  double getMatrixElement_qg() const;
  double getMatrixElement_qg_woHtoTauTauDecay() const;

 private:
  /// Get matrix element (auxiliary function for Higgs -> tautau decay)
  double getMatrixElement_HtoTauTauDecay() const;

  /// flag to enable/disable narrow-width approximation
  bool applyNWA_;

  /// flag to include/not include Higgs -> tautau decay in matrix element
  bool includeHtoTauTauDecay_;

  /// center-of-mass energy
  double s_;
  bool s_isInitialized_;

  /// Higgs mass, width and decay braching fraction
  double mH_;
  double mH2_;
  double mH8_;
  bool mH_isInitialized_;
  double GammaH_;
  double GammaH2_;
  bool GammaH_isInitialized_;
  double br_;
  bool br_isInitialized_;

  /// Mass of tau lepton pair
  mutable double q_;
  mutable double q2_;

  /// vector with momenta (to be changed each event)
  std::vector<const svFitMEM::LorentzVector*> momenta_;

  /// PDF (needed to access alphaS)
  LHAPDF::PDF* pdf_;
  bool pdfIsInitialized_;

  double mandelstam_s_;
  double mandelstam_abs_sMin_;
  double mandelstam_s2_;
  double mandelstam_s4_;
  double mandelstam_s8_;
  double mandelstam_t_;
  double mandelstam_abs_tMin_;
  double mandelstam_t2_;
  double mandelstam_t4_;
  double mandelstam_t8_;
  double mandelstam_u_;
  double mandelstam_abs_uMin_;
  double mandelstam_u2_;
  double mandelstam_u4_;
  double mandelstam_u8_;

  double K_gg_;
  double K_qq_;
  double K_qg_;

  double G2_;
  double g2_;
}; 

#endif 
