#ifndef TauAnalysis_SVfitMEM_me_ggH_lit_h
#define TauAnalysis_SVfitMEM_me_ggH_lit_h

/**
   \class me_ggH_lit me_ggH_lit.h "TauAnalysis/SVfitStandalone/interface/me_ggH_lit.h"
   \brief Compute leading order gg -> Higgs matrix element. 
      The formulas are taken from hep-ph/9504378 and from S. Bentvelsen, E. Laenen and P. Motylinski: "Higgs production through gluon fusion at leading order".
*/

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"

class me_ggH_lit
{
 public:
  // Constructor and destructor
  me_ggH_lit(bool = false, bool = true);
  ~me_ggH_lit();

  void setS(double s) 
  { 
    s_ = s;
    s_isInitialized_ = true;
  }
  void setHiggsMass(double mH) 
  { 
    mH_ = mH;
    mH2_ = svFitMEM::square(mH_);
    mH_isInitialized_ = true;
  }
  void setHiggsWidth(double width) 
  {
    GammaH_ = width;
    GammaH2_ = svFitMEM::square(GammaH_);
    GammaH_isInitialized_ = true;
  }
  void setBR(double br) 
  {
    br_ = br;
    br_isInitialized_ = true;
  }
  
  // Get and set momenta for matrix element evaluation
  void setMomenta(std::vector<double*>& momenta){ momenta_ = momenta; }

  // Get matrix element.
  double getMatrixElement() const;
  double getMatrixElement_woHtoTauTauDecay() const;

 private:
  // flag to enable/disable narrow-width approximation
  bool applyNWA_;

  // flag to include/not include Higgs -> tautau decay in matrix element
  bool includeHtoTauTauDecay_;

  // center-of-mass energy
  double s_;
  bool s_isInitialized_;

  // Higgs mass, width and decay braching fraction
  double mH_;
  double mH2_;
  bool mH_isInitialized_;
  double GammaH_;
  double GammaH2_;
  bool GammaH_isInitialized_;
  double br_;
  bool br_isInitialized_;

  // Mass of tau lepton pair
  mutable double q_;
  mutable double q2_;

  // vector with momenta (to be changed each event)
  std::vector<double*> momenta_;
}; 

#endif 
