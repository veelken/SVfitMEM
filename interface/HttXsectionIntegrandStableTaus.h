#ifndef TauAnalysis_SVfitMEM_HttXsectionIntegrandStableTaus_h
#define TauAnalysis_SVfitMEM_HttXsectionIntegrandStableTaus_h

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_mg5.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lit.h"

#include "LHAPDF/LHAPDF.h"

#include <TMatrixD.h>

class HttXsectionIntegrandStableTaus 
{
 public:
  HttXsectionIntegrandStableTaus(double, double, const std::string&, int, const std::string&, int);
  ~HttXsectionIntegrandStableTaus();

  /// set Higgs -> tautau decay branching fraction
  void setBR(double br) { me_lit_.setBR(br); }
  
  /// evaluate integrand for given value of integration variables x
  double Eval(const double* x) const;
    
  // get Higgs width from Madgraph
  double GammaH() const { return GammaH_; }
  
  /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS integration)
  static const HttXsectionIntegrandStableTaus* gHttXsectionIntegrandStableTaus;
    
  enum { kMadgraph, kLiterature };

 protected:
  double compProb(double, double, double, double, double, double, double) const;

  int mode_;

  double mH_;
  double mH2_;
  double GammaH_;
  
  double sqrtS_;
  double s_;
  double invSqrtS_;
  svFitMEM::Vector beamAxis_;
  
  LHAPDF::PDF* pdf_;
  bool pdfIsInitialized_;
  
  mutable me_ggH_mg5 me_madgraph_;
  bool me_madgraph_isInitialized_;
  mutable me_ggH_lit me_lit_;
  double* madgraphGluon1P4_;
  double* madgraphGluon2P4_;
  double* madgraphTau1P4_;
  double* madgraphTau2P4_;
  mutable vector<double*> madgraphMomenta_;
  
  /// verbosity level
  int verbosity_;
};

#endif
