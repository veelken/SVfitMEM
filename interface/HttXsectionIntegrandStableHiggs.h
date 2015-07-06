#ifndef TauAnalysis_SVfitMEM_HttXsectionIntegrandStableHiggs_h
#define TauAnalysis_SVfitMEM_HttXsectionIntegrandStableHiggs_h

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_mg5.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lit.h"

#include <TMatrixD.h>

namespace svFitMEM
{
  class HttXsectionIntegrandStableHiggs
  {
   public:
    HttXsectionIntegrandStableHiggs(const std::string&, double, double, const std::string&, bool, int);
    ~HttXsectionIntegrandStableHiggs();
  
    /// evaluate integrand for given value of integration variables x
    double Eval(const double* x) const;

    // get Higgs width from Madgraph
    double GammaH() const { return me_madgraph_.getHiggsWidth(); }	

    /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS integration)
    static const HttXsectionIntegrandStableHiggs* gHttXsectionIntegrandStableHiggs;

   protected:
    double mH_;
    double mH2_;

    double sqrtS_;
    double s_;
    double invSqrtS_;

    bool applyNWA_;

    Vector beamAxis_;

    static bool pdfIsInitialized_;

    mutable me_ggH_mg5 me_madgraph_;
    mutable me_ggH_lit me_lit_;
    double* madgraphGluon1P4_;
    double* madgraphGluon2P4_;
    double* madgraphTau1P4_;
    double* madgraphTau2P4_;
    mutable vector<double*> madgraphMomenta_;

    /// verbosity level
    int verbosity_;
  };
}

#endif
