#ifndef TauAnalysis_SVfitMEM_HttXsectionIntegrandStableTaus_h
#define TauAnalysis_SVfitMEM_HttXsectionIntegrandStableTaus_h

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_mg5.h"

#include <TMatrixD.h>

namespace svFitMEM
{
  class HttXsectionIntegrandStableTaus 
  {
   public:
    HttXsectionIntegrandStableTaus(const std::string&, double, double, const std::string&, int);
    ~HttXsectionIntegrandStableTaus();
  
    /// evaluate integrand for given value of integration variables x
    double Eval(const double* x) const;

    // get Higgs width from Madgraph
    double GammaH() const { return madgraph_.getHiggsWidth(); }	

    /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS integration)
    static const HttXsectionIntegrandStableTaus* gHttXsectionIntegrandStableTaus;

   protected:
    double compProb(double, double, double, double, double, double, double) const;

    double mH_;
    double mH2_;

    double s_;
    double invSqrtS_;
    Vector beamAxis_;

    static bool pdfIsInitialized_;

    mutable me_ggH_mg5 madgraph_;
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
