#ifndef TauAnalysis_SVfitMEM_SVfitIntegrand_lo_h
#define TauAnalysis_SVfitMEM_SVfitIntegrand_lo_h

#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lo_mg5.h"
#include "TauAnalysis/SVfitMEM/interface/me_ggH_lo_lit.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"

#include "LHAPDF/LHAPDF.h"

#include <TMatrixD.h>

namespace svFitMEM
{
  class SVfitIntegrand_lo
  {
   public:
    /// error codes that can be read out by SVfitMEM class
    enum ErrorCodes {
      None            = 0x00000000,
      MatrixInversion = 0x00000001,
      LeptonNumber    = 0x00000010,
      TestMass        = 0x00000100
    };

    SVfitIntegrand_lo(double, const std::string&, int, const std::string&, int);
    ~SVfitIntegrand_lo();
  
    void setIdxLeg1_X(int idx) { idxLeg1_X_ = idx; }
    void setIdxLeg1_phi(int idx) { idxLeg1_phi_ = idx; }
    void setIdxLeg1VisPtShift(int idx) { idxLeg1VisPtShift_ = idx; }
    void setIdxLeg1_mNuNu(int idx) { idxLeg1_mNuNu_ = idx; }
    void setIdxLeg2_t(int idx) { idxLeg2_t_ = idx; }
    void setIdxLeg2_phi(int idx) { idxLeg2_phi_ = idx; }
    void setIdxLeg2VisPtShift(int idx) { idxLeg2VisPtShift_ = idx; }
    void setIdxLeg2_mNuNu(int idx) { idxLeg2_mNuNu_ = idx; }
    void setNumDimensions(unsigned numDimensions) { numDimensions_ = numDimensions; }

    void setMtest(double);

    /// set transfer functions for pT of hadronic tau decays
    void setHadTauTF(const HadTauTFBase* hadTauTF) 
    { 
      delete hadTauTF1_;
      hadTauTF1_ = hadTauTF->Clone("leg1"); 
      delete hadTauTF2_;
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

    /// set correlation between hadronic tau pT and MET
    void setRhoHadTau(double rhoHadTau) 
    { 
      rhoHadTau_ = rhoHadTau;
    }

    /// set momenta of visible tau decay products and of reconstructed missing transverse energy
    void setInputs(const std::vector<svFitMEM::MeasuredTauLepton>&, double, double, const TMatrixD&);

    /// evaluate integrand for given value of integration variables x
    double Eval(const double* x) const;

    /// static pointer to this (needed for interfacing the likelihood function calls to VEGAS integration)
    static const SVfitIntegrand_lo* gSVfitIntegrand;

    enum { kMadgraph, kLiterature };

   protected:
    /// measured tau leptons
    MeasuredTauLepton measuredTauLepton1_;
    bool leg1isLep_;
    double leg1Mass_;
    double leg1Mass2_;
    double leg1eX_x_;
    double leg1eX_y_;
    double leg1eX_z_;
    double leg1eY_x_;
    double leg1eY_y_;
    double leg1eY_z_;
    double leg1eZ_x_;
    double leg1eZ_y_;
    double leg1eZ_z_;
    MeasuredTauLepton measuredTauLepton2_;
    bool leg2isLep_;
    double leg2Mass_;
    double leg2Mass2_;
    double leg2eX_x_;
    double leg2eX_y_;
    double leg2eX_z_;
    double leg2eY_x_;
    double leg2eY_y_;
    double leg2eY_z_;
    double leg2eZ_x_;
    double leg2eZ_y_;
    double leg2eZ_z_;

    int mode_;

    mutable double mVis_measured_;
    mutable double mVis2_measured_;
    double mTest_;
    double mTest2_;
    mutable double GammaH_;
    mutable double q2_;
    mutable double GammaH_times_mTest_;
    mutable double GammaH2_times_mTest2_;

    double sqrtS_;
    double s_;
    double invSqrtS_;
    Vector beamAxis_;

    /// measured MET
    double measuredMETx_;
    double measuredMETy_;

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

    double rhoHadTau_;

    int idxLeg1_X_;
    int idxLeg1_phi_;
    int idxLeg1VisPtShift_;
    int idxLeg1_mNuNu_;
    int idxLeg2_t_;
    int idxLeg2_phi_;
    int idxLeg2VisPtShift_;
    int idxLeg2_mNuNu_;
    unsigned numDimensions_;

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

    /// error code that can be passed on
    int errorCode_;

    /// verbosity level
    int verbosity_;
  };
}

#endif
