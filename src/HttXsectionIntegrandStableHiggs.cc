#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandStableHiggs.h"

#include <TMath.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  double xfx(int nset, double x, double Q, int fl);
}

/// global function pointer, needed for VEGAS integration
const HttXsectionIntegrandStableHiggs* HttXsectionIntegrandStableHiggs::gHttXsectionIntegrandStableHiggs = 0;
bool HttXsectionIntegrandStableHiggs::pdfIsInitialized_ = false;

HttXsectionIntegrandStableHiggs::HttXsectionIntegrandStableHiggs(double sqrtS, double mH, const std::string& pdfFileName, bool applyNWA, int mode, const std::string& madgraphFileName, int verbosity) 
  : mode_(mode),
    mH_(mH),
    mH2_(mH*mH), 
    GammaH_(1.e-2*mH_),
    sqrtS_(sqrtS),
    s_(square(sqrtS_)),
    invSqrtS_(1./sqrtS_),
    applyNWA_(applyNWA),
    beamAxis_(0., 0., 1.),
    me_madgraph_(applyNWA_, false),
    me_madgraph_isInitialized_(false),
    me_lit_(applyNWA_, false),
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandStableHiggs::HttXsectionIntegrandStableHiggs>:" << std::endl;
  }

  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    LHAPDF::initPDFSet(1, pdfFileName);
    pdfIsInitialized_ = true;
  }

   // initialize Madgraph
  if ( madgraphFileName != "" ) {
    me_madgraph_.initProc(madgraphFileName);
    me_madgraph_isInitialized_ = true;
  } else if ( mode_ == kMadgraph ) {
    std::cerr << "Error in <HttXsectionIntegrandStableTaus>: No param.dat file for Madgraph given !!" << std::endl;
    assert(0);
  }
  if ( mode == kMadgraph ) {
    GammaH_ = me_madgraph_.getHiggsWidth();
  }

  me_lit_.setS(s_);
  me_lit_.setHiggsMass(mH_);
  me_lit_.setHiggsWidth(GammaH_);

  madgraphGluon1P4_ = new double[4];
  madgraphGluon1P4_[1] = 0.;
  madgraphGluon1P4_[2] = 0.;
  madgraphMomenta_.push_back(madgraphGluon1P4_);
  madgraphGluon2P4_ = new double[4];
  madgraphGluon2P4_[1] = 0.;
  madgraphGluon2P4_[2] = 0.;
  madgraphMomenta_.push_back(madgraphGluon2P4_);
  madgraphTau1P4_ = new double[4];
  madgraphMomenta_.push_back(madgraphTau1P4_);
  madgraphTau2P4_ = new double[4];
  madgraphMomenta_.push_back(madgraphTau2P4_);

  // set global function pointer to this
  gHttXsectionIntegrandStableHiggs = this;
}

HttXsectionIntegrandStableHiggs::~HttXsectionIntegrandStableHiggs()
{
  std::cout << "<HttXsectionIntegrandStableHiggs::~HttXsectionIntegrandStableHiggs>:" << std::endl;
  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphTau1P4_;
  delete [] madgraphTau2P4_;
}

double
HttXsectionIntegrandStableHiggs::Eval(const double* x) const 
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<HttXsectionIntegrandStableTaus:Eval(const double*)>:" << std::endl;
  }
  
  double tk = ( applyNWA_ ) ? 0. : x[1];

  if ( me_madgraph_isInitialized_ && TMath::Abs(me_madgraph_.getHiggsMass() - mH_) > 1.e-3*mH_ ) {
    std::cerr << "Error: Higgs mass defined in Madgraph = " << me_madgraph_.getHiggsMass() << " does not match mH = " << mH_ << " !!" << std::endl;
    assert(0);
  }
  double q2 = mH2_ + mH_*GammaH_*TMath::Tan(tk);
  if ( q2 <= 0. ) return 0.;
  double q = TMath::Sqrt(q2);
  double jacobiFactor = (square(q2 - mH2_) + mH2_*square(GammaH_))/(mH_*GammaH_); // dq2/dtk, taken from Eq. (8) in arXiv:1010.2263v3, with Pi factor removed from denominator after checking with Mathematica

  double tauP  = 0.5*q;
  double tauEn = TMath::Sqrt(square(tauP) + tauLeptonMass2);
  double sHat = square(2.*tauEn);
    
  // compute Bjorken-x of incoming protons and evaluate PDF factor
  //if ( !(q > (0.70*mH_) && q < (1.30e+3*mH_)) ) return 0.;
  double xa = x[0];
  double xb = sHat/(s_*xa);
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  //double Q = mH_;
  double Q = sqrt(sHat);
  assert(pdfIsInitialized_);
  double fa = LHAPDF::xfx(1, xa, Q, 0)/xa;
  double fb = LHAPDF::xfx(1, xb, Q, 0)/xb;
  double prob_PDF = (fa*fb);

  // perform boost into MEM frame and evaluate LO matrix element, 
  // computed by Madgraph  
  madgraphGluon1P4_[0] =  0.5*tauEn; 
  madgraphGluon1P4_[3] = +0.5*tauEn;
  madgraphGluon2P4_[0] =  0.5*tauEn;
  madgraphGluon2P4_[3] = -0.5*tauEn;
  madgraphTau1P4_[0] = tauEn;
  madgraphTau1P4_[1] = 0.;
  madgraphTau1P4_[2] = 0.;
  madgraphTau1P4_[3] = +tauP;
  madgraphTau2P4_[0] = tauEn;
  madgraphTau2P4_[1] = 0.;
  madgraphTau2P4_[2] = 0.;
  madgraphTau2P4_[3] = -tauP;
  double prob_ME_madgraph = -1.;
  if ( me_madgraph_isInitialized_ ) {
    me_madgraph_.setMomenta(madgraphMomenta_);
    me_madgraph_.sigmaKin();
    prob_ME_madgraph = me_madgraph_.getMatrixElements()[0];
  }
  me_lit_.setMomenta(madgraphMomenta_);
  double prob_ME_lit = me_lit_.getMatrixElement();
  if ( verbosity_ >= 2 ) {
    std::cout << "prob_ME: madgraph = " << prob_ME_madgraph << ", lit = " << prob_ME_lit << std::endl;
  }
  double prob_ME = 0.;
  if ( mode_ == kMadgraph ) {
    prob_ME = prob_ME_madgraph;
  } else if ( mode_ == kLiterature ) {
    prob_ME = prob_ME_lit;
  } else {
    assert(0);
  }
  assert(prob_ME >= 0.);

  // compute gg -> Higgs, Higgs -> tautau cross-section
  // according to Eq. (43) in http://www.itp.phys.ethz.ch/education/fs10/aft/Thesis_MB.pdf

  const double hbar_c = 0.1973; // GeV fm
  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  double constFactor = conversionFactor*TMath::Pi();
  double prob_const = constFactor/(sHat*s_*xa);

  double prob = prob_PDF*prob_ME*prob_const;
  if ( !applyNWA_ ) {
    prob *= jacobiFactor;
  }
  assert(prob >= 0.);
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: PDF = " << prob_PDF << ", ME = " << prob_ME << ", const = " << prob_const;
    if ( !applyNWA_ ) {
      std::cout << ", Jacobi = " << jacobiFactor;
    }
    std::cout << " --> returning " << prob << std::endl;
  }
  return prob;
}

