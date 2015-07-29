#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandStableTaus.h"

#include <TMath.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

/// global function pointer, needed for VEGAS integration
const HttXsectionIntegrandStableTaus* HttXsectionIntegrandStableTaus::gHttXsectionIntegrandStableTaus = 0;

HttXsectionIntegrandStableTaus::HttXsectionIntegrandStableTaus(double sqrtS, double mH, const std::string& pdfName, int mode, const std::string& madgraphFileName, int verbosity) 
  : mode_(mode),
    mH_(mH),
    mH2_(mH*mH),
    GammaH_(1.e-2*mH_),
    sqrtS_(sqrtS),
    s_(square(sqrtS_)),
    invSqrtS_(1./sqrtS_),
    beamAxis_(0., 0., 1.),
    pdf_(0),
    pdfIsInitialized_(false),
    me_madgraph_(false, true),
    me_madgraph_isInitialized_(false),
    me_lit_(pdfName, false, true),
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandStableTaus::HttXsectionIntegrandStableTaus>:" << std::endl;
  }

  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    pdf_ = LHAPDF::mkPDF(pdfName.data(), 0);
    pdfIsInitialized_ = true;
  }

  // initialize Madgraph
  if ( madgraphFileName != "" ) {
    me_madgraph_.initProc(madgraphFileName);
    me_madgraph_isInitialized_ = true;
  } else if ( mode == kMadgraph ) {
    std::cerr << "Error in <HttXsectionIntegrandStableTaus>: No param.dat file for Madgraph given !!" << std::endl;
    assert(0);
  }
  if ( mode_ == kMadgraph ) {
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
  gHttXsectionIntegrandStableTaus = this;
}

HttXsectionIntegrandStableTaus::~HttXsectionIntegrandStableTaus()
{
  std::cout << "<HttXsectionIntegrandStableTaus::~HttXsectionIntegrandStableTaus>:" << std::endl;
  delete pdf_;

  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphTau1P4_;
  delete [] madgraphTau2P4_;
}

namespace
{
  double compDgDp2z(double tau1Pt2, double tau1Pz, double tau2Pz)
  {
    double tau1Pz2 = square(tau1Pz);
    double tau2Pz2 = square(tau2Pz);
    if ( tau2Pz == 0. ) return 0.;
    else return -2.*tau1Pz + 2.*TMath::Sqrt((tau1Pt2 + tau1Pz2)/(tau1Pt2 + tau2Pz2))*tau2Pz;
  }
}

double
HttXsectionIntegrandStableTaus::Eval(const double* x) const 
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<HttXsectionIntegrandStableTaus:Eval(const double*)>:" << std::endl;
  }

  double tau1Px = x[0];
  double tau1Py = x[1];
  double tau1Pz = x[2];
  double tau1Px2 = square(tau1Px);
  double tau1Py2 = square(tau1Py);
  double tau1Pt2 = tau1Px2 + tau1Py2;
  //if ( tau1Pt2 > (0.25*mH2_) ) return 0.;
  double tau1Pz2 = square(tau1Pz);
  double tau1P2 = tau1Px2 + tau1Py2 + tau1Pz2;
  double tau1P = TMath::Sqrt(tau1P2);
  double tk = x[3];

  if ( me_madgraph_isInitialized_ && TMath::Abs(me_madgraph_.getHiggsMass() - mH_) > 1.e-3*mH_ ) {
    std::cerr << "Error: Higgs mass defined in Madgraph = " << me_madgraph_.getHiggsMass() << " does not match mH = " << mH_ << " !!" << std::endl;
    assert(0);
  }

  double q2 = mH2_ + mH_*GammaH_*TMath::Tan(tk);
  if ( q2 <= 0. ) return 0.;
  double jacobiFactor = (square(q2 - mH2_) + mH2_*square(GammaH_))/(mH_*GammaH_); // dq2/dtk, taken from Eq. (8) in arXiv:1010.2263v3, with Pi factor removed from denominator after checking with Mathematica

  double tau2Px = -tau1Px;
  double tau2Py = -tau1Py; 
  double term1 = q2 - 2.*tau1Pt2;
  double term22 = q2*(q2 - 4.*tau1Pt2);
  if ( term22 <= 0. || tau1Pt2 <= 0. ) return 0.;
  double term2 = TMath::Sqrt(term22);
  double tau2Pz_p = (tau1Pz*term1 + tau1P*term2)/(2.*tau1Pt2);
  double absDgDp2z_p = TMath::Abs(compDgDp2z(tau1Pt2, tau1Pz, tau2Pz_p));
  double tau2Pz_m = (tau1Pz*term1 - tau1P*term2)/(2.*tau1Pt2);
  double absDgDp2z_m = TMath::Abs(compDgDp2z(tau1Pt2, tau1Pz, tau2Pz_m));
  double prob = 0.;
  if ( TMath::Abs(tau2Pz_p) <= sqrtS_ && absDgDp2z_p != 0. ) { 
    double prob_p = compProb(tau1Px, tau1Py, tau1Pz, tau2Px, tau2Py, tau2Pz_p, jacobiFactor/absDgDp2z_p);
    if ( verbosity_ >= 3 ) {
      std::cout << "solution+: tau2Pz = " << tau2Pz_p << ", prob = " << prob_p << std::endl;
    }
    prob += prob_p;
  }
  if ( TMath::Abs(tau2Pz_m) <= sqrtS_ && absDgDp2z_m != 0. ) { 
    double prob_m = compProb(tau1Px, tau1Py, tau1Pz, tau2Px, tau2Py, tau2Pz_m, jacobiFactor/absDgDp2z_m);
    if ( verbosity_ >= 3 ) {
      std::cout << "solution-: tau2Pz = " << tau2Pz_m << ", prob = " << prob_m << std::endl;
    }
    prob += prob_m;
  }

  return prob;
}

double
HttXsectionIntegrandStableTaus::compProb(double tau1Px, double tau1Py, double tau1Pz, double tau2Px, double tau2Py, double tau2Pz, double jacobiFactor) const
{
  double tau1En = TMath::Sqrt(tau1Px*tau1Px + tau1Py*tau1Py + tau1Pz*tau1Pz + tauLeptonMass2);
  double tau2En = TMath::Sqrt(tau2Px*tau2Px + tau2Py*tau2Py + tau2Pz*tau2Pz + tauLeptonMass2);

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  double ditauPx = tau1Px + tau2Px;
  double ditauPy = tau1Py + tau2Py;
  double ditauPz = tau1Pz + tau2Pz;
  double ditauEn = tau1En + tau2En;
  double ditauMass = TMath::Sqrt(square(ditauEn) - (square(ditauPx) + square(ditauPy) + square(ditauPz)));
  //if ( !(ditauMass > 0.70*mH_ && ditauMass < 1.30*mH_) ) return 0.;
  double xa = invSqrtS_*(ditauEn + ditauPz);
  double xb = invSqrtS_*(ditauEn - ditauPz);   
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  //double Q = mH_;
  double Q = ditauMass;
  assert(pdfIsInitialized_);
  double fa = pdf_->xfxQ(21, xa, Q)/xa; // gluon distribution
  double fb = pdf_->xfxQ(21, xb, Q)/xb;
  double prob_PDF = (fa*fb);

  // evaluate flux factor
  double prob_flux = (1./(2.*s_*xa*xb));  

  // perform boost into MEM frame and evaluate LO matrix element, 
  // computed by Madgraph
  double memFramePx = tau1Px + tau2Px;
  double memFramePy = tau1Py + tau2Py;
  double memFrameEn = tau1En + tau2En;
  Vector boost(-memFramePx/memFrameEn, -memFramePy/memFrameEn, 0.);
  //std::cout << "boost: Px = " << boost.x() << ", Py = " << boost.y() << ", Pz = " << boost.z() << std::endl;
  LorentzVector tau1P4(tau1Px, tau1Py, tau1Pz, tau1En);
  LorentzVector tau1P4_mem = ROOT::Math::VectorUtil::boost(tau1P4, boost); 
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);
  LorentzVector tau2P4_mem = ROOT::Math::VectorUtil::boost(tau2P4, boost); 
  if ( verbosity_ >= 3 ) {
    std::cout << "lab:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4.pt() << ", eta = " << tau1P4.eta() << ", phi = " << tau1P4.phi() << ", mass = " << tau1P4.mass() << std::endl;
    std::cout << "      (En = " << tau1P4.energy() << ", Px = " << tau1P4.px() << ", Py = " << tau1P4.py() << ", Pz = " << tau1P4.pz() << ")" << std::endl;    
    std::cout << " tau2: Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << ", mass = " << tau2P4.mass() << std::endl;
    std::cout << "      (En = " << tau2P4.energy() << ", Px = " << tau2P4.px() << ", Py = " << tau2P4.py() << ", Pz = " << tau2P4.pz() << ")" << std::endl;  
    LorentzVector ditauP4 = tau1P4 + tau2P4;
    std::cout << " ditau: Pt = " << ditauP4.pt() << ", eta = " << ditauP4.eta() << ", phi = " << ditauP4.phi() << ", mass = " << ditauP4.mass() << std::endl;
    std::cout << "mem:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4_mem.pt() << ", eta = " << tau1P4_mem.eta() << ", phi = " << tau1P4_mem.phi() << ", mass = " << tau1P4_mem.mass() << std::endl;
    std::cout << "      (En = " << tau1P4_mem.energy() << ", Px = " << tau1P4_mem.px() << ", Py = " << tau1P4_mem.py() << ", Pz = " << tau1P4_mem.pz() << ")" << std::endl;  
    std::cout << " tau2: Pt = " << tau2P4_mem.pt() << ", eta = " << tau2P4_mem.eta() << ", phi = " << tau2P4_mem.phi() << ", mass = " << tau2P4_mem.mass() << std::endl;
    std::cout << "      (En = " << tau2P4_mem.energy() << ", Px = " << tau2P4_mem.px() << ", Py = " << tau2P4_mem.py() << ", Pz = " << tau2P4_mem.pz() << ")" << std::endl;  
    LorentzVector ditauP4_mem = tau1P4_mem + tau2P4_mem;
    std::cout << " ditau: Pt = " << ditauP4_mem.pt() << ", eta = " << ditauP4_mem.eta() << ", phi = " << ditauP4_mem.phi() << ", mass = " << ditauP4_mem.mass() << std::endl;
  }
  madgraphGluon1P4_[0] =  0.5*xa*sqrtS_; 
  madgraphGluon1P4_[3] = +0.5*xa*sqrtS_;
  madgraphGluon2P4_[0] =  0.5*xb*sqrtS_;
  madgraphGluon2P4_[3] = -0.5*xb*sqrtS_;
  madgraphTau1P4_[0] = tau1P4_mem.energy();
  madgraphTau1P4_[1] = tau1P4_mem.px();
  madgraphTau1P4_[2] = tau1P4_mem.py();
  madgraphTau1P4_[3] = tau1P4_mem.pz();
  madgraphTau2P4_[0] = tau2P4_mem.energy();
  madgraphTau2P4_[1] = tau2P4_mem.px();
  madgraphTau2P4_[2] = tau2P4_mem.py();
  madgraphTau2P4_[3] = tau2P4_mem.pz();
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
  
  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  const double constFactor = conversionFactor/(8.*square(TMath::Pi()));
  double prob_PS = constFactor/(s_*tau1En*tau2En);

  double prob = prob_flux*prob_PDF*prob_ME*prob_PS*jacobiFactor;
  assert(prob >= 0.);
  //if ( verbosity_ >= 1 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS = " << prob_PS << ", Jacobi = " << jacobiFactor
	      << " --> returning " << prob << std::endl;
  //}
  return prob;
}
