#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandStableTaus.h"

#include <TMath.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  double xfx(int nset, double x, double Q, int fl);
}

/// global function pointer, needed for VEGAS integration
const HttXsectionIntegrandStableTaus* HttXsectionIntegrandStableTaus::gHttXsectionIntegrandStableTaus = 0;
bool HttXsectionIntegrandStableTaus::pdfIsInitialized_ = false;

HttXsectionIntegrandStableTaus::HttXsectionIntegrandStableTaus(const std::string& madgraphFileName, double sqrtS, double mH, const std::string& pdfFileName, int verbosity) 
  : mH_(mH),
    mH2_(mH*mH),
    GammaH_(1.e-2*mH_),
    sqrtS_(sqrtS),
    s_(square(sqrtS_)),
    invSqrtS_(1./sqrtS_),
    beamAxis_(0., 0., 1.),
    me_madgraph_(false),
    me_lit_(false),
    verbosity_(verbosity)
{
  //if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandStableTaus::HttXsectionIntegrandStableTaus>:" << std::endl;
  //}

  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    LHAPDF::initPDFSet(1, pdfFileName);
    pdfIsInitialized_ = true;
  }

  // initialize Madgraph
  me_madgraph_.initProc(madgraphFileName);
  //GammaH_ = me_madgraph_.getHiggsWidth();
  //std::cout << "mH = " << mH_ << ": GammaH = " << GammaH_ << std::endl;

  me_lit_.setS(s_);
  me_lit_.setHiggsMass(mH_);
  me_lit_.setHiggsWidth(GammaH_);
  me_lit_.setBR(1.);
  //const double v2 = square(246.22); // GeV^2
  //me_lit_.setBR(2.*tauLeptonMass2/v2)*square(mH_)*(1. - 4.*tauLeptonMass2/square(mH_))); // square of Higgs-tau Yukawa coupling

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

  //GammaH_ = me_madgraph_.getHiggsWidth();
  if ( TMath::Abs(me_madgraph_.getHiggsMass() - mH_) > 1.e-3*mH_ ) {
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
    if ( verbosity_ >= 2 ) {
      std::cout << "solution+: tau2Pz = " << tau2Pz_p << ", prob = " << prob_p << std::endl;
    }
    prob += prob_p;
  }
  if ( TMath::Abs(tau2Pz_m) <= sqrtS_ && absDgDp2z_m != 0. ) { 
    double prob_m = compProb(tau1Px, tau1Py, tau1Pz, tau2Px, tau2Py, tau2Pz_m, jacobiFactor/absDgDp2z_m);
    if ( verbosity_ >= 2 ) {
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
  double fa = LHAPDF::xfx(1, xa, Q, 0)/xa;
  double fb = LHAPDF::xfx(1, xb, Q, 0)/xb;
  double prob_PDF = (fa*fb);

  // evaluate flux factor
  double prob_flux = (1./(2.*s_*xa*xb));  
  
  // perform boost into MEM frame and evaluate LO matrix element, 
  // computed by Madgraph
  LorentzVector tau1P4(tau1Px, tau1Py, tau1Pz, tau1En);
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);
  LorentzVector comSystem = tau1P4 + tau2P4;
  Vector boost = comSystem.BoostToCM();
  LorentzVector tau1P4_rf = ROOT::Math::VectorUtil::boost(tau1P4, boost); 
  LorentzVector tau2P4_rf = ROOT::Math::VectorUtil::boost(tau2P4, boost); 
  if ( verbosity_ >= 2 ) {
    std::cout << "lab:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4.pt() << ", eta = " << tau1P4.eta() << ", phi = " << tau1P4.phi() << ", mass = " << tau1P4.mass() << std::endl;
    std::cout << "      (En = " << tau1P4.energy() << ", Px = " << tau1P4.px() << ", Py = " << tau1P4.py() << ", Pz = " << tau1P4.pz() << ")" << std::endl;    
    std::cout << " tau2: Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << ", mass = " << tau2P4.mass() << std::endl;
    std::cout << "      (En = " << tau2P4.energy() << ", Px = " << tau2P4.px() << ", Py = " << tau2P4.py() << ", Pz = " << tau2P4.pz() << ")" << std::endl;  
    LorentzVector ditauP4 = tau1P4 + tau2P4;
    std::cout << " ditau: Pt = " << ditauP4.pt() << ", eta = " << ditauP4.eta() << ", phi = " << ditauP4.phi() << ", mass = " << ditauP4.mass() << std::endl;
    std::cout << "rf:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4_rf.pt() << ", eta = " << tau1P4_rf.eta() << ", phi = " << tau1P4_rf.phi() << ", mass = " << tau1P4_rf.mass() << std::endl;
    std::cout << "      (En = " << tau1P4_rf.energy() << ", Px = " << tau1P4_rf.px() << ", Py = " << tau1P4_rf.py() << ", Pz = " << tau1P4_rf.pz() << ")" << std::endl;  
    std::cout << " tau2: Pt = " << tau2P4_rf.pt() << ", eta = " << tau2P4_rf.eta() << ", phi = " << tau2P4_rf.phi() << ", mass = " << tau2P4_rf.mass() << std::endl;
    std::cout << "      (En = " << tau2P4_rf.energy() << ", Px = " << tau2P4_rf.px() << ", Py = " << tau2P4_rf.py() << ", Pz = " << tau2P4_rf.pz() << ")" << std::endl;  
    LorentzVector ditauP4_rf = tau1P4_rf + tau2P4_rf;
    std::cout << " ditau: Pt = " << ditauP4_rf.pt() << ", eta = " << ditauP4_rf.eta() << ", phi = " << ditauP4_rf.phi() << ", mass = " << ditauP4_rf.mass() << std::endl;
  }
  madgraphGluon1P4_[0] =  0.5*Q; 
  madgraphGluon1P4_[3] = +0.5*Q;
  madgraphGluon2P4_[0] =  0.5*Q;
  madgraphGluon2P4_[3] = -0.5*Q;
  madgraphTau1P4_[0] = tau1P4_rf.energy();
  madgraphTau1P4_[1] = tau1P4_rf.px();
  madgraphTau1P4_[2] = tau1P4_rf.py();
  madgraphTau1P4_[3] = tau1P4_rf.pz();
  madgraphTau2P4_[0] = tau2P4_rf.energy();
  madgraphTau2P4_[1] = tau2P4_rf.px();
  madgraphTau2P4_[2] = tau2P4_rf.py();
  madgraphTau2P4_[3] = tau2P4_rf.pz();
  me_madgraph_.setMomenta(madgraphMomenta_);
  me_madgraph_.sigmaKin();
  //double prob_ME = me_madgraph_.getMatrixElements()[0];
  me_lit_.setHiggsWidth(GammaH_);
  me_lit_.setMomenta(madgraphMomenta_);
  double prob_ME = me_lit_.getMatrixElement();
  //std::cout << "prob_ME: madgraph = " << me_madgraph_.getMatrixElements()[0] << ", lit = " << me_lit_.getMatrixElement() << std::endl;
  assert(prob_ME >= 0.);

  const double hbar_c = 0.1973; // GeV fm
  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  const double constFactor = conversionFactor/(8.*square(TMath::Pi()));
  double prob_PS = constFactor/(s_*tau1En*tau2En);

  double prob = prob_flux*prob_PDF*prob_ME*prob_PS*jacobiFactor;
  assert(prob >= 0.);
  if ( verbosity_ >= 1 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS = " << prob_PS << ", Jacobi = " << jacobiFactor
	      << " --> returning " << prob << std::endl;
  }
  return prob;
}
