#include "TauAnalysis/SVfitMEM/interface/SVfitIntegrand.h"

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVEGAS.h"

#include <TMath.h>
#include <TH1.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  double xfx(int nset, double x, double Q, int fl);
}

/// global function pointer, needed for VEGAS integration
const SVfitIntegrand* SVfitIntegrand::gSVfitIntegrand = 0;
bool SVfitIntegrand::pdfIsInitialized_ = false;

SVfitIntegrand::SVfitIntegrand(double sqrtS, const std::string& pdfFileName, int mode, const std::string& madgraphFileName, int verbosity) 
  : mode_(mode),
    mTest_(0.),
    mTest2_(0.),
    GammaH_(1.e-2*mTest_),
    GammaH_times_mTest_(GammaH_*mTest_),
    GammaH2_times_mTest2_(square(GammaH_times_mTest_)),
    sqrtS_(sqrtS),
    s_(square(sqrtS_)),
    invSqrtS_(1./sqrtS_),
    beamAxis_(0., 0., 1.),
    invCovMET_(2,2),
    shiftVisPt_(false),
    leg1lutVisPtRes_(0),
    leg2lutVisPtRes_(0),
    idxLeg1_X_(-1),
    idxLeg1_phi_(-1),
    idxLeg1VisPtShift_(-1),
    idxLeg1_mNuNu_(-1),
    idxLeg2_t_(-1),
    idxLeg2_phi_(-1),
    idxLeg2VisPtShift_(-1),
    idxLeg2_mNuNu_(-1),
    me_madgraph_(false, true),
    me_madgraph_isInitialized_(false),
    me_lit_(false, true),
    errorCode_(0),
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<SVfitIntegrand::SVfitIntegrand>:" << std::endl;
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
  } else if ( mode == kMadgraph ) {
    std::cerr << "Error in <SVfitIntegrand>: No param.dat file for Madgraph given !!" << std::endl;
    assert(0);
  }
  if ( mode_ == kMadgraph ) {
    GammaH_ = me_madgraph_.getHiggsWidth();
  }

  me_lit_.setS(s_);
  // set Higgs -> tautau decay branching fraction 
  // to common value for all test-mass hypotheses
  me_lit_.setBR(1.e-1); 

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
  gSVfitIntegrand = this;
}

SVfitIntegrand::~SVfitIntegrand()
{
  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphTau1P4_;
  delete [] madgraphTau2P4_;
}

void 
SVfitIntegrand::setMtest(double mTest) 
{ 
  // reset 'TestMass' error code
  errorCode_ &= (errorCode_ ^ TestMass);
  
  if ( mTest < TMath::Sqrt(mVis2_) ) {
    std::cerr << "Error: Cannot have mTest < mVis = " << TMath::Sqrt(mVis2_) << " !!" << std::endl;
    errorCode_ |= TestMass;
    return;
  }
  mTest_ = mTest; 
  mTest2_ = square(mTest); 
  GammaH_ = 1.e-2*mTest_;
  GammaH_times_mTest_ = GammaH_*mTest_;
  GammaH2_times_mTest2_ = square(GammaH_times_mTest_);
}

void 
SVfitIntegrand::shiftVisPt(bool value, const TH1* leg1lutVisPtRes, const TH1* leg2lutVisPtRes)
{
  shiftVisPt_ = value;
  if ( shiftVisPt_ ) {
    leg1lutVisPtRes_ = leg1lutVisPtRes;
    leg2lutVisPtRes_ = leg2lutVisPtRes;
  }
}

namespace
{
  double norm(const Vector& v)
  {
    return TMath::Sqrt(v.mag2());
  }
}

void 
SVfitIntegrand::setInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons, double measuredMETx, double measuredMETy, const TMatrixD& covMET) 
{
  if ( verbosity_ ) {
    std::cout << "<SVfitIntegrand::setInputs>:" << std::endl;
  }

  // reset 'LeptonNumber' and 'MatrixInversion' error codes
  errorCode_ &= (errorCode_ ^ LeptonNumber);
  errorCode_ &= (errorCode_ ^ MatrixInversion);

  if ( measuredTauLeptons.size() != 2 ) {
    std::cerr << "Error: Number of MeasuredTauLeptons is not equal to two !!" << std::endl;
    errorCode_ |= LeptonNumber;
  }
  measuredTauLepton1_ = measuredTauLeptons[0];
  leg1isLep_ = measuredTauLepton1_.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton1_.type() == MeasuredTauLepton::kTauToMuDecay;
  leg1Mass_ = measuredTauLepton1_.mass();
  leg1Mass2_ = square(leg1Mass_);
  Vector eZ1 = normalize(measuredTauLepton1_.p3());
  Vector eY1 = normalize(compCrossProduct(eZ1, beamAxis_));
  Vector eX1 = normalize(compCrossProduct(eY1, eZ1));
  if ( verbosity_ >= 2 ) {
    std::cout << "eX1: theta = " << eX1.theta() << ", phi = " << eX1.phi() << ", norm = " << norm(eX1) << std::endl;
    std::cout << "eY1: theta = " << eY1.theta() << ", phi = " << eY1.phi() << ", norm = " << norm(eY1) << std::endl;
    std::cout << "eZ1: theta = " << eZ1.theta() << ", phi = " << eZ1.phi() << ", norm = " << norm(eZ1) << std::endl;
    std::cout << "(eX1 x eY1 = " << norm(compCrossProduct(eX1, eY1)) << ", eX1 x eZ1 = " << norm(compCrossProduct(eY1, eZ1)) << ", eY1 x eZ1 = " << norm(compCrossProduct(eY1, eZ1)) << ")" << std::endl;
  }
  leg1eX_x_ = eX1.x();
  leg1eX_y_ = eX1.y();
  leg1eX_z_ = eX1.z();
  leg1eY_x_ = eY1.x();
  leg1eY_y_ = eY1.y();
  leg1eY_z_ = eY1.z();
  leg1eZ_x_ = eZ1.x();
  leg1eZ_y_ = eZ1.y();
  leg1eZ_z_ = eZ1.z();

  measuredTauLepton2_ = measuredTauLeptons[1];
  leg2isLep_ = measuredTauLepton2_.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton2_.type() == MeasuredTauLepton::kTauToMuDecay;
  leg2Mass_ = measuredTauLepton2_.mass();
  leg2Mass2_ = square(leg2Mass_);
  Vector eZ2 = normalize(measuredTauLepton2_.p3());
  Vector eY2 = normalize(compCrossProduct(eZ2, beamAxis_));
  Vector eX2 = normalize(compCrossProduct(eY2, eZ2));
  if ( verbosity_ >= 2 ) {
    std::cout << "eX2: theta = " << eX2.theta() << ", phi = " << eX2.phi() << ", norm = " << norm(eX2) << std::endl;
    std::cout << "eY2: theta = " << eY2.theta() << ", phi = " << eY2.phi() << ", norm = " << norm(eY2) << std::endl;
    std::cout << "eZ2: theta = " << eZ2.theta() << ", phi = " << eZ2.phi() << ", norm = " << norm(eZ2) << std::endl;
    //std::cout << "(eX2 x eY2 = " << norm(compCrossProduct(eX2, eY2)) << ", eX2 x eZ2 = " << norm(compCrossProduct(eY2, eZ2)) << ", eY2 x eZ2 = " << norm(compCrossProduct(eY2, eZ2)) << ")" << std::endl;
  }
  leg2eX_x_ = eX2.x();
  leg2eX_y_ = eX2.y();
  leg2eX_z_ = eX2.z();
  leg2eY_x_ = eY2.x();
  leg2eY_y_ = eY2.y();
  leg2eY_z_ = eY2.z();
  leg2eZ_x_ = eZ2.x();
  leg2eZ_y_ = eZ2.y();
  leg2eZ_z_ = eZ2.z();

  mVis_  = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).mass();
  if ( verbosity_ >= 2 ) {
    std::cout << "mVis = " << mVis_ << std::endl;
  }
  mVis2_ = square(mVis_);

  measuredMETx_ = measuredMETx;
  measuredMETy_ = measuredMETy;

  // determine transfer matrix for MET
  invCovMET_ = covMET;
  double covDet = invCovMET_.Determinant();
  const_MET_ = 0.;
  if ( covDet != 0 ) { 
    invCovMET_.Invert(); 
    invCovMETxx_ = invCovMET_(0,0);
    invCovMETxy_ = invCovMET_(0,1);
    invCovMETyx_ = invCovMET_(1,0);
    invCovMETyy_ = invCovMET_(1,1);
    const_MET_ = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));
  } else{
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!" << std::endl;
    errorCode_ |= MatrixInversion;
  }
}

namespace
{
  double extractProbFromLUT(double x, const TH1* lut)
  {    
    assert(lut);
    int bin = (const_cast<TH1*>(lut))->FindBin(x);
    int numBins = lut->GetNbinsX();
    if ( bin < 1       ) bin = 1;
    if ( bin > numBins ) bin = numBins;
    return lut->GetBinContent(bin);
  }

  double compCosThetaNuNu(double visEn, double visP, double visMass2, double nunuEn, double nunuP, double nunuMass2)
  {
    double cosThetaNuNu = (visEn*nunuEn - 0.5*(tauLeptonMass2 - (visMass2 + nunuMass2)))/(visP*nunuP);
    return cosThetaNuNu;
  }

  double compPSfactor_tauToLepDecay(double x, double visEn, double visP, double visMass, double nunuEn, double nunuP, double nunuMass)
  {
    //std::cout << "<compPSfactor_tauToLepDecay>:" << std::endl;
    //std::cout << " x = " << x << std::endl;
    //std::cout << " visEn = " << visEn << std::endl;
    //std::cout << " visP = " << visP << std::endl;
    //std::cout << " visMass = " << visMass << std::endl;
    //std::cout << " nunuEn = " << nunuEn << std::endl;
    //std::cout << " nunuP = " << nunuP << std::endl;
    //std::cout << " nunuMass = " << nunuMass << std::endl;
    double visMass2 = square(visMass);
    double nunuMass2 = square(nunuMass);
    if ( x >= (visMass2/tauLeptonMass2) && x <= 1. && nunuMass2 < ((1. - x)*tauLeptonMass2) ) { // physical solution
      const double GFfactor = square(GF)/square(TMath::Pi());
      double tauEn_rf = (tauLeptonMass2 + nunuMass2 - visMass2)/(2.*nunuMass);
      double visEn_rf = tauEn_rf - nunuMass;
      if ( !(tauEn_rf >= tauLeptonMass && visEn_rf >= visMass) ) return 0.;
      double I = GFfactor*nunuMass2*(2.*tauEn_rf*visEn_rf - (2./3.)*TMath::Sqrt((square(tauEn_rf) - tauLeptonMass2)*(square(visEn_rf) - visMass2)));
      double cosThetaNuNu = compCosThetaNuNu(visEn, visP, visMass2, nunuEn, nunuP, nunuMass2);
      const double epsilon = 1.e-3;
      if ( !(cosThetaNuNu >= (-1. + epsilon) && cosThetaNuNu <= +1.) ) return 0.;
      double PSfactor = (visEn + nunuEn)*I/(8.*visP*square(x)*TMath::Sqrt(square(visP) + square(nunuP) + 2.*visP*nunuP*cosThetaNuNu + tauLeptonMass2));
      //-------------------------------------------------------------------------
      // CV: fudge factor to reproduce literature value for cross-section times branching fraction
      PSfactor *= 2.;
      //-------------------------------------------------------------------------
      return PSfactor;
    } else {
      return 0.;
    }
  }

  double compPSfactor_tauToHadDecay(double x, double visEn, double visP, double visMass, double nuEn, double nuP)
  {
    //std::cout << "<compPSfactor_tauToHadDecay>:" << std::endl;
    //std::cout << " x = " << x << std::endl;
    //std::cout << " visEn = " << visEn << std::endl;
    //std::cout << " visP = " << visP << std::endl;
    //std::cout << " visMass = " << visMass << std::endl;
    //std::cout << " nuEn = " << nuEn << std::endl;
    //std::cout << " nuP = " << nuP << std::endl;
    double visMass2 = square(visMass);
    if ( x >= (visMass2/tauLeptonMass2) && x <= 1. ) { // physical solution
      double cosThetaNu = compCosThetaNuNu(visEn, visP, visMass2, nuEn, nuP, 0.);
      //std::cout << "cosThetaNu = " << cosThetaNu << std::endl;
      const double epsilon = 1.e-3;
      if ( !(cosThetaNu >= (-1. + epsilon) && cosThetaNu <= +1.) ) return 0.;
      double PSfactor = (visEn + nuEn)/(8.*visP*square(x)*TMath::Sqrt(square(visP) + square(nuP) + 2.*visP*nuP*cosThetaNu + tauLeptonMass2));
      //-------------------------------------------------------------------------
      // CV: multiply by constant matrix element, 
      //     chosen such that the branching fraction of the tau to decay into hadrons is reproduced      
      const double M2 = 16.*TMath::Pi()*cube(tauLeptonMass)*GammaTauToHad/(tauLeptonMass2 - visMass2);
      PSfactor *= M2;
      //-------------------------------------------------------------------------
      return PSfactor;
    } else {
      return 0.;
    }
  }
}

double
SVfitIntegrand::Eval(const double* x) const 
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<SVfitIntegrand::Eval(const double*)>:" << std::endl;
  }

  // in case of initialization errors don't start to do anything
  if ( errorCode_ != 0 ) { 
    return 0.;
  } 

  double visPtShift1 = ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) ? 1./x[idxLeg1VisPtShift_] : 1.;
  double visPtShift2 = ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) ? 1./x[idxLeg2VisPtShift_] : 1.;

  // compute four-vector of visible decay products for first tau
  double vis1Px = visPtShift1*measuredTauLepton1_.px();
  double vis1Py = visPtShift1*measuredTauLepton1_.py();
  double vis1Pz = visPtShift1*measuredTauLepton1_.pz();
  double vis1En = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);
  double vis1P = vis1P4.P();
  double vis1Theta = vis1P4.theta();
  double vis1Phi = vis1P4.phi();

  // compute four-vector of visible decay products for second tau
  double vis2Px = visPtShift2*measuredTauLepton2_.px();
  double vis2Py = visPtShift2*measuredTauLepton2_.py();
  double vis2Pz = visPtShift2*measuredTauLepton2_.pz();
  double vis2En = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);
  double vis2P = vis2P4.P();
  double vis2Theta = vis2P4.theta();
  double vis2Phi = vis2P4.phi();

  // compute visible energy fractions for both taus
  assert(idxLeg1_X_ != -1);
  double x1 = x[idxLeg1_X_];
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;
  
  assert(idxLeg2_t_ != -1);  
  double tk = x[idxLeg2_t_];
  q2_ = mTest2_ + GammaH_times_mTest_*TMath::Tan(tk);
  if ( q2_ <= 0. ) return 0.;

  mVis_ = (vis1P4 + vis2P4).mass();
  mVis2_ = square(mVis_);
  //if ( mVis2_ < 1.e-3*mTest2_ ) return 0.;
  double x2 = mVis2_/(q2_*x1);
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;

  // compute neutrino and tau lepton four-vector for first tau
  double nu1En = vis1En*(1. - x1)/x1;
  double nu1Mass = ( idxLeg1_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg1_mNuNu_]) : 0.;
  double nu1P = TMath::Sqrt(TMath::Max(0., nu1En*nu1En - square(nu1Mass)));
  assert(idxLeg1_phi_ != -1);  
  double phiNu1 = x[idxLeg1_phi_];
  double cosThetaNu1 = compCosThetaNuNu(vis1En, vis1P, leg1Mass2_, nu1En, nu1P, square(nu1Mass));
  if ( !(cosThetaNu1 >= -1. && cosThetaNu1 <= +1.) ) return 0.;
  double thetaNu1 = TMath::ACos(cosThetaNu1);
  double nu1Px_local = nu1P*TMath::Cos(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Py_local = nu1P*TMath::Sin(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Pz_local = nu1P*TMath::Cos(thetaNu1);
  double nu1Px = nu1Px_local*leg1eX_x_ + nu1Py_local*leg1eY_x_ + nu1Pz_local*leg1eZ_x_;
  double nu1Py = nu1Px_local*leg1eX_y_ + nu1Py_local*leg1eY_y_ + nu1Pz_local*leg1eZ_y_;
  double nu1Pz = nu1Px_local*leg1eX_z_ + nu1Py_local*leg1eY_z_ + nu1Pz_local*leg1eZ_z_;
  LorentzVector nu1P4(nu1Px, nu1Py, nu1Pz, nu1En);

  double tau1En = vis1En + nu1En;
  double tau1Px = vis1P4.px() + nu1Px;
  double tau1Py = vis1P4.py() + nu1Py;
  double tau1Pz = vis1P4.pz() + nu1Pz;
  LorentzVector tau1P4(tau1Px, tau1Py, tau1Pz, tau1En);

  // compute neutrino and tau lepton four-vector for second tau
  double nu2En = vis2En*(1. - x2)/x2;
  double nu2Mass = ( idxLeg2_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg2_mNuNu_]) : 0.;
  double nu2P = TMath::Sqrt(TMath::Max(0., nu2En*nu2En - square(nu2Mass)));
  assert(idxLeg2_phi_ != -2);  
  double phiNu2 = x[idxLeg2_phi_];
  double cosThetaNu2 = compCosThetaNuNu(vis2En, vis2P, leg2Mass2_, nu2En, nu2P, square(nu2Mass));
  if ( !(cosThetaNu2 >= -1. && cosThetaNu2 <= +1.) ) return 0.;
  double thetaNu2 = TMath::ACos(cosThetaNu2);
  double nu2Px_local = nu2P*TMath::Cos(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Py_local = nu2P*TMath::Sin(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Pz_local = nu2P*TMath::Cos(thetaNu2);
  double nu2Px = nu2Px_local*leg2eX_x_ + nu2Py_local*leg2eY_x_ + nu2Pz_local*leg2eZ_x_;
  double nu2Py = nu2Px_local*leg2eX_y_ + nu2Py_local*leg2eY_y_ + nu2Pz_local*leg2eZ_y_;
  double nu2Pz = nu2Px_local*leg2eX_z_ + nu2Py_local*leg2eY_z_ + nu2Pz_local*leg2eZ_z_;
  LorentzVector nu2P4(nu2Px, nu2Py, nu2Pz, nu2En);

  double tau2En = vis2En + nu2En;
  double tau2Px = vis2P4.px() + nu2Px;
  double tau2Py = vis2P4.py() + nu2Py;
  double tau2Pz = vis2P4.pz() + nu2Pz;
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);

  // evaluate transfer function for MET/hadronic recoil
  double sumNuPx = nu1Px + nu2Px;
  double sumNuPy = nu1Py + nu2Py;
  double residualX = measuredMETx_ - sumNuPx;
  double residualY = measuredMETy_ - sumNuPy;
  double pull2 = residualX*(invCovMETxx_*residualX + invCovMETxy_*residualY) + residualY*(invCovMETyx_*residualX + invCovMETyy_*residualY);
  double prob_TF = const_MET_*TMath::Exp(-0.5*pull2);

  // evaluate transfer functions for tau energy reconstruction
  if ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) {
    prob_TF *= extractProbFromLUT(x[idxLeg1VisPtShift_], leg1lutVisPtRes_);
  }
  if ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) {
    prob_TF *= extractProbFromLUT(x[idxLeg2VisPtShift_], leg2lutVisPtRes_);
  }

  // perform boost into MEM frame 
  double memFramePx = tau1Px + tau2Px;
  double memFramePy = tau1Py + tau2Py;
  double memFrameEn = tau1En + tau2En;
  Vector boost(-memFramePx/memFrameEn, -memFramePy/memFrameEn, 0.);
  //std::cout << "boost: Px = " << boost.x() << ", Py = " << boost.y() << ", Pz = " << boost.z() << std::endl;
  LorentzVector tau1P4_mem = ROOT::Math::VectorUtil::boost(tau1P4, boost); 
  LorentzVector tau2P4_mem = ROOT::Math::VectorUtil::boost(tau2P4, boost); 
  if ( verbosity_ >= 3 ) {
    std::cout << "lab:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4.pt() << ", eta = " << tau1P4.eta() << ", phi = " << tau1P4.phi() << ", mass = " << tau1P4.mass() << std::endl;
    std::cout << "      (En = " << tau1P4.energy() << ", Px = " << tau1P4.px() << ", Py = " << tau1P4.py() << ", Pz = " << tau1P4.pz() << ")" << std::endl;    
    std::cout << " tau2: Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << ", mass = " << tau2P4.mass() << std::endl;
    std::cout << "      (En = " << tau2P4.energy() << ", Px = " << tau2P4.px() << ", Py = " << tau2P4.py() << ", Pz = " << tau2P4.pz() << ")" << std::endl;  
    LorentzVector ditauP4 = tau1P4 + tau2P4;
    std::cout << " ditau: Pt = " << ditauP4.pt() << ", eta = " << ditauP4.eta() << ", phi = " << ditauP4.phi() << ", mass = " << ditauP4.mass() << std::endl;
    std::cout << "      (En = " << ditauP4.energy() << ", Px = " << ditauP4.px() << ", Py = " << ditauP4.py() << ", Pz = " << ditauP4.pz() << ")" << std::endl; 
    std::cout << "mem:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4_mem.pt() << ", eta = " << tau1P4_mem.eta() << ", phi = " << tau1P4_mem.phi() << ", mass = " << tau1P4_mem.mass() << std::endl;
    std::cout << "      (En = " << tau1P4_mem.energy() << ", Px = " << tau1P4_mem.px() << ", Py = " << tau1P4_mem.py() << ", Pz = " << tau1P4_mem.pz() << ")" << std::endl;  
    std::cout << " tau2: Pt = " << tau2P4_mem.pt() << ", eta = " << tau2P4_mem.eta() << ", phi = " << tau2P4_mem.phi() << ", mass = " << tau2P4_mem.mass() << std::endl;
    std::cout << "      (En = " << tau2P4_mem.energy() << ", Px = " << tau2P4_mem.px() << ", Py = " << tau2P4_mem.py() << ", Pz = " << tau2P4_mem.pz() << ")" << std::endl;  
    LorentzVector ditauP4_mem = tau1P4_mem + tau2P4_mem;
    std::cout << " ditau: Pt = " << ditauP4_mem.pt() << ", eta = " << ditauP4_mem.eta() << ", phi = " << ditauP4_mem.phi() << ", mass = " << ditauP4_mem.mass() << std::endl;
    std::cout << "      (En = " << ditauP4_mem.energy() << ", Px = " << ditauP4_mem.px() << ", Py = " << ditauP4_mem.py() << ", Pz = " << ditauP4_mem.pz() << ")" << std::endl; 
  }

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  LorentzVector ditauP4 = tau1P4_mem + tau2P4_mem;
  double ditauPz = ditauP4.pz();
  double ditauEn = ditauP4.E();
  double ditauMass = ditauP4.mass();
  //if ( !(ditauMass > 0.70*mTest_ && ditauMass < 1.30*mTest_) ) return 0.;
  //std::cout << "ditau: En = " << ditauEn << ", Pz = " << ditauPz << std::endl;
  double xa = invSqrtS_*(ditauEn + ditauPz);
  double xb = invSqrtS_*(ditauEn - ditauPz);   
  //std::cout << "xa = " << xa << ", xb = " << xb << std::endl;
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  //double Q = mTest_;
  double Q = ditauMass;
  assert(pdfIsInitialized_);
  double fa = LHAPDF::xfx(1, xa, Q, 0)/xa;
  double fb = LHAPDF::xfx(1, xb, Q, 0)/xb;
  double prob_PDF = (fa*fb);

  // evaluate flux factor
  double prob_flux = (1./(s_*xa*xb));  

  // evaluate LO matrix element, 
  // computed by Madgraph or taken from literature 
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
    if ( TMath::IsNaN(prob_ME_madgraph) ) {
      std::cerr << "Warning: Magraph returned NaN --> skipping event !!" << std::endl;
      std::cerr << " tau1: Pt = " << tau1P4.pt() << ", eta = " << tau1P4.eta() << ", phi = " << tau1P4.phi() << ", mass = " << tau1P4.mass() << std::endl;
      std::cerr << " tau2: Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << ", mass = " << tau2P4.mass() << std::endl;
      return 0.;
    }
  }
  me_lit_.setHiggsMass(mTest_);
  me_lit_.setHiggsWidth(GammaH_);
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
  const double constFactor = 2.*conversionFactor/eigth(2.*TMath::Pi());
  double prob_PS_and_tauDecay = constFactor/s_;
  LorentzVector vis1P4_mem = ROOT::Math::VectorUtil::boost(vis1P4, boost); 
  LorentzVector nu1P4_mem  = ROOT::Math::VectorUtil::boost(nu1P4, boost); 
  double x1_mem = vis1P4_mem.energy()/tau1P4_mem.energy();
  if ( !(x1_mem >= 1.e-5 && x1_mem <= 1.) ) return 0.;
  LorentzVector vis2P4_mem = ROOT::Math::VectorUtil::boost(vis2P4, boost); 
  LorentzVector nu2P4_mem  = ROOT::Math::VectorUtil::boost(nu2P4, boost); 
  double x2_mem = vis2P4_mem.energy()/tau2P4_mem.energy();
  if ( !(x2_mem >= 1.e-5 && x2_mem <= 1.) ) return 0.;
  double prob_tauDecay_leg1 = 0.;
  if ( leg1isLep_ ) {
    prob_tauDecay_leg1 = compPSfactor_tauToLepDecay(x1_mem, vis1P4_mem.E(), vis1P4_mem.P(), leg1Mass_, nu1P4_mem.E(), nu1P4_mem.P(), nu1Mass);
  } else {
    prob_tauDecay_leg1 = compPSfactor_tauToHadDecay(x1_mem, vis1P4_mem.E(), vis1P4_mem.P(), leg1Mass_, nu1P4_mem.E(), nu1P4_mem.P());
  }  
  prob_PS_and_tauDecay *= prob_tauDecay_leg1;
  double prob_tauDecay_leg2 = 0.;
  if ( leg2isLep_ ) {
    prob_tauDecay_leg2 = compPSfactor_tauToLepDecay(x2_mem, vis2P4_mem.E(), vis2P4_mem.P(), leg2Mass_, nu2P4_mem.E(), nu2P4_mem.P(), nu2Mass);
  } else {
    prob_tauDecay_leg2 = compPSfactor_tauToHadDecay(x2_mem, vis2P4_mem.E(), vis2P4_mem.P(), leg2Mass_, nu2P4_mem.E(), nu2P4_mem.P());
  }
  prob_PS_and_tauDecay *= prob_tauDecay_leg2;
  // CV: multiply matrix element by factor (Pi/(mTau GammaTau))^2 from Luca's write-up
  prob_PS_and_tauDecay *= square(TMath::Pi()/(tauLeptonMass*GammaTau));

  double jacobiFactor = (square(q2_ - mTest2_) + GammaH2_times_mTest2_)/GammaH_times_mTest_; // parametrization of q^2 by tk (Eq. 8 of arXiv:1010.2263 without factor 1/pi, as agreed with Luca and Andrew)
  double prob = prob_flux*prob_PDF*prob_ME*prob_PS_and_tauDecay*prob_TF*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
	      << " TF = " << prob_TF << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
  }

  return prob;
}
