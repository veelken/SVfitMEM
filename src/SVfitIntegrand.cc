#include "TauAnalysis/SVfitMEM/interface/SVfitIntegrand.h"

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"

#include <TMath.h>
#include <TH1.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  double xfx(int nset, double x, double Q, int fl);
}

const double GF = 1.166e-5; // in units of GeV^-2, taken from http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf

/// global function pointer, needed for VEGAS integration
const SVfitIntegrand* SVfitIntegrand::gSVfitIntegrand = 0;
bool SVfitIntegrand::pdfIsInitialized_ = false;

SVfitIntegrand::SVfitIntegrand(const std::string& madgraphFileName, double sqrtS, const std::string& pdfFileName, int verbosity) 
  : mTest2_(0.),
    beamAxis_(0., 0., 1.),
    invCovMET_(2,2),
    shiftVisPt_(false),
    leg1lutVisPtRes_(0),
    leg2lutVisPtRes_(0),
    idxLeg1_t_(-1),
    idxLeg1_phi_(-1),
    idxLeg1VisPtShift_(-1),
    idxLeg1_mNuNu_(-1),
    idxLeg2_X_(-1),
    idxLeg2_phi_(-1),
    idxLeg2VisPtShift_(-1),
    idxLeg2_mNuNu_(-1),
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
  madgraph_.initProc(madgraphFileName);

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

  // set center-of-mass energy
  s_ = square(sqrtS);
  invSqrtS_ = 1./sqrtS;

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
  //std::cout << "mTest = " << mTest_ << std::endl;
  mTest2_ = square(mTest); 
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
  leg1eX_z_ = eZ1.z();

  measuredTauLepton2_ = measuredTauLeptons[1];
  leg2isLep_ = measuredTauLepton2_.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton2_.type() == MeasuredTauLepton::kTauToMuDecay;
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
  leg2eX_z_ = eZ2.z();

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

  double compCosThetaNu(double visEn, double visP, double visMass2, double nuEn, double nuP, double nuMass2)
  {
    double cosThetaNu = (visEn*nuEn - 0.5*(tauLeptonMass2 - (visMass2 + nuMass2)))/(visP*nuP);
    return cosThetaNu;
  }

  double compPSfactor_tauToLepDecay(double x, double visEn, double visP, double visMass, double nuEn, double nuP, double nuMass)
  {
    if ( x >= 0. && x <= 1. && nuMass < TMath::Sqrt((1. - x)*tauLeptonMass2) ) { // physical solution
      const double GFfactor = square(GF)/square(TMath::Pi());
      double nuMass2 = square(nuMass);
      double visMass2 = square(visMass);
      double tauEn_rf = (tauLeptonMass2 + nuMass2 - visMass2)/(2.*nuMass);
      double visEn_rf = tauEn_rf - nuMass;
      double I = GFfactor*nuMass2*(2.*tauEn_rf*visEn_rf + (2./3.)*TMath::Sqrt((square(tauEn_rf) - tauLeptonMass2)*(square(visEn_rf) - visMass2)));
      double cosThetaNu = compCosThetaNu(visEn, visP, visMass2, nuEn, nuP, nuMass2);
      double PSfactor = (visEn + nuEn)*I/(8.*visP*square(x)*TMath::Sqrt(square(visP) + square(nuP) + 2.*visP*nuP*cosThetaNu + tauLeptonMass2));
      //std::cout << "<compPSfactor_tauToLepDecay>: PSfactor = " << PSfactor << " (I = " << I << ")" << std::endl;
      return PSfactor;
    } else {
      return 0.;
    }
  }

  double compPSfactor_tauToHadDecay(double x, double visEn, double visP, double visMass, double nuEn, double nuP)
  {
    if ( x >= (square(visMass)/tauLeptonMass2) && x <= 1. ) { // physical solution
      double cosThetaNu = compCosThetaNu(visEn, visP, square(visMass), nuEn, nuP, 0.);
      double PSfactor = (visEn + nuEn)/(8.*visP*square(x)*TMath::Sqrt(square(visP) + square(nuP) + 2.*visP*nuP*cosThetaNu + tauLeptonMass2));
      //std::cout << "<compPSfactor_tauToHadDecay>: PSfactor = " << PSfactor << std::endl;
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

  assert(idxLeg2_X_ != -1);
  double x2 = x[idxLeg2_X_];
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;

  assert(idxLeg1_t_ != -1);  
  double x1 = (mVis2_*visPtShift1*visPtShift2)/((mTest2_ + 1.e-2*mTest2_*TMath::Tan(x[idxLeg1_t_]))*x2);
  //std::cout << "x2 = " << x2 << ", t1 = " << x[idxLeg1_t_] << " --> x1 = " << x1 << std::endl;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  // compute neutrino and tau lepton four-vector for first tau
  double vis1En   = visPtShift1*measuredTauLepton1_.energy();
  double vis1Mass = measuredTauLepton1_.mass();
  double vis1P    = TMath::Sqrt(TMath::Max(0., vis1En*vis1En - square(vis1Mass)));
  double vis1Px   = vis1P*measuredTauLepton1_.cosPhi_sinTheta();
  double vis1Py   = vis1P*measuredTauLepton1_.sinPhi_sinTheta();
  double vis1Pz   = vis1P*measuredTauLepton1_.cosTheta();

  double nu1En    = vis1En*(1. - x1)/x1;
  //std::cout << "vis1En = " << vis1En << ", x1 = " << x1 << " --> nu1En = " << nu1En << std::endl;
  double nu1Mass  = ( idxLeg1_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg1_mNuNu_]) : 0.;
  //std::cout << "nu1Mass = " << nu1Mass << std::endl;
  double nu1P     = TMath::Sqrt(TMath::Max(0., nu1En*nu1En - square(nu1Mass)));
  assert(idxLeg1_phi_ != -1);  
  double phiNu1 = x[idxLeg1_phi_];
  double cosThetaNu1 = compCosThetaNu(vis1En, vis1P, square(vis1Mass), nu1En, nu1P, square(nu1Mass));
  //std::cout << "cosThetaNu1 = " << cosThetaNu1 << std::endl;
  if ( !(cosThetaNu1 >= -1. && cosThetaNu1 <= +1.) ) return 0.;
  double thetaNu1 = TMath::ACos(cosThetaNu1);
  double nu1Px_local = nu1P*TMath::Cos(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Py_local = nu1P*TMath::Sin(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Pz_local = nu1P*TMath::Cos(thetaNu1);
  double nu1Px = nu1Px_local*leg1eX_x_ + nu1Py_local*leg1eY_x_ + nu1Pz_local*leg1eZ_x_;
  double nu1Py = nu1Px_local*leg1eX_y_ + nu1Py_local*leg1eY_y_ + nu1Pz_local*leg1eZ_y_;
  double nu1Pz = nu1Px_local*leg1eX_z_ + nu1Py_local*leg1eY_z_ + nu1Pz_local*leg1eZ_z_;

  double tau1En = vis1En + nu1En;
  double tau1Px = vis1Px + nu1Px;
  double tau1Py = vis1Py + nu1Py;
  double tau1Pz = vis1Pz + nu1Pz;
  
  // compute neutrino and tau lepton four-vector for second tau
  double vis2En   = visPtShift2*measuredTauLepton2_.energy();
  double vis2Mass = measuredTauLepton2_.mass();
  double vis2P    = TMath::Sqrt(TMath::Max(0., vis2En*vis2En - square(vis2Mass)));
  double vis2Px   = vis2P*measuredTauLepton2_.cosPhi_sinTheta();
  double vis2Py   = vis2P*measuredTauLepton2_.sinPhi_sinTheta();
  double vis2Pz   = vis2P*measuredTauLepton2_.cosTheta();

  double nu2En    = vis2En*(1. - x2)/x2;
  //std::cout << "vis2En = " << vis2En << ", x2 = " << x2 << " --> nu2En = " << nu2En << std::endl;
  double nu2Mass  = ( idxLeg2_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg2_mNuNu_]) : 0.;
  //std::cout << "nu2Mass = " << nu2Mass << std::endl;
  double nu2P     = TMath::Sqrt(TMath::Max(0., nu2En*nu2En - square(nu2Mass)));
  assert(idxLeg2_phi_ != -2);  
  double phiNu2 = x[idxLeg2_phi_];
  double cosThetaNu2 = compCosThetaNu(vis2En, vis2P, square(vis2Mass), nu2En, nu2P, square(nu2Mass));
  //std::cout << "cosThetaNu2 = " << cosThetaNu2 << std::endl;
  if ( !(cosThetaNu2 >= -1. && cosThetaNu2 <= +1.) ) return 0.;
  double thetaNu2 = TMath::ACos(cosThetaNu2);
  double nu2Px_local = nu2P*TMath::Cos(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Py_local = nu2P*TMath::Sin(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Pz_local = nu2P*TMath::Cos(thetaNu2);
  double nu2Px = nu2Px_local*leg2eX_x_ + nu2Py_local*leg2eY_x_ + nu2Pz_local*leg2eZ_x_;
  double nu2Py = nu2Px_local*leg2eX_y_ + nu2Py_local*leg2eY_y_ + nu2Pz_local*leg2eZ_y_;
  double nu2Pz = nu2Px_local*leg2eX_z_ + nu2Py_local*leg2eY_z_ + nu2Pz_local*leg2eZ_z_;

  double tau2En = vis2En + nu2En;
  double tau2Px = vis2Px + nu2Px;
  double tau2Py = vis2Py + nu2Py;
  double tau2Pz = vis2Pz + nu2Pz;

  // evaluate transfer function for MET/hadronic recoil
  double sumNuPx = nu1Px + nu2Px;
  double sumNuPy = nu1Py + nu2Py;
  double residualX = measuredMETx_ - sumNuPx;
  double residualY = measuredMETy_ - sumNuPy;
  double pull2 = residualX*(invCovMETxx_*residualX + invCovMETxy_*residualY) + residualY*(invCovMETyx_*residualX + invCovMETyy_*residualY);
  double prob_TF = const_MET_*TMath::Exp(-0.5*pull2);
  
  // evaluate transfer functions for tau energy reconstruction
  if ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) prob_TF *= extractProbFromLUT(x[idxLeg1VisPtShift_], leg1lutVisPtRes_);
  if ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) prob_TF *= extractProbFromLUT(x[idxLeg2VisPtShift_], leg2lutVisPtRes_);

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  double ditauPz = tau1Pz + tau2Pz;
  double ditauEn = tau1En + tau2En;
  double xa = invSqrtS_*(ditauEn + ditauPz);
  double xb = invSqrtS_*(ditauEn - ditauPz);    
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  double Q = mTest_;
  assert(pdfIsInitialized_);
  double fa = LHAPDF::xfx(1, xa, Q, 0)/xa;
  double fb = LHAPDF::xfx(1, xb, Q, 0)/xb;
  double prob_PDF = (fa*fb);

  // evaluate flux factor
  double prob_flux = (1./(s_*xa*xb));  
  
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
    std::cout << " tau2: Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << ", mass = " << tau2P4.mass() << std::endl;
    std::cout << "rf:" << std::endl;
    std::cout << " tau1: Pt = " << tau1P4_rf.pt() << ", eta = " << tau1P4_rf.eta() << ", phi = " << tau1P4_rf.phi() << ", mass = " << tau1P4_rf.mass() << std::endl;
    std::cout << " tau2: Pt = " << tau2P4_rf.pt() << ", eta = " << tau2P4_rf.eta() << ", phi = " << tau2P4_rf.phi() << ", mass = " << tau2P4_rf.mass() << std::endl;
  }
  madgraphGluon1P4_[0] =  0.5*mTest_;
  madgraphGluon1P4_[3] = +0.5*mTest_;
  madgraphGluon2P4_[0] =  0.5*mTest_;
  madgraphGluon2P4_[3] = -0.5*mTest_;
  madgraphTau1P4_[0] = tau1P4_rf.energy();
  madgraphTau1P4_[1] = tau1P4_rf.px();
  madgraphTau1P4_[2] = tau1P4_rf.py();
  madgraphTau1P4_[3] = tau1P4_rf.pz();
  madgraphTau2P4_[0] = tau2P4_rf.energy();
  madgraphTau2P4_[1] = tau2P4_rf.px();
  madgraphTau2P4_[2] = tau2P4_rf.py();
  madgraphTau2P4_[3] = tau2P4_rf.pz();
  //madgraph_.setHiggsMass(mTest_);
  //madgraph_.setHiggsWidth(1.e-2*mTest_);
  madgraph_.setMomenta(madgraphMomenta_);
  madgraph_.sigmaKin();
  double prob_ME = madgraph_.getMatrixElements()[0];
  
  const double twoPiFactor = 1./sixth(2.*TMath::Pi());
  double prob_PS_and_tauDecay = twoPiFactor*twoPiFactor;
  if ( leg1isLep_ ) {
    prob_PS_and_tauDecay *= compPSfactor_tauToLepDecay(x1, vis1En, vis1P, vis1Mass, nu1En, nu1P, nu1Mass);
  } else {
    prob_PS_and_tauDecay *= compPSfactor_tauToHadDecay(x1, vis1En, vis1P, vis1Mass, nu1En, nu1P);
  }  
  if ( leg2isLep_ ) {
    prob_PS_and_tauDecay *= compPSfactor_tauToLepDecay(x2, vis2En, vis2P, vis2Mass, nu2En, nu2P, nu2Mass);
  } else {
    prob_PS_and_tauDecay *= compPSfactor_tauToHadDecay(x2, vis2En, vis2P, vis2Mass, nu2En, nu2P);
  }

  double Gamma_times_m = 1.e-2*mTest_*mTest_;
  double q2 = mTest2_ + 1.e-2*mTest2_*TMath::Tan(x[idxLeg1_t_]);
  double jacobiFactor = ((2.*mVis2_*visPtShift1*visPtShift2)/(square(q2)*x2))*(square(comSystem.M2() - mTest2_) + square(Gamma_times_m))/Gamma_times_m; 
  if ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) jacobiFactor *= measuredTauLepton1_.pt();
  if ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) jacobiFactor *= measuredTauLepton2_.pt();

  double prob = prob_flux*prob_PDF*prob_ME*prob_PS_and_tauDecay*prob_TF*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
	      << " TF = " << prob_TF << " --> returning " << prob << std::endl;
  }
  return prob;
}
