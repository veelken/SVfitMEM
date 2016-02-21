#include "TauAnalysis/SVfitMEM/interface/SVfitIntegrand_nlo.h"

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVEGAS.h"

#include <TMath.h>
#include <TH1.h>
#include <Math/VectorUtil.h>

#include <math.h>

using namespace svFitMEM;

/// global function pointer, needed for VEGAS integration
const SVfitIntegrand_nlo* SVfitIntegrand_nlo::gSVfitIntegrand = 0;

SVfitIntegrand_nlo::SVfitIntegrand_nlo(double sqrtS, const std::string& pdfName, int verbosity) 
  : minLeg1Pt_(0.),
    maxLeg1AbsEta_(1.e+3),
    minLeg2Pt_(0.),
    maxLeg2AbsEta_(1.e+3),
    useAcceptanceCuts_(false),
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
    hadTauTF1_(0),
    hadTauTF2_(0),
    useHadTauTF_(false),
    rhoHadTau_(0.),
    idxLeg1_X_(-1),
    idxLeg1_phi_(-1),
    idxLeg1VisPtShift_(-1),
    idxLeg1_mNuNu_(-1),
    idxLeg2_t_(-1),
    idxLeg2_phi_(-1),
    idxLeg2VisPtShift_(-1),
    idxLeg2_mNuNu_(-1),
    numDimensions_(0),
    pdf_(0),
    pdfIsInitialized_(false),
    me_lit_(pdfName, false, true),
    errorCode_(0),
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<SVfitIntegrand_nlo::SVfitIntegrand_nlo>:" << std::endl;
  }

  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    pdf_ = LHAPDF::mkPDF(pdfName.data(), 0);
    pdfIsInitialized_ = true;
  }

  me_lit_.setS(s_);
  // set Higgs -> tautau decay branching fraction 
  // to common value for all test-mass hypotheses
  me_lit_.setBR(1.e-1); 

  momenta_.push_back(&parton1P4_incoming_);
  momenta_.push_back(&parton2P4_incoming_);
  momenta_.push_back(&tau1P4_outgoing_);
  momenta_.push_back(&tau2P4_outgoing_);
  momenta_.push_back(&jetP4_outgoing_);

  // set global function pointer to this
  gSVfitIntegrand = this;
}

SVfitIntegrand_nlo::~SVfitIntegrand_nlo()
{
  std::cout << "<SVfitIntegrand_nlo::~SVfitIntegrand_nlo>:" << std::endl;
  
  delete hadTauTF1_;
  delete hadTauTF2_;

  delete pdf_;
}

void 
SVfitIntegrand_nlo::setMtest(double mTest) 
{ 
  // reset 'TestMass' error code
  errorCode_ &= (errorCode_ ^ TestMass);
  
  if ( mTest < TMath::Sqrt(mVis2_measured_) ) {
    std::cerr << "Error: Cannot have mTest < mVis = " << TMath::Sqrt(mVis2_measured_) << " !!" << std::endl;
    errorCode_ |= TestMass;
    return;
  }
  mTest_ = mTest; 
  mTest2_ = square(mTest); 
  GammaH_ = 1.e-2*mTest_;
  GammaH_times_mTest_ = GammaH_*mTest_;
  GammaH2_times_mTest2_ = square(GammaH_times_mTest_);
}

namespace
{
  double norm(const Vector& v)
  {
    return TMath::Sqrt(v.mag2());
  }
}

void 
SVfitIntegrand_nlo::setInputs(const std::vector<MeasuredTauLepton>& measuredTauLeptons, double measuredHadRecoilPx, double measuredHadRecoilPy, double measuredHadRecoilPz, const TMatrixD& covHadRecoil) 
{
  if ( verbosity_ ) {
    std::cout << "<SVfitIntegrand_nlo::setInputs>:" << std::endl;
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

  mVis_measured_ = (measuredTauLepton1_.p4() + measuredTauLepton2_.p4()).mass();
  if ( verbosity_ >= 2 ) {
    std::cout << "mVis = " << mVis_measured_ << std::endl;
  }
  mVis2_measured_ = square(mVis_measured_);

  measuredHadRecoilPx_ = measuredHadRecoilPx;
  measuredHadRecoilPy_ = measuredHadRecoilPy;
  measuredHadRecoilPz_ = measuredHadRecoilPz;

  measuredMETx_ = -(measuredTauLepton1_.px() + measuredTauLepton2_.px() + measuredHadRecoilPx);
  measuredMETy_ = -(measuredTauLepton1_.py() + measuredTauLepton2_.py() + measuredHadRecoilPy);

  // determine transfer matrix for hadronic recoil/MET
  invCovMET_ = covHadRecoil;
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

  if ( useHadTauTF_ ) {
    if ( measuredTauLepton1_.type() == MeasuredTauLepton::kTauToHadDecay ) { 
      assert(hadTauTF1_);
      hadTauTF1_->setDecayMode(measuredTauLepton1_.decayMode());
    }
    if ( measuredTauLepton2_.type() == MeasuredTauLepton::kTauToHadDecay ) { 
      assert(hadTauTF2_);
      hadTauTF2_->setDecayMode(measuredTauLepton2_.decayMode());
    }
  }
}

double
SVfitIntegrand_nlo::Eval(const double* x) const 
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<SVfitIntegrand_nlo::Eval(const double*)>:" << std::endl;
    std::cout << " x = { ";
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << x[iDimension];
      if ( iDimension < (numDimensions_ - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  // in case of initialization errors don't start to do anything
  if ( errorCode_ != 0 ) { 
    return 0.;
  } 

  double visPtShift1 = ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) ? (1./x[idxLeg1VisPtShift_]) : 1.;
  double visPtShift2 = ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) ? (1./x[idxLeg2VisPtShift_]) : 1.;
  if ( visPtShift1 < 1.e-2 || visPtShift2 < 1.e-2 ) return 0.;

  // compute four-vector of visible decay products for first tau
  double vis1Px = visPtShift1*measuredTauLepton1_.px();
  double vis1Py = visPtShift1*measuredTauLepton1_.py();
  double vis1Pz = visPtShift1*measuredTauLepton1_.pz();
  double vis1En = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  //std::cout << "vis1: En = " << vis1En << ", Pt = " << TMath::Sqrt(vis1Px*vis1Px + vis1Py*vis1Py) << std::endl;
  LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);
  double vis1P = vis1P4.P();

  // compute four-vector of visible decay products for second tau
  double vis2Px = visPtShift2*measuredTauLepton2_.px();
  double vis2Py = visPtShift2*measuredTauLepton2_.py();
  double vis2Pz = visPtShift2*measuredTauLepton2_.pz();
  double vis2En = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  //std::cout << "vis2: En = " << vis2En << ", Pt = " << TMath::Sqrt(vis2Px*vis2Px + vis2Py*vis2Py) << std::endl;
  LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);
  double vis2P = vis2P4.P();

  // compute visible energy fractions for both taus
  assert(idxLeg1_X_ != -1);
  double x1_dash = x[idxLeg1_X_];
  double x1 = x1_dash/visPtShift1;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;
  
  assert(idxLeg2_t_ != -1);  
  double tk = x[idxLeg2_t_];
  q2_ = mTest2_ + GammaH_times_mTest_*TMath::Tan(tk);
  if ( q2_ <= 0. ) return 0.;

  double x2_dash = mVis2_measured_/(q2_*x1_dash);
  double x2 = x2_dash/visPtShift2;
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
  //std::cout << "nu1: En = " << nu1En << ", Pt = " << TMath::Sqrt(nu1Px*nu1Px + nu1Py*nu1Py) << std::endl;
  LorentzVector nu1P4(nu1Px, nu1Py, nu1Pz, nu1En);

  double tau1En = vis1En + nu1En;
  double tau1Px = vis1P4.px() + nu1Px;
  double tau1Py = vis1P4.py() + nu1Py;
  double tau1Pz = vis1P4.pz() + nu1Pz;
  //std::cout << "tau1: En = " << tau1En << ", Pt = " << TMath::Sqrt(tau1Px*tau1Px + tau1Py*tau1Py) << std::endl;
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
  //std::cout << "nu2: En = " << nu2En << ", Pt = " << TMath::Sqrt(nu2Px*nu2Px + nu2Py*nu2Py) << std::endl;
  LorentzVector nu2P4(nu2Px, nu2Py, nu2Pz, nu2En);

  double tau2En = vis2En + nu2En;
  double tau2Px = vis2P4.px() + nu2Px;
  double tau2Py = vis2P4.py() + nu2Py;
  double tau2Pz = vis2P4.pz() + nu2Pz;
  //std::cout << "tau2: En = " << tau2En << ", Pt = " << TMath::Sqrt(tau2Px*tau2Px + tau2Py*tau2Py) << std::endl;
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);

  // evaluate transfer function for MET/hadronic recoil
  double sumNuPx = nu1Px + nu2Px;
  double sumNuPy = nu1Py + nu2Py;
  double residualX = measuredMETx_ - sumNuPx;
  double residualY = measuredMETy_ - sumNuPy;
  if ( rhoHadTau_ != 0. ) {
    residualX += (rhoHadTau_*((visPtShift1 - 1.)*measuredTauLepton1_.px() + (visPtShift2 - 1.)*measuredTauLepton2_.px()));
    residualY += (rhoHadTau_*((visPtShift1 - 1.)*measuredTauLepton1_.py() + (visPtShift2 - 1.)*measuredTauLepton2_.py()));
  }
  double pull2 = residualX*(invCovMETxx_*residualX + invCovMETxy_*residualY) + residualY*(invCovMETyx_*residualX + invCovMETyy_*residualY);
  double prob_TF_met = const_MET_*TMath::Exp(-0.5*pull2);

  double prob_TF = prob_TF_met;

  // evaluate transfer functions for tau energy reconstruction
  if ( useHadTauTF_ && idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) {
    double genVis1Pt = vis1P4.pt();
    double genVis1Eta = vis1P4.eta();
    double prob_TF_leg1 = (*hadTauTF1_)(measuredTauLepton1_.pt(), genVis1Pt, genVis1Eta);
    if ( verbosity_ >= 2 ) {
      std::cout << "TF(leg1): recPt = " << measuredTauLepton1_.pt() << ", genPt = " << genVis1Pt << ", genEta = " << genVis1Eta << " --> prob = " << prob_TF_leg1 << std::endl;    
    }
    prob_TF *= prob_TF_leg1;
    if ( useAcceptanceCuts_ ) {
      double epsilon_leg1 = 1.;
      if ( TMath::Abs(genVis1Eta) < maxLeg1AbsEta_ ) {
	epsilon_leg1 = hadTauTF1_->integral(minLeg1Pt_, 2.*genVis1Pt, genVis1Pt, genVis1Eta);
      } else {
	epsilon_leg1 = 0.;
      }
      prob_TF *= epsilon_leg1;
    }
  }
  if ( useHadTauTF_ && idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) {
    double genVis2Pt = vis2P4.pt();
    double genVis2Eta = vis2P4.eta();
    double prob_TF_leg2 = (*hadTauTF2_)(measuredTauLepton2_.pt(), genVis2Pt, genVis2Eta);
    if ( verbosity_ >= 2 ) {
      std::cout << "TF(leg2): recPt = " << measuredTauLepton2_.pt() << ", genPt = " << genVis2Pt << ", genEta = " << genVis2Eta << " --> prob = " << prob_TF_leg2 << std::endl;
    }
    prob_TF *= prob_TF_leg2;
    if ( useAcceptanceCuts_ ) {
      double epsilon_leg2 = 1.;
      if ( TMath::Abs(genVis2Eta) < maxLeg2AbsEta_ ) {
	epsilon_leg2 = hadTauTF2_->integral(minLeg2Pt_, 2.*genVis2Pt, genVis2Pt, genVis2Eta);
      } else {
	epsilon_leg2 = 0.;
      }
      prob_TF *= epsilon_leg2;
    }
  }

  // compute jet four-vector 
  double jetPx = -(tau1Px + tau2Px);
  double jetPy = -(tau1Py + tau2Py);
  double jetPz = measuredHadRecoilPz_;
  double jetEn = TMath::Sqrt(jetPx*jetPx + jetPy*jetPy + jetPz*jetPz);
  LorentzVector jetP4(jetPx, jetPy, jetPz, jetEn);

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  LorentzVector ditauP4 = tau1P4 + tau2P4;
  double ditauPz = ditauP4.pz();
  double ditauEn = ditauP4.E();
  double ditauMass = ditauP4.mass();
  double xa = invSqrtS_*(ditauEn + jetEn + (ditauPz + jetPz));
  double xb = invSqrtS_*(ditauEn + jetEn - (ditauPz + jetPz));
  //std::cout << "xa = " << xa << ", xb = " << xb << std::endl;
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  double Q = ditauMass;
  assert(pdfIsInitialized_);
  double fa_g = pdf_->xfxQ(21, xa, Q)/xa; // gluon 
  double fa_u = pdf_->xfxQ( 2, xa, Q)/xa; // up-quark
  double fa_d = pdf_->xfxQ( 1, xa, Q)/xa; // down-quark
  double fb_g = pdf_->xfxQ(21, xb, Q)/xb;
  double fb_u = pdf_->xfxQ( 2, xb, Q)/xb;
  double fb_d = pdf_->xfxQ( 1, xb, Q)/xb;

  // evaluate flux factor
  double prob_flux = (1./(s_*xa*xb));  

  // evaluate LO matrix element for gg -> Higgs+jet, 
  // taken from literature (arXiv:hep-ph/0201114)
  parton1P4_incoming_.SetPxPyPzE(0., 0., +0.5*xa*sqrtS_, 0.5*xa*sqrtS_);
  parton2P4_incoming_.SetPxPyPzE(0., 0., -0.5*xb*sqrtS_, 0.5*xb*sqrtS_);
  tau1P4_outgoing_ = tau1P4;
  tau2P4_outgoing_ = tau2P4;
  jetP4_outgoing_ = jetP4;
  me_lit_.setHiggsMass(mTest_);
  me_lit_.setHiggsWidth(GammaH_);
  me_lit_.setMomenta(momenta_);
  double prob_ME_lit_gg = me_lit_.getMatrixElement_gg();
  double prob_ME_lit_qg = me_lit_.getMatrixElement_qg();
  if ( verbosity_ >= 2 ) {
    std::cout << "prob_ME_lit: gg = " << prob_ME_lit_gg << ", qg = " << prob_ME_lit_qg << std::endl;
  }
  double prob_PDF_times_ME = fa_g*fb_g*prob_ME_lit_gg + ((fa_u + fa_d)*fb_g + fa_g*(fb_u + fb_d))*prob_ME_lit_qg;
  assert(prob_PDF_times_ME >= 0.);
  
  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  const double constFactor = 2.*conversionFactor/eigth(2.*TMath::Pi());
  double prob_PS_and_tauDecay = constFactor/s_;
  double prob_tauDecay_leg1 = 0.;
  if ( leg1isLep_ ) {
    prob_tauDecay_leg1 = compPSfactor_tauToLepDecay(x1, vis1P4.E(), vis1P4.P(), leg1Mass_, nu1P4.E(), nu1P4.P(), nu1Mass);
  } else {
    prob_tauDecay_leg1 = compPSfactor_tauToHadDecay(x1, vis1P4.E(), vis1P4.P(), leg1Mass_, nu1P4.E(), nu1P4.P());
  }  
  prob_PS_and_tauDecay *= prob_tauDecay_leg1;
  double prob_tauDecay_leg2 = 0.;
  if ( leg2isLep_ ) {
    prob_tauDecay_leg2 = compPSfactor_tauToLepDecay(x2, vis2P4.E(), vis2P4.P(), leg2Mass_, nu2P4.E(), nu2P4.P(), nu2Mass);
  } else {
    prob_tauDecay_leg2 = compPSfactor_tauToHadDecay(x2, vis2P4.E(), vis2P4.P(), leg2Mass_, nu2P4.E(), nu2P4.P());
  }
  prob_PS_and_tauDecay *= prob_tauDecay_leg2;
  // CV: multiply matrix element by factor (Pi/(mTau GammaTau))^2 from Luca's write-up
  prob_PS_and_tauDecay *= square(TMath::Pi()/(tauLeptonMass*GammaTau));

  if ( q2_ <= 0. || x1_dash <= 0. ) return 0.;
  double jacobiFactor = 1./(visPtShift1*visPtShift2);                                  // product of derrivatives dx1/dx1' and dx2/dx2' for parametrization of x1, x2 by x1', x2'
  jacobiFactor *= mVis2_measured_/(square(q2_)*x1_dash);                               // derrivative dx2'/dq^2 for parametrization of x2' by q^2
  jacobiFactor *= (square(q2_ - mTest2_) + GammaH2_times_mTest2_)/GammaH_times_mTest_; // parametrization of q^2 by tk (Eq. 8 of arXiv:1010.2263 without factor 1/pi, as agreed with Luca and Andrew)
  double prob = prob_flux*prob_PDF_times_ME*prob_PS_and_tauDecay*prob_TF*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF*ME = " << prob_PDF_times_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
	      << " TF = " << prob_TF << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
  }
  if ( TMath::IsNaN(prob) ) {
    std::cerr << "Warning: prob = " << prob << " (flux = " << prob_flux << ", PDF*ME = " << prob_PDF_times_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
	      << " TF = " << prob_TF << ", Jacobi = " << jacobiFactor << ") --> setting prob = 0 !!" << std::endl;
    prob = 0.;
  }

  return prob;
}
