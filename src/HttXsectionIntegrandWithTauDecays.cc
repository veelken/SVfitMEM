#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandWithTauDecays.h"

#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

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
const HttXsectionIntegrandWithTauDecays* HttXsectionIntegrandWithTauDecays::gHttXsectionIntegrandWithTauDecays = 0;
bool HttXsectionIntegrandWithTauDecays::pdfIsInitialized_ = false;

HttXsectionIntegrandWithTauDecays::HttXsectionIntegrandWithTauDecays(const std::string& madgraphFileName, double sqrtS, const std::string& pdfFileName, int verbosity) 
  : mTest2_(0.),
    GammaH_(0.),
    s_(square(sqrtS)),
    invSqrtS_(1./sqrtS),
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
    me_madgraph_(false),
    me_lit_(false),
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandWithTauDecays::HttXsectionIntegrandWithTauDecays>:" << std::endl;
  }

  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    LHAPDF::initPDFSet(1, pdfFileName);
    pdfIsInitialized_ = true;
  }

  // initialize Madgraph
  me_madgraph_.initProc(madgraphFileName);

  me_lit_.setS(s_);
  me_lit_.setHiggsMass(mTest_);
  me_lit_.setHiggsWidth(GammaH_);
  me_lit_.setBR(1.);

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
  gHttXsectionIntegrandWithTauDecays = this;
}

HttXsectionIntegrandWithTauDecays::~HttXsectionIntegrandWithTauDecays()
{
  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphTau1P4_;
  delete [] madgraphTau2P4_;
}

void 
HttXsectionIntegrandWithTauDecays::setMtest(double mTest) 
{ 
  mTest_ = mTest; 
  mTest2_ = square(mTest); 
  
  //GammaH_ = 1.e-2*mTest;
  GammaH_ = 1.;
}

void 
HttXsectionIntegrandWithTauDecays::shiftVisPt(bool value, const TH1* leg1lutVisPtRes, const TH1* leg2lutVisPtRes)
{
  shiftVisPt_ = value;
  if ( shiftVisPt_ ) {
    leg1lutVisPtRes_ = leg1lutVisPtRes;
    leg2lutVisPtRes_ = leg2lutVisPtRes;
  }
}

void 
HttXsectionIntegrandWithTauDecays::setInputs(int tau1Type, double tau1Mass, int tau2Type, double tau2Mass, const TMatrixD& covMET)
{
  if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandWithTauDecays::setInputs>:" << std::endl;
  }

  leg1isLep_ = tau1Type == MeasuredTauLepton::kTauToElecDecay || tau1Type == MeasuredTauLepton::kTauToMuDecay;
  leg1Mass_ = tau1Mass;
  leg1Mass2_ = square(leg1Mass_);
  
  leg2isLep_ = tau2Type == MeasuredTauLepton::kTauToElecDecay || tau2Type == MeasuredTauLepton::kTauToMuDecay;
  leg2Mass_ = tau2Mass;
  leg2Mass2_ = square(leg2Mass_);

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
  } else {
    std::cerr << "Warning: MET covariance matrix cannot be inverted !!" << std::endl;
    assert(0.);
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
      return PSfactor;
    } else {
      return 0.;
    }
  }
}

double
HttXsectionIntegrandWithTauDecays::Eval(const double* x) const 
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<HttXsectionIntegrandWithTauDecays::Eval(const double*)>:" << std::endl;
  }

  double visPtShift1 = ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) ? 1./x[idxLeg1VisPtShift_] : 1.;
  double visPtShift2 = ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) ? 1./x[idxLeg2VisPtShift_] : 1.;

  // compute four-vector of visible decay products for first tau
  double vis1Theta = x[1];
  double vis1Phi   = x[2];
  Vector eZ1 = Vector(TMath::Cos(vis1Phi)*TMath::Sin(vis1Theta), TMath::Sin(vis1Phi)*TMath::Sin(vis1Theta), TMath::Cos(vis1Theta));
  Vector eY1 = normalize(compCrossProduct(eZ1, beamAxis_));
  Vector eX1 = normalize(compCrossProduct(eY1, eZ1));
  leg1eX_x_ = eX1.x();
  leg1eX_y_ = eX1.y();
  leg1eX_z_ = eX1.z();
  leg1eY_x_ = eY1.x();
  leg1eY_y_ = eY1.y();
  leg1eY_z_ = eY1.z();
  leg1eZ_x_ = eZ1.x();
  leg1eZ_y_ = eZ1.y();
  leg1eZ_z_ = eZ1.z();
  double vis1En    = visPtShift1*x[0];
  double vis1P     = TMath::Sqrt(TMath::Max(0., vis1En*vis1En - leg1Mass2_));
  double vis1Px    = vis1P*eZ1.x();
  double vis1Py    = vis1P*eZ1.y();
  double vis1Pz    = vis1P*eZ1.z();
  LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);

  // compute four-vector of visible decay products for second tau
  double vis2Theta = x[4];
  double vis2Phi   = x[5];
  Vector eZ2 = Vector(TMath::Cos(vis2Phi)*TMath::Sin(vis2Theta), TMath::Sin(vis2Phi)*TMath::Sin(vis2Theta), TMath::Cos(vis2Theta));
  Vector eY2 = normalize(compCrossProduct(eZ2, beamAxis_));
  Vector eX2 = normalize(compCrossProduct(eY2, eZ2));
  leg2eX_x_ = eX2.x();
  leg2eX_y_ = eX2.y();
  leg2eX_z_ = eX2.z();
  leg2eY_x_ = eY2.x();
  leg2eY_y_ = eY2.y();
  leg2eY_z_ = eY2.z();
  leg2eZ_x_ = eZ2.x();
  leg2eZ_y_ = eZ2.y();
  leg2eZ_z_ = eZ2.z();
  double vis2En    = visPtShift2*(x[3]/(2.*vis1En*TMath::Max(1.e-2, 1. - compScalarProduct(eZ1, eZ2))));
  double vis2P     = TMath::Sqrt(TMath::Max(0., vis2En*vis2En - leg2Mass2_));
  double vis2Px    = vis2P*eZ2.x();
  double vis2Py    = vis2P*eZ2.y();
  double vis2Pz    = vis2P*eZ2.z();
  LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);

  double mVis = (vis1P4 + vis2P4).mass();
  double mVis2 = square(mVis);

  assert(idxLeg2_X_ != -1);
  double x2 = x[idxLeg2_X_];
  //std::cout << "x2 = " << x2 << std::endl;
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;
  
  //GammaH_ = me_madgraph_.getHiggsWidth();
  double GammaH_times_mH = GammaH_*mTest_;
  double tk = x[idxLeg1_t_];
  double q2 = mTest2_ + GammaH_times_mH*TMath::Tan(tk);
  //std::cout << "tk = " << tk << ": q2 = " << q2 << " (mVis2 = " << mVis2 << ", visPtShift1 = " << visPtShift1 << ", visPtShift2 = " << visPtShift2 << ")" << std::endl;
  if ( q2 <= 0. ) return 0.;
  
  assert(idxLeg1_t_ != -1);  
  double x1 = (mVis2*visPtShift1*visPtShift2)/(q2*x2);
  //std::cout << "x1 = " << x1 << std::endl;
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  // compute neutrino and tau lepton four-vector for first tau
  double nu1En    = vis1En*(1. - x1)/x1;
  double nu1Mass  = ( idxLeg1_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg1_mNuNu_]) : 0.;
  double nu1P     = TMath::Sqrt(TMath::Max(0., nu1En*nu1En - square(nu1Mass)));
  //std::cout << "nu1P = " << nu1P << std::endl;
  assert(idxLeg1_phi_ != -1);  
  double phiNu1 = x[idxLeg1_phi_];
  double cosThetaNu1 = compCosThetaNu(vis1En, vis1P, leg1Mass2_, nu1En, nu1P, square(nu1Mass));
  //std::cout << "cosThetaNu1 = " << cosThetaNu1 << std::endl;
  if ( !(cosThetaNu1 >= -1. && cosThetaNu1 <= +1.) ) return 0.;
  double thetaNu1 = TMath::ACos(cosThetaNu1);
  double nu1Px_local = nu1P*TMath::Cos(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Py_local = nu1P*TMath::Sin(phiNu1)*TMath::Sin(thetaNu1);
  double nu1Pz_local = nu1P*TMath::Cos(thetaNu1);
  //std::cout << "nu1_local: Px = " << nu1Px_local << ", Py = " << nu1Py_local << ", Pz = " << nu1Pz_local << std::endl;
  double nu1Px = nu1Px_local*leg1eX_x_ + nu1Py_local*leg1eY_x_ + nu1Pz_local*leg1eZ_x_;
  double nu1Py = nu1Px_local*leg1eX_y_ + nu1Py_local*leg1eY_y_ + nu1Pz_local*leg1eZ_y_;
  double nu1Pz = nu1Px_local*leg1eX_z_ + nu1Py_local*leg1eY_z_ + nu1Pz_local*leg1eZ_z_;
  //std::cout << "nu1: Px = " << nu1Px << ", Py = " << nu1Py << ", Pz = " << nu1Pz << std::endl;

  double tau1En = vis1En + nu1En;
  double tau1Px = vis1Px + nu1Px;
  double tau1Py = vis1Py + nu1Py;
  double tau1Pz = vis1Pz + nu1Pz;
  //std::cout << "tau1: Px = " << tau1Px << ", Py = " << tau1Py << ", Pz = " << tau1Pz << std::endl;
  
  // compute neutrino and tau lepton four-vector for second tau
  double nu2En    = vis2En*(1. - x2)/x2;
  double nu2Mass  = ( idxLeg2_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg2_mNuNu_]) : 0.;
  double nu2P     = TMath::Sqrt(TMath::Max(0., nu2En*nu2En - square(nu2Mass)));
  //std::cout << "nu2P = " << nu2P << std::endl;
  assert(idxLeg2_phi_ != -2);  
  double phiNu2 = x[idxLeg2_phi_];
  double cosThetaNu2 = compCosThetaNu(vis2En, vis2P, leg2Mass2_, nu2En, nu2P, square(nu2Mass));
  if ( !(cosThetaNu2 >= -1. && cosThetaNu2 <= +1.) ) return 0.;
  //std::cout << "cosThetaNu2 = " << cosThetaNu2 << std::endl;
  double thetaNu2 = TMath::ACos(cosThetaNu2);
  double nu2Px_local = nu2P*TMath::Cos(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Py_local = nu2P*TMath::Sin(phiNu2)*TMath::Sin(thetaNu2);
  double nu2Pz_local = nu2P*TMath::Cos(thetaNu2);
  //std::cout << "nu2_local: Px = " << nu2Px_local << ", Py = " << nu2Py_local << ", Pz = " << nu2Pz_local << std::endl;
  double nu2Px = nu2Px_local*leg2eX_x_ + nu2Py_local*leg2eY_x_ + nu2Pz_local*leg2eZ_x_;
  double nu2Py = nu2Px_local*leg2eX_y_ + nu2Py_local*leg2eY_y_ + nu2Pz_local*leg2eZ_y_;
  double nu2Pz = nu2Px_local*leg2eX_z_ + nu2Py_local*leg2eY_z_ + nu2Pz_local*leg2eZ_z_;
  //std::cout << "nu2: Px = " << nu2Px << ", Py = " << nu2Py << ", Pz = " << nu2Pz << std::endl;

  double tau2En = vis2En + nu2En;
  double tau2Px = vis2Px + nu2Px;
  double tau2Py = vis2Py + nu2Py;
  double tau2Pz = vis2Pz + nu2Pz;
  //std::cout << "tau2: Px = " << tau2Px << ", Py = " << tau2Py << ", Pz = " << tau2Pz << std::endl;

  // evaluate transfer function for MET/hadronic recoil
  double residualX = x[6];
  double residualY = x[7];
  double sumNuPx = nu1Px + nu2Px;
  double sumNuPy = nu1Py + nu2Py;
  double measuredMETx = residualX + sumNuPx;
  double measuredMETy = residualY + sumNuPy;
  double pull2 = residualX*(invCovMETxx_*residualX + invCovMETxy_*residualY) + residualY*(invCovMETyx_*residualX + invCovMETyy_*residualY);
  double prob_TF = const_MET_*TMath::Exp(-0.5*pull2);

  double prob_acceptance = 1.;
  if ( acceptance_ ) {
    prob_acceptance *= (*acceptance_)(vis1P4, vis2P4, measuredMETx, measuredMETy);
  }

  // evaluate transfer functions for tau energy reconstruction
  if ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) prob_TF *= extractProbFromLUT(x[idxLeg1VisPtShift_], leg1lutVisPtRes_);
  if ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) prob_TF *= extractProbFromLUT(x[idxLeg2VisPtShift_], leg2lutVisPtRes_);

  // compute Bjorken-x of incoming protons and evaluate PDF factor
  double ditauPx = tau1Px + tau2Px;
  double ditauPy = tau1Py + tau2Py;
  double ditauPz = tau1Pz + tau2Pz;
  double ditauEn = tau1En + tau2En;
  double ditauMass = TMath::Sqrt(square(ditauEn) - (square(ditauPx) + square(ditauPy) + square(ditauPz)));
  //if ( !(ditauMass > 0.70*mH_ && ditauMass < 1.30*mH_) ) return 0.;
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

  // perform boost into MEM frame and evaluate LO matrix element, 
  // computed by Madgraph
  LorentzVector tau1P4(tau1Px, tau1Py, tau1Pz, tau1En);
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);
  LorentzVector comSystem = tau1P4 + tau2P4;
  Vector boost = comSystem.BoostToCM();
  LorentzVector tau1P4_rf = ROOT::Math::VectorUtil::boost(tau1P4, boost);
  tau1P4_rf = LorentzVector(tau1P4_rf.px(), tau1P4_rf.py(), tau1P4_rf.pz(), TMath::Sqrt(square(tau1P4_rf.P()) + tauLeptonMass2));
  LorentzVector tau2P4_rf = ROOT::Math::VectorUtil::boost(tau2P4, boost);
  tau2P4_rf = LorentzVector(tau2P4_rf.px(), tau2P4_rf.py(), tau2P4_rf.pz(), TMath::Sqrt(square(tau2P4_rf.P()) + tauLeptonMass2));
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
  me_madgraph_.setMomenta(madgraphMomenta_);
  me_madgraph_.sigmaKin();
  //double prob_ME = me_madgraph_.getMatrixElements()[0];
  //if ( TMath::IsNaN(prob_ME) ) {
  //  std::cerr << "Warning: Magraph returned NaN --> skipping event !!" << std::endl;
  //  std::cerr << " tau1: Pt = " << tau1P4.pt() << ", eta = " << tau1P4.eta() << ", phi = " << tau1P4.phi() << ", mass = " << tau1P4.mass() << std::endl;
  //  std::cerr << " tau2: Pt = " << tau2P4.pt() << ", eta = " << tau2P4.eta() << ", phi = " << tau2P4.phi() << ", mass = " << tau2P4.mass() << std::endl;
  //  return 0.;
  //}
  me_lit_.setHiggsMass(mTest_);
  me_lit_.setHiggsWidth(GammaH_);
  me_lit_.setMomenta(madgraphMomenta_);
  double prob_ME = me_lit_.getMatrixElement();
  //std::cout << "prob_ME: madgraph = " << me_madgraph_.getMatrixElements()[0] << ", lit = " << (1./16.)*me_lit_.getMatrixElement() << std::endl;
  assert(prob_ME >= 0.);
  
  const double hbar_c = 0.1973; // GeV fm
  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  const double constFactor = conversionFactor/eigth(2.*TMath::Pi());
  double prob_PS_and_tauDecay = constFactor;
  if ( leg1isLep_ ) {
    prob_PS_and_tauDecay *= compPSfactor_tauToLepDecay(x1, vis1En, vis1P, leg1Mass_, nu1En, nu1P, nu1Mass);
  } else {
    prob_PS_and_tauDecay *= compPSfactor_tauToHadDecay(x1, vis1En, vis1P, leg1Mass_, nu1En, nu1P);
  }  
  if ( leg2isLep_ ) {
    prob_PS_and_tauDecay *= compPSfactor_tauToLepDecay(x2, vis2En, vis2P, leg2Mass_, nu2En, nu2P, nu2Mass);
  } else {
    prob_PS_and_tauDecay *= compPSfactor_tauToHadDecay(x2, vis2En, vis2P, leg2Mass_, nu2En, nu2P);
  }
  prob_PS_and_tauDecay *= (1./(4.*tau1En*tau2En));

  double jacobiFactor = ((2.*mVis2*visPtShift1*visPtShift2)/(square(q2)*x2))*(square(comSystem.M2() - mTest2_) + square(GammaH_times_mH))/GammaH_times_mH; 
  jacobiFactor *= square(vis1P)*TMath::Sin(vis1Theta);
  jacobiFactor *= square(vis2P)*TMath::Sin(vis2Theta);
  jacobiFactor *= (1./(2.*vis1En*vis2En*TMath::Max(1.e-2, 1. - compScalarProduct(eZ1, eZ2))));
  if ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) jacobiFactor *= vis1P4.pt();
  if ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) jacobiFactor *= vis2P4.pt();

  double prob = prob_flux*prob_PDF*prob_ME*prob_PS_and_tauDecay*prob_TF*prob_acceptance*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
	      << " TF = " << prob_TF << ", acceptance = " << prob_acceptance << " --> returning " << prob << std::endl;
  }
  return prob;
}
