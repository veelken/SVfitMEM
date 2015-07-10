#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandWithTauDecays.h"

#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include <TMath.h>
#include <TH1.h>
#include <TString.h>
#include <TFile.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  double xfx(int nset, double x, double Q, int fl);
}

/// global function pointer, needed for VEGAS integration
const HttXsectionIntegrandWithTauDecays* HttXsectionIntegrandWithTauDecays::gHttXsectionIntegrandWithTauDecays = 0;
bool HttXsectionIntegrandWithTauDecays::pdfIsInitialized_ = false;

HttXsectionIntegrandWithTauDecays::HttXsectionIntegrandWithTauDecays(double sqrtS, double mH, const std::string& pdfFileName, int mode, const std::string& madgraphFileName, int verbosity) 
  : mode_(mode),
    applyMEtTF_(true),
    mH_(mH),
    mH2_(mH*mH),
    GammaH_(1.e-2*mH_),
    sqrtS_(sqrtS),
    s_(square(sqrtS_)),
    invSqrtS_(1./sqrtS_),
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
    me_madgraph_(false, true),
    me_madgraph_isInitialized_(false),
    me_lit_(false, true),
    isFirstCall_(true),
    histogramLogIntegrand_(0),    
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

  histogramLogIntegrand_ = new TH1D("histogramLogIntegrand", "histogramLogIntegrand", 2000, -1.e+2, +1.e+2);
  
  // set global function pointer to this
  gHttXsectionIntegrandWithTauDecays = this;
}

HttXsectionIntegrandWithTauDecays::~HttXsectionIntegrandWithTauDecays()
{
  std::cout << "<HttXsectionIntegrandWithTauDecays::~HttXsectionIntegrandWithTauDecays>:" << std::endl;
  delete [] madgraphGluon1P4_;
  delete [] madgraphGluon2P4_;
  delete [] madgraphTau1P4_;
  delete [] madgraphTau2P4_;

  std::cout << "integrand: min = " << minIntegrand_ << ", max = " << maxIntegrand_ << std::endl;
  outputFileName_ = Form("HttXsectionIntegrandWithTauDecays_mH%03.0f", mH_);
  if ( acceptance_ ) outputFileName_.append("_wAcc");
  else outputFileName_.append("_woAcc");
  outputFileName_.append(".root");
  std::cout << "saving histogramLogIntegrand to file = " << outputFileName_ << std::endl;
  TFile* outputFile = new TFile(outputFileName_.data(), "RECREATE");
  histogramLogIntegrand_->Write();
  delete outputFile;

  delete histogramLogIntegrand_;
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
HttXsectionIntegrandWithTauDecays::setInputs(int tau1Type, double vis1Mass, int tau2Type, double vis2Mass, const TMatrixD& covMET)
{
  if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandWithTauDecays::setInputs>:" << std::endl;
  }

  leg1isLep_ = (tau1Type == MeasuredTauLepton::kTauToElecDecay || tau1Type == MeasuredTauLepton::kTauToMuDecay);
  leg1Mass_ = vis1Mass;
  leg1Mass2_ = square(leg1Mass_);
  
  leg2isLep_ = (tau2Type == MeasuredTauLepton::kTauToElecDecay || tau2Type == MeasuredTauLepton::kTauToMuDecay);
  leg2Mass_ = vis2Mass;
  leg2Mass2_ = square(leg2Mass_);

  // determine transfer matrix for MET
  double covDet = covMET.Determinant();
  const_MET_ = 0.;
  if ( covDet != 0 ) { 
    invCovMET_ = covMET;
    invCovMET_.Invert(); 
    //std::cout << "invCovMET:" << std::endl;
    //invCovMET_.Print();
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
    //std::cout << "<extractProbFromLUT>:" << std::endl;
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
    //std::cout << "<compPSfactor_tauToLepDecay>:" << std::endl;
    //std::cout << " x = " << x << std::endl;
    //std::cout << " visEn = " << visEn << std::endl;
    //std::cout << " visP = " << visP << std::endl;
    //std::cout << " visMass = " << visMass << std::endl;
    //std::cout << " nuEn = " << nuEn << std::endl;
    //std::cout << " nuP = " << nuP << std::endl;
    //std::cout << " nuMass = " << nuMass << std::endl;
    if ( x >= 0. && x <= 1. && nuMass < TMath::Sqrt((1. - x)*tauLeptonMass2) ) { // physical solution
      const double GFfactor = square(GF)/(2.*square(TMath::Pi()));
      double nuMass2 = square(nuMass);
      double visMass2 = square(visMass);
      double tauEn_rf = (tauLeptonMass2 + nuMass2 - visMass2)/(2.*nuMass);
      double visEn_rf = tauEn_rf - nuMass;
      double I = GFfactor*nuMass2*(2.*tauEn_rf*visEn_rf - (2./3.)*TMath::Sqrt((square(tauEn_rf) - tauLeptonMass2)*(square(visEn_rf) - visMass2)));
      double cosThetaNu = compCosThetaNu(visEn, visP, visMass2, nuEn, nuP, nuMass2);
      //std::cout << "cosThetaNu = " << cosThetaNu << std::endl;
      if ( !(cosThetaNu >= -1. && cosThetaNu <= +1.) ) return 0.;
      double PSfactor = (visEn + nuEn)*I/(8.*visP*square(x)*TMath::Sqrt(square(visP) + square(nuP) + 2.*visP*nuP*cosThetaNu + tauLeptonMass2));
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
      double cosThetaNu = compCosThetaNu(visEn, visP, visMass2, nuEn, nuP, 0.);
      //std::cout << "cosThetaNu = " << cosThetaNu << std::endl;
      if ( !(cosThetaNu >= -1. && cosThetaNu <= +1.) ) return 0.;
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

  void fillWithOverFlow(TH1* histogram, double x)
  {
    TAxis* xAxis = histogram->GetXaxis();
    int bin = xAxis->FindBin(x);
    int numBins = xAxis->GetNbins();
    if ( bin < 1       ) bin = 1;
    if ( bin > numBins ) bin = numBins;
    double x_bin = histogram->GetBinCenter(bin);
    histogram->Fill(x_bin);
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
  //std::cout << "visPtShift1 = " << visPtShift1 << ", visPtShift2 = " << visPtShift2 << std::endl;

  // compute four-vector of visible decay products for first tau
  //-----------------------------------------------------------------------------
  // CV: build vis1P4 for the case that integral over momenta of visible tau decay products
  //     is computed using Cartesian coordinates
  double vis1Px    = visPtShift1*x[0];
  double vis1Py    = visPtShift1*x[1];
  double vis1Pz    = visPtShift1*x[2];
  double vis1En    = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);
  double vis1P     = vis1P4.P();
  double vis1Theta = vis1P4.theta();
  double vis1Phi   = vis1P4.phi();
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // CV: build vis1P4 for the case that integral over momenta of visible tau decay products
  //     is computed using polar coordinates
  //double vis1P     = visPtShift1*x[0];
  //double vis1Theta = x[1];
  //double vis1Phi   = x[2];
  //double vis1Px    = vis1P*TMath::Cos(vis1Phi)*TMath::Sin(vis1Theta);
  //double vis1Py    = vis1P*TMath::Sin(vis1Phi)*TMath::Sin(vis1Theta);
  //double vis1Pz    = vis1P*TMath::Cos(vis1Theta);
  //double vis1En    = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  //LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);
  //-----------------------------------------------------------------------------
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

  // compute four-vector of visible decay products for second tau
  //-----------------------------------------------------------------------------
  // CV: build vis2P4 for the case that integral over momenta of visible tau decay products
  //     is computed using Cartesian coordinates
  //double vis2Px    = visPtShift2*x[3];
  //double vis2Py    = visPtShift2*x[4];
  //double vis2Pz    = visPtShift2*x[5];
  //double vis2En    = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  //LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);
  //double vis2P     = vis2P4.P();
  //double vis2Theta = vis2P4.theta();
  //double vis2Phi   = vis2P4.phi();
  //Vector eZ2 = Vector(TMath::Cos(vis2Phi)*TMath::Sin(vis2Theta), TMath::Sin(vis2Phi)*TMath::Sin(vis2Theta), TMath::Cos(vis2Theta));
  //Vector eY2 = normalize(compCrossProduct(eZ2, beamAxis_));
  //Vector eX2 = normalize(compCrossProduct(eY2, eZ2));
  //double mVis = (vis1P4 + vis2P4).mass();
  //double mVis2 = square(mVis);
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // CV: build vis2P4 for the case that integral over momenta of visible tau decay products
  //     is computed using polar coordinates
  double mVis2     = x[3];
  double mVis      = TMath::Sqrt(mVis2);
  double vis2Theta = x[4];
  double vis2Phi   = x[5];
  Vector eZ2 = Vector(TMath::Cos(vis2Phi)*TMath::Sin(vis2Theta), TMath::Sin(vis2Phi)*TMath::Sin(vis2Theta), TMath::Cos(vis2Theta));
  Vector eY2 = normalize(compCrossProduct(eZ2, beamAxis_));
  Vector eX2 = normalize(compCrossProduct(eY2, eZ2));
  double vis2P     = mVis2/(2.*vis1P*TMath::Max(1.e-2, 1. - compScalarProduct(eZ1, eZ2)));
  double vis2Px    = vis2P*TMath::Cos(vis2Phi)*TMath::Sin(vis2Theta);
  double vis2Py    = vis2P*TMath::Sin(vis2Phi)*TMath::Sin(vis2Theta);
  double vis2Pz    = vis2P*TMath::Cos(vis2Theta);
  double vis2En    = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);
  //-----------------------------------------------------------------------------
  leg2eX_x_ = eX2.x();
  leg2eX_y_ = eX2.y();
  leg2eX_z_ = eX2.z();
  leg2eY_x_ = eY2.x();
  leg2eY_y_ = eY2.y();
  leg2eY_z_ = eY2.z();
  leg2eZ_x_ = eZ2.x();
  leg2eZ_y_ = eZ2.y();
  leg2eZ_z_ = eZ2.z();

  assert(idxLeg2_X_ != -1);
  double x2 = x[idxLeg2_X_];
  if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;
  
  double GammaH_times_mH = GammaH_*mH_;
  assert(idxLeg1_t_ != -1);  
  double tk = x[idxLeg1_t_];
  double q2 = mH2_ + GammaH_times_mH*TMath::Tan(tk);
  if ( q2 <= 0. ) return 0.;
  double x1 = mVis2/(q2*x2);
  if ( !(x1 >= 1.e-5 && x1 <= 1.) ) return 0.;

  // compute neutrino and tau lepton four-vector for first tau
  double nu1En    = vis1En*(1. - x1)/x1;
  double nu1Mass  = ( idxLeg1_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg1_mNuNu_]) : 0.;
  double nu1P     = TMath::Sqrt(TMath::Max(0., nu1En*nu1En - square(nu1Mass)));
  assert(idxLeg1_phi_ != -1);  
  double phiNu1 = x[idxLeg1_phi_];
  double cosThetaNu1 = compCosThetaNu(vis1En, vis1P, leg1Mass2_, nu1En, nu1P, square(nu1Mass));
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
  double tau1Px = vis1Px + nu1Px;
  double tau1Py = vis1Py + nu1Py;
  double tau1Pz = vis1Pz + nu1Pz;
  LorentzVector tau1P4(tau1Px, tau1Py, tau1Pz, tau1En);

  // compute neutrino and tau lepton four-vector for second tau
  double nu2En    = vis2En*(1. - x2)/x2;
  double nu2Mass  = ( idxLeg2_mNuNu_ != -1 ) ? TMath::Sqrt(x[idxLeg2_mNuNu_]) : 0.;
  double nu2P     = TMath::Sqrt(TMath::Max(0., nu2En*nu2En - square(nu2Mass)));
  assert(idxLeg2_phi_ != -2);  
  double phiNu2 = x[idxLeg2_phi_];
  double cosThetaNu2 = compCosThetaNu(vis2En, vis2P, leg2Mass2_, nu2En, nu2P, square(nu2Mass));
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
  double tau2Px = vis2Px + nu2Px;
  double tau2Py = vis2Py + nu2Py;
  double tau2Pz = vis2Pz + nu2Pz;
  LorentzVector tau2P4(tau2Px, tau2Py, tau2Pz, tau2En);

  double prob_TF = 1.;

  // evaluate transfer function for MET/hadronic recoil
  double residualX = x[6];
  double residualY = x[7];
  double sumNuPx = nu1Px + nu2Px;
  double sumNuPy = nu1Py + nu2Py;
  double measuredMETx = residualX + sumNuPx;
  double measuredMETy = residualY + sumNuPy;
  if ( applyMEtTF_ ) {
    double pull2 = residualX*(invCovMETxx_*residualX + invCovMETxy_*residualY) + residualY*(invCovMETyx_*residualX + invCovMETyy_*residualY);
    //std::cout << "residualX = " << residualX << ", residualY = " << residualY << ": pull2 = " << pull2 << std::endl;
    prob_TF *= (const_MET_*TMath::Exp(-0.5*pull2));
  }

  double prob_acceptance = 1.;
  if ( acceptance_ ) {
    prob_acceptance *= (*acceptance_)(vis1P4, vis2P4, measuredMETx, measuredMETy);
  }

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
  //if ( !(ditauMass > 0.70*mH_ && ditauMass < 1.30*mH_) ) return 0.;
  //std::cout << "ditau: En = " << ditauEn << ", Pz = " << ditauPz << std::endl;
  double xa = invSqrtS_*(ditauEn + ditauPz);
  double xb = invSqrtS_*(ditauEn - ditauPz);   
  //std::cout << "xa = " << xa << ", xb = " << xb << std::endl;
  if ( xa <= 0. || xa >= 1. ) return 0.;
  if ( xb <= 0. || xb >= 1. ) return 0.;
  //double Q = mH_;
  double Q = ditauMass;
  assert(pdfIsInitialized_);
  double fa = LHAPDF::xfx(1, xa, Q, 0)/xa;
  double fb = LHAPDF::xfx(1, xb, Q, 0)/xb;
  double prob_PDF = (fa*fb);

  // evaluate flux factor
  double prob_flux = (1./(s_*xa*xb));  

  // evaluate LO matrix element, 
  // computed by Madgraph
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
  
  // CV: The LO cross-section assumes that the Higgs has zero zero transverse momentum,
  //     corresponding to a delta-function 
  //      delta(tau1Px - tau2Px) * delta(tau1Px - tau2Px) 
  //    = delta(visPtShift1*vis1Px/x1 - visPtShift2*vis2Px/x2) * delta(visPtShift1*vis1Py/x1 - visPtShift2*vis2Py/x2)
  //    = (x2/visPtShift2)^2 * delta(visPtShift1*vis1Px*x2/(x1*visPtShift2) - vis2Px) * delta(visPtShift1*vis1Py*x2/(x1*visPtShift2) - vis2Py)
  //     where the delta-function rule: delta(alpha x) = 1/|alpha| * delta(x) has been used, cf. https://en.wikipedia.org/wiki/Dirac_delta_function 
  //
  //     Instead of including the delta-functions into the numeric integration,
  //     we correct the integral by the following means:
  //      1) in HttXsectionIntegrandWithTauDecays::Eval we multiply the integrand by the factor (x2/visPtShift2)^2
  //      2) in HttXsectionWithTauDecays::integrate we multiply the value of the integral by a factor 1/((xu_vis2Px - xl_vis2Px)*(xu_vis2Py - xl_vis2Py))
  //
  //     Note that as we boost all four-vectors into the MEM frame, the integrand is constant for all values of vis2Px and vis2Py
  //
  //-----------------------------------------------------------------------------
  // CV: factor for the case that integral over momenta of visible tau decay products
  //     is computed using Cartesian coordinates
  double memFrameFactor = ( visPtShift1 > 0. ) ? square(x1/visPtShift1) : 0.;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // CV: factor for the case that integral over momenta of visible tau decay products
  //     is computed using Cartesian coordinates
  //double memFrameFactor = square(2.*vis1P*TMath::Max(1.e-2, 1. - compScalarProduct(eZ1, eZ2))*x2)/(TMath::Cos(vis2Phi)*TMath::Sin(vis2Theta)*TMath::Sin(vis2Phi)*TMath::Sin(vis2Theta));
  //-----------------------------------------------------------------------------

  const double conversionFactor = 1.e+10*square(hbar_c); // conversion factor from GeV^-2 to picobarn = 10^-40m
  const double constFactor = 2.*conversionFactor/sixth(2.*TMath::Pi());
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

  double mVis_mem = (vis1P4_mem + vis2P4_mem).mass();
  double jacobiFactor = (square(mVis_mem)/(square(q2)*x2_mem))*(square(q2 - mH2_) + square(GammaH_times_mH))/GammaH_times_mH;
  //-----------------------------------------------------------------------------
  // CV: extra Jacobi factors if integrating over momenta of visible tau decay products
  //     using polar coordinates
  //jacobiFactor *= square(vis1P)*TMath::Sin(vis1Theta);
  jacobiFactor *= square(vis2P)*TMath::Sin(vis2Theta);
  jacobiFactor *= (1./(2.*vis1P*TMath::Max(1.e-2, 1. - compScalarProduct(eZ1, eZ2))));
  //-----------------------------------------------------------------------------
  //if ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) jacobiFactor *= vis1P4.pt();
  //if ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) jacobiFactor *= vis2P4.pt();

  double prob = prob_flux*prob_PDF*prob_ME*memFrameFactor*prob_PS_and_tauDecay*prob_TF*prob_acceptance*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", MEM-frame = " << memFrameFactor << ", PS+decay = " << prob_PS_and_tauDecay << "," 
	      << " TF = " << prob_TF << ", acceptance = " << prob_acceptance << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
  }

  if ( prob < minIntegrand_ || isFirstCall_ ) {
    minIntegrand_ = prob;
  }
  if ( prob > maxIntegrand_ || isFirstCall_ ) {
    maxIntegrand_ = prob;
  }
  if ( prob > 0. ) {
    fillWithOverFlow(histogramLogIntegrand_, TMath::Log10(prob));
  }
  isFirstCall_ = false;
  //if ( prob < 0. || prob > 1.e+6 ) {
  //  std::cout << "CHECK prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", MEM-frame = " << memFrameFactor << ", PS+decay = " << prob_PS_and_tauDecay << "," 
  //	        << " TF = " << prob_TF << ", acceptance = " << prob_acceptance << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
  //  std::cout << "xa = " << xa << ", xb = " << xb << ": q2 = " << q2 << std::endl;
  //  std::cout << "vis1: Pt = " << vis1P4.pt() << ", eta = " << vis1P4.eta() << ", phi = " << vis1P4.phi() << ", mass = " << vis1P4.mass() << std::endl;
  //  std::cout << "vis2: Pt = " << vis2P4.pt() << ", eta = " << vis2P4.eta() << ", phi = " << vis2P4.phi() << ", mass = " << vis2P4.mass() << std::endl;
  //  std::cout << "visPtShift: leg1 = " << visPtShift1 << ", leg2 = " << visPtShift2 << std::endl;
  //  std::cout << "mVis = " << mVis << " (in MEM-frame = " << mVis_mem << ")" << std::endl;
  //  std::cout << "compScalarProduct(eZ1, eZ2) = " << compScalarProduct(eZ1, eZ2) << std::endl;
  //  std::cout << "x1 = " << x1 << " (in MEM-frame = " << x1_mem << "), x2 = " << x2 << " (in MEM-frame = " << x2_mem << ")" << std::endl;
  //  std::cout << "prob_tauDecay: leg1 = " << prob_tauDecay_leg1 << ", leg2 = " << prob_tauDecay_leg2 << std::endl;
  //  std::cout << "ditau: En = " << ditauEn << ", Pz = " << ditauPz << std::endl;
  //}

  return prob;
}
