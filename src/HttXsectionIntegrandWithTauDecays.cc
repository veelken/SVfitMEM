#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandWithTauDecays.h"

#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include <TMath.h>
#include <TH1.h>
#include <TString.h>
#include <TFile.h>
#include <Math/VectorUtil.h>

using namespace svFitMEM;

/// global function pointer, needed for VEGAS integration
const HttXsectionIntegrandWithTauDecays* HttXsectionIntegrandWithTauDecays::gHttXsectionIntegrandWithTauDecays = 0;
int HttXsectionIntegrandWithTauDecays::gNumInstances = 0;

HttXsectionIntegrandWithTauDecays::HttXsectionIntegrandWithTauDecays(double sqrtS, double mH, const std::string& pdfName, int mode, const std::string& madgraphFileName, int verbosity) 
  : mode_(mode),
    applyMEtTF_(true),    
    mH_(mH),
    mH2_(mH*mH),
    GammaH_(1.e-2*mH_),
    GammaH_times_mH_(GammaH_*mH_),
    GammaH2_times_mH2_(square(GammaH_times_mH_)),
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
    pdf_(0),
    pdfIsInitialized_(false),
    me_madgraph_(false, true),
    me_madgraph_isInitialized_(false),
    me_lit_(pdfName, false, true),
    isFirstCall_(true),
    histogramLogIntegrand_(0),    
    verbosity_(verbosity)
{
  if ( verbosity_ ) {
    std::cout << "<HttXsectionIntegrandWithTauDecays::HttXsectionIntegrandWithTauDecays>:" << std::endl;
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
    std::cerr << "Error in <HttXsectionIntegrandWithTauDecays>: No param.dat file for Madgraph given !!" << std::endl;
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

  std::string histogramNameLogIntegrand = Form("histogramLogIntegrand%i", gNumInstances);
  histogramLogIntegrand_ = new TH1D(histogramNameLogIntegrand.data(), histogramNameLogIntegrand.data(), 2000, -1.e+2, +1.e+2);
  
  // set global function pointer to this
  gHttXsectionIntegrandWithTauDecays = this;
  ++gNumInstances;
}

HttXsectionIntegrandWithTauDecays::~HttXsectionIntegrandWithTauDecays()
{
  std::cout << "<HttXsectionIntegrandWithTauDecays::~HttXsectionIntegrandWithTauDecays>:" << std::endl;
  delete pdf_;

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

void HttXsectionIntegrandWithTauDecays::updateNumDimensions()
{
  numDimensions_ = -1;
  if ( idxLeg1_X_         >= numDimensions_ ) numDimensions_ = idxLeg1_X_;
  if ( idxLeg1_phi_       >= numDimensions_ ) numDimensions_ = idxLeg1_phi_;
  if ( idxLeg1VisPtShift_ >= numDimensions_ ) numDimensions_ = idxLeg1VisPtShift_;
  if ( idxLeg1_mNuNu_     >= numDimensions_ ) numDimensions_ = idxLeg1_mNuNu_;
  if ( idxLeg2_t_         >= numDimensions_ ) numDimensions_ = idxLeg2_t_;
  if ( idxLeg2_phi_       >= numDimensions_ ) numDimensions_ = idxLeg2_phi_;
  if ( idxLeg2VisPtShift_ >= numDimensions_ ) numDimensions_ = idxLeg2VisPtShift_;
  if ( idxLeg2_mNuNu_     >= numDimensions_ ) numDimensions_ = idxLeg2_mNuNu_;
  ++numDimensions_; // number of dimensions = max(idx) + 1
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
  double compDgDp2z(double vis1Pz, double vis1P, double vis2Pz, double vis2P)
  {
    if ( vis2P <= 0. ) return 0.;
    else return (2.*vis1P*vis2Pz - 2.*vis1Pz*vis2P)/vis2P;
  }

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

  template <typename T>
  std::string format_T_array(const T* array, unsigned d)
  {
    std::ostringstream os;
    
    os << "{ ";
    
    for ( unsigned i = 0; i < d; ++i ) {
      os << array[i];
      if ( i < (d - 1) ) os << ", ";
    }
    
    os << " }";
    
    return os.str();
  }
  
  std::string format_double_array(const double* array, unsigned d)
  {
    return format_T_array(array, d);
  }
}

double
HttXsectionIntegrandWithTauDecays::Eval(const double* x) const 
{
  if ( verbosity_ >= 2 ) {
    std::cout << "<HttXsectionIntegrandWithTauDecays::Eval(const double*)>:" << std::endl;
    std::cout << " x = " << format_double_array(x, numDimensions_) << std::endl;
  }

  double visPtShift1 = ( idxLeg1VisPtShift_ != -1 && !leg1isLep_ ) ? 1./x[idxLeg1VisPtShift_] : 1.;
  //double visPtShift2 = ( idxLeg2VisPtShift_ != -1 && !leg2isLep_ ) ? 1./x[idxLeg2VisPtShift_] : 1.;
  //std::cout << "visPtShift1 = " << visPtShift1 << ", visPtShift2 = " << visPtShift2 << std::endl;

  // compute four-vector of visible decay products for first tau
  double vis1Px = visPtShift1*x[0];
  double vis1Py = visPtShift1*x[1];
  double vis1Pz = visPtShift1*TMath::Power(x[2], 3.);
  double vis1En = TMath::Sqrt(square(vis1Px) + square(vis1Py) + square(vis1Pz) + leg1Mass2_);
  LorentzVector vis1P4(vis1Px, vis1Py, vis1Pz, vis1En);
  double vis1P = vis1P4.P();
  double vis1Theta = vis1P4.theta();
  double vis1Phi = vis1P4.phi();
  eZ1_ = Vector(TMath::Cos(vis1Phi)*TMath::Sin(vis1Theta), TMath::Sin(vis1Phi)*TMath::Sin(vis1Theta), TMath::Cos(vis1Theta));
  eY1_ = normalize(compCrossProduct(eZ1_, beamAxis_));
  eX1_ = normalize(compCrossProduct(eY1_, eZ1_));
  leg1eX_x_ = eX1_.x();
  leg1eX_y_ = eX1_.y();
  leg1eX_z_ = eX1_.z();
  leg1eY_x_ = eY1_.x();
  leg1eY_y_ = eY1_.y();
  leg1eY_z_ = eY1_.z();
  leg1eZ_x_ = eZ1_.x();
  leg1eZ_y_ = eZ1_.y();
  leg1eZ_z_ = eZ1_.z();

  // compute four-vector of visible decay products for second tau
  assert(idxLeg1_X_ != -1);
  double x1 = x[idxLeg1_X_];
  if ( !(x1 >= 1.e-3 && x1 <= 1.) ) return 0.;
  
  assert(idxLeg2_t_ != -1);  
  double tk = x[idxLeg2_t_];
  q2_ = mH2_ + GammaH_times_mH_*TMath::Tan(tk);
  if ( q2_ <= 0. ) return 0.;

  double mVis2 = x[3];
  //if ( mVis2 < 1.e-3*mH2_ ) return 0.;
  double x2 = mVis2/(q2_*x1);
  if ( !(x2 >= 1.e-3 && x2 <= 1.) ) return 0.;

  double vis2Px = -x2*vis1Px/x1;
  double vis2Py = -x2*vis1Py/x1;
  double vis1Pt2 = square(vis1P4.pt());
  if ( verbosity_ >= 3 ) {
    std::cout << "vis1Pt = " << TMath::Sqrt(vis1Pt2) << std::endl;
  }
  double term1 = 2.*vis1Pt2;
  double term2 = vis1Pz*(mVis2 - 2.*(x2/x1)*vis1Pt2);
  double vis1P2 = square(vis1P);
  double term3_2 = vis1P2*mVis2*(mVis2 - 4.*(x2/x1)*vis1Pt2);
  if ( verbosity_ >= 3 ) {
    std::cout << "term1 = " << term1 << ", term3_2 = " << term3_2 << std::endl;
  }
  if ( term3_2 <= 0. || term1 <= 0. ) return 0.;
  double term3 = TMath::Sqrt(term3_2);
  double vis2Pz_p = (1./term1)*(term2 + term3);
  double vis2En_p = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz_p) + leg2Mass2_);
  LorentzVector vis2P4_p(vis2Px, vis2Py, vis2Pz_p, vis2En_p);
  double absDgDp2z_p = TMath::Abs(compDgDp2z(vis1Pz, vis1P, vis2Pz_p, vis2P4_p.P()));
  double vis2Pz_m = (1./term1)*(term2 - term3);
  double vis2En_m = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz_m) + leg2Mass2_);
  LorentzVector vis2P4_m(vis2Px, vis2Py, vis2Pz_m, vis2En_m);
  double absDgDp2z_m = TMath::Abs(compDgDp2z(vis1Pz, vis1P, vis2Pz_m, vis2P4_m.P()));
  double prob = 0.;
  if ( verbosity_ >= 3 ) {
    std::cout << " vis2Pz_p = " << vis2Pz_p << ", absDgDp2z_p = " << absDgDp2z_p << " (sqrtS = " << sqrtS_ << ")" << std::endl;
  }
  if ( TMath::Abs(vis2Pz_p) <= sqrtS_ && absDgDp2z_p != 0. ) { 
    double prob_p = compProb(x, vis1P4, x1, vis2P4_p, x2, absDgDp2z_p);
    if ( verbosity_ >= 3 ) {
      std::cout << "solution+: vis2Pz = " << vis2Pz_p << ", prob = " << prob_p << std::endl;
    }
    prob += prob_p;
  }
  if ( verbosity_ >= 3 ) {
    std::cout << " vis2Pz_m = " << vis2Pz_m << ", absDgDp2z_m = " << absDgDp2z_m << std::endl;
  }
  if ( TMath::Abs(vis2Pz_m) <= sqrtS_ && absDgDp2z_m != 0. ) { 
    double prob_m = compProb(x, vis1P4, x1, vis2P4_m, x2, absDgDp2z_m);
    if ( verbosity_ >= 3 ) {
      std::cout << "solution-: vis2Pz = " << vis2Pz_m << ", prob = " << prob_m << std::endl;
    }
    prob += prob_m;
  }

  //double vis2Px = -vis1Px; // CV: this is only approximate, assuming that pT balance is between visible tau decay products and not between taus
  //double vis2Py = -vis1Py;
  //double vis2Pz = x[3];
  //double vis2En = TMath::Sqrt(square(vis2Px) + square(vis2Py) + square(vis2Pz) + leg2Mass2_);
  //LorentzVector vis2P4(vis2Px, vis2Py, vis2Pz, vis2En);
  //
  //double mVis2 = (vis1P4 + vis2P4).mass();
  //double x2 = mVis2/(q2_*x1);
  //if ( !(x2 >= 1.e-5 && x2 <= 1.) ) return 0.;
  //
  //double prob = compProb(x, vis1P4, x1, vis2P4, x2, 1.);

  return prob;
}

double
HttXsectionIntegrandWithTauDecays::compProb(const double* x, const LorentzVector& vis1P4, double x1, const LorentzVector& vis2P4, double x2, double absDgDp2z) const
{
  //std::cout << "<compProb>:" << std::endl;
  //std::cout << " absDgDp2z = " << absDgDp2z << std::endl;

  double vis1En = vis1P4.E();
  double vis1P = vis1P4.P();

  double vis2En = vis2P4.E();
  double vis2P = vis2P4.P();
  double vis2Theta = vis2P4.theta();
  double vis2Phi = vis2P4.phi();
  eZ2_ = Vector(TMath::Cos(vis2Phi)*TMath::Sin(vis2Theta), TMath::Sin(vis2Phi)*TMath::Sin(vis2Theta), TMath::Cos(vis2Theta));
  eY2_ = normalize(compCrossProduct(eZ2_, beamAxis_));
  eX2_ = normalize(compCrossProduct(eY2_, eZ2_));
  leg2eX_x_ = eX2_.x();
  leg2eX_y_ = eX2_.y();
  leg2eX_z_ = eX2_.z();
  leg2eY_x_ = eY2_.x();
  leg2eY_y_ = eY2_.y();
  leg2eY_z_ = eY2_.z();
  leg2eZ_x_ = eZ2_.x();
  leg2eZ_y_ = eZ2_.y();
  leg2eZ_z_ = eZ2_.z();

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

  double prob_TF = 1.;

  // evaluate transfer function for MET/hadronic recoil
  double residualX = x[4];
  double residualY = x[5];
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
  double fa = pdf_->xfxQ(21, xa, Q)/xa; // gluon distribution
  double fb = pdf_->xfxQ(21, xb, Q)/xb;
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

  double jacobiFactor = 1./absDgDp2z;
  double mVis2_mem = square((vis1P4_mem + vis2P4_mem).mass());
  double q4 = square(q2_);
  jacobiFactor *= (mVis2_mem/(q4*x1_mem));                                    // transformation from x2 to q^2
  jacobiFactor *= (square(q2_ - mH2_) + GammaH2_times_mH2_)/GammaH_times_mH_; // parametrization of q^2 by tk (Eq. 8 of arXiv:1010.2263 without factor 1/pi, as agreed with Luca and Andrew)
  jacobiFactor *= square(x2_mem);                                             // from delta-functions (tau1Px - tau2Px)*(tau1Py - tau2Py) = (vis1Px/x1 - vis2Px/x2)*(vis1Py/x1 - vis2Py/x2) 
  jacobiFactor *= 3.*square(x[2]);                                            // CV: reparametrization of vis1Pz by vis1Pz^3
  double prob = prob_flux*prob_PDF*prob_ME*prob_PS_and_tauDecay*prob_TF*prob_acceptance*jacobiFactor;
  if ( verbosity_ >= 2 ) {
    std::cout << "prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
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
  //  std::cout << "CHECK prob: flux = " << prob_flux << ", PDF = " << prob_PDF << ", ME = " << prob_ME << ", PS+decay = " << prob_PS_and_tauDecay << "," 
  //	        << " TF = " << prob_TF << ", acceptance = " << prob_acceptance << ", Jacobi = " << jacobiFactor << " --> returning " << prob << std::endl;
  //  std::cout << "xa = " << xa << ", xb = " << xb << ": q2 = " << q2_ << std::endl;
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
