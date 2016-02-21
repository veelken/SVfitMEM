#include "TauAnalysis/SVfitMEM/interface/me_ggH_nlo_lit.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

const double mtop = 172.5;
const double mtop2 = svFitMEM::square(mtop);

const double N = 4.;
const double C_A = N;
const double C_F = (N*N - 1.)/(2.*N);

using namespace svFitMEM;

me_ggH_nlo_lit::me_ggH_nlo_lit(const std::string& pdfName, bool applyNWA, bool includeHtoTauTauDecay)
  : applyNWA_(applyNWA),
    includeHtoTauTauDecay_(includeHtoTauTauDecay),
    s_(0.),
    s_isInitialized_(false),
    mH_(0.),
    mH2_(0.),
    mH8_(0.),
    mH_isInitialized_(false),
    GammaH_(0.),
    GammaH2_(0.),
    GammaH_isInitialized_(false),
    br_(0.),
    br_isInitialized_(false),
    pdf_(0),
    pdfIsInitialized_(false),
    mandelstam_s_(0.),
    mandelstam_abs_sMin_(5.),
    mandelstam_t_(0.),
    mandelstam_abs_tMin_(5.),
    mandelstam_u_(0.),
    mandelstam_abs_uMin_(5.)
{
  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    pdf_ = LHAPDF::mkPDF(pdfName.data(), 0);
    pdfIsInitialized_ = true;
  }

  K_gg_ = 1./(4.*square(N*N - 1.));
  K_qq_ = 1./(4.*N*N);
  K_qg_ = 1./(4.*N*(N*N - 1.));
}
 
me_ggH_nlo_lit::~me_ggH_nlo_lit() 
{
  delete pdf_;
}

void me_ggH_nlo_lit::setS(double s) 
{ 
  s_ = s;
  s_isInitialized_ = true;
}

void me_ggH_nlo_lit::setHiggsMass(double mH) 
{ 
  mH_ = mH;
  mH2_ = square(mH_);
  mH8_ = fourth(mH2_);
  mH_isInitialized_ = true;
}

void me_ggH_nlo_lit::setHiggsWidth(double width) 
{
  GammaH_ = width;
  GammaH2_ = square(GammaH_);
  GammaH_isInitialized_ = true;
}

void me_ggH_nlo_lit::setBR(double br) 
{
  br_ = br;
  br_isInitialized_ = true;
}

void me_ggH_nlo_lit::setMomenta(std::vector<const LorentzVector*>& momenta)
{ 
  //std::cout << "<me_ggH_nlo_lit::setMomenta>:" << std::endl;

  momenta_ = momenta;

  assert(momenta_.size() == 5);
  const LorentzVector& parton1P4_incoming = (*momenta_[0]);
  const LorentzVector& parton2P4_incoming = (*momenta_[1]);
  const LorentzVector& tau1P4_outgoing = (*momenta_[2]);
  const LorentzVector& tau2P4_outgoing = (*momenta_[3]);
  LorentzVector higgsP4_outgoing = tau1P4_outgoing + tau2P4_outgoing;
  //const LorentzVector& jetP4_outgoing = (*momenta_[4]);
  
  q_ = (tau1P4_outgoing + tau2P4_outgoing).mass();
  q2_ = square(q_);
  //std::cout << " q2 = " << q2_ << std::endl;

  double Q2 = q2_;
  assert(pdfIsInitialized_);
  double alphaS = pdf_->alphasQ2(Q2);
  //std::cout << " alphaS = " << alphaS << std::endl;

  const double GF = 1.166e-5; // in units of GeV^-2, taken from http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf
  double tau = 4.*mtop2/q2_;
  double tau2 = tau*tau;
  double Re_f = 0.;
  double Im_f = 0.;
  if ( tau < 1. ) {
    double beta = TMath::Sqrt(1. - tau);
    double logBetaFactor = TMath::Log((1. + beta)/(1. - beta));
    const double Pi2 = square(TMath::Pi());
    Re_f = -0.25*(square(logBetaFactor) - Pi2);
    Im_f = 0.50*TMath::Pi()*logBetaFactor;
  } else {
    double arg = 1./TMath::Sqrt(tau);
    assert(arg >= -1. && arg <= +1);
    Re_f = square(TMath::ASin(arg));
    Im_f = 0.;
  }
  double F2 = square(1. + (1. - tau)*Re_f) + square((1. - tau)*Im_f);
  const double four_Pi = 4.*TMath::Pi();
  double C = 1. + (alphaS/four_Pi)*11.;
  double C2 = C*C;
  G2_ = 4.*TMath::Sqrt(2.)*square(alphaS/four_Pi)*GF*tau2*F2*C2;
  
  double g = 4.*TMath::Pi()*alphaS;
  g2_ = g*g;

  mandelstam_s_ = (parton1P4_incoming + parton2P4_incoming).mass();
  double mandelstam_abs_s = TMath::Abs(mandelstam_s_);
  if ( mandelstam_abs_s < mandelstam_abs_sMin_ ) mandelstam_s_ = mandelstam_abs_sMin_*TMath::Sign(1., mandelstam_s_);
  mandelstam_s2_ = mandelstam_s_*mandelstam_s_;
  mandelstam_s4_ = mandelstam_s2_*mandelstam_s2_;
  mandelstam_s8_ = mandelstam_s4_*mandelstam_s4_;
  mandelstam_t_ = (parton1P4_incoming - higgsP4_outgoing).mass();
  double mandelstam_abs_t = TMath::Abs(mandelstam_t_);
  if ( mandelstam_abs_t < mandelstam_abs_tMin_ ) mandelstam_t_ = mandelstam_abs_tMin_*TMath::Sign(1., mandelstam_t_);
  mandelstam_t2_ = mandelstam_t_*mandelstam_t_;
  mandelstam_t4_ = mandelstam_t2_*mandelstam_t2_;
  mandelstam_t8_ = mandelstam_t4_*mandelstam_t4_;
  mandelstam_u_ = (parton2P4_incoming - higgsP4_outgoing).mass();
  double mandelstam_abs_u = TMath::Abs(mandelstam_u_);
  if ( mandelstam_abs_u < mandelstam_abs_uMin_ ) mandelstam_u_ = mandelstam_abs_uMin_*TMath::Sign(1., mandelstam_u_);
  mandelstam_u2_ = mandelstam_u_*mandelstam_u_;
  mandelstam_u4_ = mandelstam_u2_*mandelstam_u2_;
  mandelstam_u8_ = mandelstam_u4_*mandelstam_u4_;
  //std::cout << "Mandelstam variables: s = " << mandelstam_s_ << ", t = " << mandelstam_t_ << ", u = " << mandelstam_u_ << std::endl; 
}

double me_ggH_nlo_lit::getMatrixElement_gg() const
{
  //std::cout << "<me_ggH_lit::getMatrixElement_gg>:" << std::endl;
  double me = getMatrixElement_gg_woHtoTauTauDecay();
  me *= getMatrixElement_HtoTauTauDecay();
  //std::cout << "me = " << me << std::endl;
  return me;
}

double me_ggH_nlo_lit::getMatrixElement_gg_woHtoTauTauDecay() const
{
  double const_factor = N*(N*N - 1.);
  double me = K_gg_*G2_*g2_*const_factor*(1./(mandelstam_s_*mandelstam_t_*mandelstam_u_))*(mandelstam_s4_ + mandelstam_t4_ + mandelstam_u4_ + mH8_);
  //std::cout << "<me_ggH_nlo_lit::getMatrixElement_gg_woHtoTauTauDecay>: me = " << me << std::endl;
  return me;
}

double me_ggH_nlo_lit::getMatrixElement_qq() const
{
  //std::cout << "<me_ggH_lit::getMatrixElement_qq>:" << std::endl;
  double me = getMatrixElement_gg_woHtoTauTauDecay();
  me *= getMatrixElement_HtoTauTauDecay();
  //std::cout << "me = " << me << std::endl;
  return me;
}

double me_ggH_nlo_lit::getMatrixElement_qq_woHtoTauTauDecay() const
{
  double const_factor = C_A*C_F;
  double me = K_qq_*G2_*g2_*const_factor*(1./mandelstam_s_)*(mandelstam_t2_ + mandelstam_u2_);
  return me;
}

double me_ggH_nlo_lit::getMatrixElement_qg() const
{
  //std::cout << "<me_ggH_lit::getMatrixElement_qg>:" << std::endl;
  double me = getMatrixElement_qg_woHtoTauTauDecay();
  me *= getMatrixElement_HtoTauTauDecay();
  //std::cout << "me = " << me << std::endl;
  return me;
}

double me_ggH_nlo_lit::getMatrixElement_qg_woHtoTauTauDecay() const
{
  double const_factor = C_A*C_F;
  double me = K_qg_*G2_*g2_*const_factor*(1./mandelstam_u_)*(-(mandelstam_s2_ + mandelstam_t2_));
  //std::cout << "<me_ggH_nlo_lit::getMatrixElement_qg_woHtoTauTauDecay>: me = " << me << std::endl;
  return me;
}

double me_ggH_nlo_lit::getMatrixElement_HtoTauTauDecay() const
{
  if ( !(mH_isInitialized_ && GammaH_isInitialized_) ) {
    std::cerr << "Error in <me_ggH_nlo_lit::getMatrixElement_HtoTauTauDecay>: Higgs mass and width have not been initialized !!" << std::endl;
    assert(0);
  }
  double me = 1.;
  const double one_over_Pi = 1./TMath::Pi();
  double GammaH_times_mH = GammaH_*mH_;
  if ( includeHtoTauTauDecay_ ) {
    if ( !br_isInitialized_ ) {
      std::cerr << "Error in <me_ggH_nlo_lit::getMatrixElement_HtoTauTauDecay>: Branching ratio for decay of Higgs into taus has not been initialized !!" << std::endl;
      assert(0);
    }
    double meHtoTauTau_q2 = 2.*(tauLeptonMass2/v2)*mH2_*(1. - 4.*tauLeptonMass2/mH2_);    
    me *= meHtoTauTau_q2;    
    if ( !applyNWA_ ) {
      me *= (one_over_Pi*GammaH_times_mH/(square(q2_ - mH2_) + square(GammaH_times_mH)));
    }
    double meHtoTauTau_mH = 2.*(tauLeptonMass2/v2)*mH2_*(1. - 4.*tauLeptonMass2/mH2_);   
    double GammaHtoTauTau = meHtoTauTau_mH*TMath::Sqrt(1. - 4.*tauLeptonMass2/q2_)/(16.*TMath::Pi()*mH_);
    double GammaH_times_mH_fromBR = (GammaHtoTauTau/br_)*mH_;
    me /= (one_over_Pi*GammaH_times_mH_fromBR);
  } else {
    if ( !applyNWA_ ) {
      me *= (one_over_Pi*GammaH_times_mH/(square(q2_ - mH2_) + square(GammaH_times_mH)));
    }
  }
  return me;
}
