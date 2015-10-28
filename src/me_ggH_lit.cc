#include "TauAnalysis/SVfitMEM/interface/me_ggH_lit.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

const double mtop = 172.5;
const double mtop2 = svFitMEM::square(mtop);

using namespace svFitMEM;

me_ggH_lit::me_ggH_lit(const std::string& pdfName, bool applyNWA, bool includeHtoTauTauDecay)
  : applyNWA_(applyNWA),
    includeHtoTauTauDecay_(includeHtoTauTauDecay),
    s_(0.),
    s_isInitialized_(false),
    mH_(0.),
    mH2_(0.),
    mH_isInitialized_(false),
    GammaH_(0.),
    GammaH2_(0.),
    GammaH_isInitialized_(false),
    br_(0.),
    br_isInitialized_(false),
    pdf_(0),
    pdfIsInitialized_(false)
{
  // initialize PDF set
  if ( !pdfIsInitialized_ ) {
    pdf_ = LHAPDF::mkPDF(pdfName.data(), 0);
    pdfIsInitialized_ = true;
  }
}
 
me_ggH_lit::~me_ggH_lit() 
{
  delete pdf_;
}

double me_ggH_lit::getMatrixElement() const
{
  //std::cout << "<me_ggH_lit::getMatrixElement>:" << std::endl;
  double me = getMatrixElement_woHtoTauTauDecay();
  const double one_over_Pi = 1./TMath::Pi();
  double GammaH_times_mH = GammaH_*mH_;
  if ( includeHtoTauTauDecay_ ) {
    if ( !br_isInitialized_ ) {
      std::cerr << "Error in <me_ggH_lit::getMatrixElement_woBWandBR>: Branching ratio for decay of Higgs into taus has not been initialized !!" << std::endl;
      assert(0);
    }
    //double meHtoTauTau_q2 = 2.*(tauLeptonMass2/v2)*q2_*(1. - 4.*tauLeptonMass2/q2_);    
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

double me_ggH_lit::getMatrixElement_woHtoTauTauDecay() const
{
  if ( !(mH_isInitialized_ && GammaH_isInitialized_) ) {
    std::cerr << "Error in <me_ggH_lit::getMatrixElement_woBWandBR>: Higgs mass and width have not been initialized !!" << std::endl;
    assert(0);
  }

  assert(momenta_.size() == 4);
  const double* madgraphTau1P4 = momenta_[2];
  LorentzVector tau1P4(madgraphTau1P4[1], madgraphTau1P4[2], madgraphTau1P4[3], madgraphTau1P4[0]);
  const double* madgraphTau2P4 = momenta_[3];
  LorentzVector tau2P4(madgraphTau2P4[1], madgraphTau2P4[2], madgraphTau2P4[3], madgraphTau2P4[0]);
  q_ = (tau1P4 + tau2P4).mass();
  q2_ = square(q_);

  double tau = 4.*mtop2/q2_;
  
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
  //std::cout << "Re_f = " << Re_f << ", Im_f = " << Im_f << std::endl;

  const double GF = 1.166e-5; // in units of GeV^-2, taken from http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf
  const double constFactor = TMath::Sqrt(2)*GF/(256.*square(TMath::Pi()));
  double Q2 = q2_;
  assert(pdfIsInitialized_);
  double alphaS = pdf_->alphasQ2(Q2);
  //std::cout << "Q = " << Q << ": alphaS = " << alphaS << std::endl;
  double me = constFactor*square(alphaS)*square(tau)*square(q2_)*(square(1. + (1. - tau)*Re_f) + square((1. - tau)*Im_f));
  //std::cout << "me = " << me << std::endl;
  return me;
}
