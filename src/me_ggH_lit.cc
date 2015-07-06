#include "TauAnalysis/SVfitMEM/interface/me_ggH_lit.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

const double mtop = 172.5;
const double mtop2 = svFitMEM::square(mtop);

using namespace svFitMEM;

namespace LHAPDF {
  double alphasPDF(int nset, double Q);
}

me_ggH_lit::me_ggH_lit(bool applyNWA)
  : applyNWA_(applyNWA),
    s_(0.),
    s_isInitialized_(false),
    mH_(0.),
    mH2_(0.),
    mH_isInitialized_(false),
    GammaH_(0.),
    GammaH2_(0.),
    GammaH_isInitialized_(false),
    br_(0.),
    br_isInitialized_(false)
{}
 
me_ggH_lit::~me_ggH_lit() 
{}

double me_ggH_lit::getMatrixElement() const
{
  //std::cout << "<me_ggH_lit::getMatrixElement>:" << std::endl;
  double me = getMatrixElement_woBWandBR();
  const double one_over_Pi = 1./TMath::Pi();
  double GammaH_div_mH = GammaH_/mH_;
  if ( !applyNWA_ ) {
    //me *= (one_over_Pi*sHat_*GammaH_div_mH/(square(sHat_ - mH2_) + square(sHat_*GammaH_div_mH)));
    me *= (one_over_Pi*mH2_*GammaH_div_mH/(square(sHat_ - mH2_) + square(mH2_*GammaH_div_mH)));
  }
  //std::cout << "br = " << br_ << std::endl;
  me *= br_;
  return me;
}

double me_ggH_lit::getMatrixElement_woBWandBR() const
{
  if ( !(mH_isInitialized_ && GammaH_isInitialized_ && br_isInitialized_) ) {
    std::cerr << "Error in <me_ggH_lit::getMatrixElement_woBWandBR>: Higgs mass, width and branching fraction have not been initialized !!" << std::endl;
    assert(0);
  }

  assert(momenta_.size() == 4);
  const double* madgraphTau1P4 = momenta_[2];
  LorentzVector tau1P4(madgraphTau1P4[1], madgraphTau1P4[2], madgraphTau1P4[3], madgraphTau1P4[0]);
  const double* madgraphTau2P4 = momenta_[3];
  LorentzVector tau2P4(madgraphTau2P4[1], madgraphTau2P4[2], madgraphTau2P4[3], madgraphTau2P4[0]);
  sHat_ = square((tau1P4 + tau2P4).mass());

  //double tau = 4.*mtop2/sHat_;
  double tau = 4.*mtop2/mH2_;
  //std::cout << "mH2 = " << mH2_ << ", sHat = " << sHat_ << std::endl;
  
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
  double Q = mH_;
  double alphaS = LHAPDF::alphasPDF(1, Q);
  //std::cout << "Q = " << Q << ": alphaS = " << alphaS << std::endl;
  double me = constFactor*square(alphaS)*square(tau)*square(sHat_)*(square(1. + (1. - tau)*Re_f) + square((1. - tau)*Im_f));
  //double f = square(TMath::ASin(1./TMath::Sqrt(tau)));
  //double me = (TMath::Sqrt(2)*GF*square(alphaS)*square(tau)*square(sHat_)/(256.*square(TMath::Pi())))*square(1. + (1. - tau)*f);
  //std::cout << "me = " << me << std::endl;
  return me;
}
