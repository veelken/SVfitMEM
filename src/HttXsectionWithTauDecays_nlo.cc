#include "TauAnalysis/SVfitMEM/interface/HttXsectionWithTauDecays_nlo.h"

#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVEGAS.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVAMP.h"

#include <TH1.h>

#include <algorithm>

using namespace svFitMEM;

namespace 
{
  double g_C(double* x, size_t dim, void* param)
  {    
    //std::cout << "<g_C>:" << std::endl;
    double retVal = HttXsectionIntegrandWithTauDecays_nlo::gHttXsectionIntegrandWithTauDecays->Eval(x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }

  double g_Fortran(double** x, size_t dim, void** param)
  {    
    //std::cout << "<g_Fortran>:" << std::endl;
    double retVal = HttXsectionIntegrandWithTauDecays_nlo::gHttXsectionIntegrandWithTauDecays->Eval(*x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }
}

HttXsectionWithTauDecays_nlo::HttXsectionWithTauDecays_nlo(double sqrtS, double mH, const std::string& pdfName, int verbosity) 
  : applyMEtTF_(false),
    applyAcceptanceCuts_(false),
    integrand_(0),
    sqrtS_(sqrtS),
    mH_(mH),
    mH2_(mH_*mH_),
    //intMode_(kMarkovChain),
    intMode_(kVEGAS),
    intAlgo_(0),
    maxObjFunctionCalls_(100000),
    numDimensions_(0),
    xl_(0),
    xu_(0),
    useHadTauTF_(false),
    clock_(0),
    verbosity_(verbosity)
{ 
  integrand_ = new HttXsectionIntegrandWithTauDecays_nlo(sqrtS_, mH_, pdfName, verbosity_);
  integrand_->setApplyMEtTF(applyMEtTF_);

  clock_ = new TBenchmark();
}

HttXsectionWithTauDecays_nlo::~HttXsectionWithTauDecays_nlo() 
{
  delete integrand_;

  delete clock_;
}

void
HttXsectionWithTauDecays_nlo::integrate(int tau1Type, int tau1DecayMode, double vis1Mass, int tau2Type, int tau2DecayMode, double vis2Mass, const TMatrixD& covMET)
{
  //if ( verbosity_ >= 1 ) {
    std::cout << "<HttXsectionWithTauDecays_nlo::integrate>:" << std::endl;
    clock_->Start("<HttXsectionWithTauDecays_nlo::integrate>");
  //}

//--- determine dimension of integration space 
  int idxLeg1_X = -1;
  int idxLeg1_phi = -1;
  int idxLeg1VisPtShift = -1;
  int idxLeg1_mNuNu = -1;

  int idxLeg2_t = -1;
  int idxLeg2_phi = -1;
  int idxLeg2VisPtShift = -1;
  int idxLeg2_mNuNu = -1;

  numDimensions_ = 7; // vis1Px, vis1Py, vis1Pz and mVis^2, hadRecoilPx, hadRecoilPy, hadRecoilPz
  
  idxLeg1_X = numDimensions_;
  numDimensions_ += 1;
  idxLeg1_phi = numDimensions_;
  numDimensions_ += 1;
  if ( tau1Type == MeasuredTauLepton::kTauToHadDecay ) { 
    if ( useHadTauTF_ ) {
      idxLeg1VisPtShift = numDimensions_;
      ++numDimensions_;
    }
  } else {
    idxLeg1_mNuNu = numDimensions_;
    numDimensions_ += 1;
  }
  idxLeg2_t = numDimensions_;
  numDimensions_ += 1;
  idxLeg2_phi = numDimensions_;
  numDimensions_ += 1;
  if ( tau2Type == MeasuredTauLepton::kTauToHadDecay ) { 
    if ( useHadTauTF_ ) {
      idxLeg2VisPtShift = numDimensions_;
      ++numDimensions_;
    }
  } else {
    idxLeg2_mNuNu = numDimensions_;
    numDimensions_ += 1;
  }

  integrand_->setInputs(tau1Type, vis1Mass, tau2Type, vis2Mass, covMET);
  if ( useHadTauTF_ ) integrand_->enableHadTauTF();
  else integrand_->disableHadTauTF();
  integrand_->setIdxLeg1_X(idxLeg1_X);
  integrand_->setIdxLeg1_phi(idxLeg1_phi);
  integrand_->setIdxLeg1VisPtShift(idxLeg1VisPtShift);
  integrand_->setIdxLeg1_mNuNu(idxLeg1_mNuNu);
  integrand_->setIdxLeg2_t(idxLeg2_t);
  integrand_->setIdxLeg2_phi(idxLeg2_phi);
  integrand_->setIdxLeg2VisPtShift(idxLeg2VisPtShift);
  integrand_->setIdxLeg2_mNuNu(idxLeg2_mNuNu);
  HttXsectionIntegrandWithTauDecays_nlo::gHttXsectionIntegrandWithTauDecays = integrand_;

  if ( intMode_ == kMarkovChain ) {
    //unsigned numChains = TMath::Nint(maxObjFunctionCalls_/100000.);
    unsigned numChains = 1;
    unsigned numIterBurnin = TMath::Nint(0.10*maxObjFunctionCalls_/numChains);
    unsigned numIterSampling = TMath::Nint(0.90*maxObjFunctionCalls_/numChains);
    unsigned numIterSimAnnealingPhase1 = TMath::Nint(0.20*numIterBurnin);
    unsigned numIterSimAnnealingPhase2 = TMath::Nint(0.60*numIterBurnin);
    std::string treeFileName = Form("SVfitIntegratorMarkovChain_HttXsectionWithTauDecays_mH%03.0f", mH_);
    if ( applyAcceptanceCuts_ ) treeFileName.append("_wAcc");
    else treeFileName.append("_woAcc");
    treeFileName.append(".root");
    intAlgo_ = new SVfitIntegratorMarkovChain(
      "uniform", 
      numIterBurnin, numIterSampling, numIterSimAnnealingPhase1, numIterSimAnnealingPhase2,
      15., 1. - 1./(0.1*numIterBurnin),
      numChains, 100, 
      1.e-2, 0.71,
      treeFileName.data());
  } else if ( intMode_ == kVEGAS ) {
    unsigned numCallsGridOpt = TMath::Nint(0.20*maxObjFunctionCalls_);
    unsigned numCallsIntEval = TMath::Nint(0.80*maxObjFunctionCalls_);
    intAlgo_ = new SVfitIntegratorVEGAS(
      numCallsGridOpt, numCallsIntEval, 
      2., 1);
  } else if ( intMode_ == kVAMP ) {
    unsigned numCallsGridOpt = TMath::Nint(0.20*maxObjFunctionCalls_);
    unsigned numCallsIntEval = TMath::Nint(0.80*maxObjFunctionCalls_);
    intAlgo_ = new SVfitIntegratorVAMP(
      numCallsGridOpt, numCallsIntEval);
  } else {
    std::cerr << "<HttXsectionWithTauDecays_nlo::integrate>: Invalid Configuration Parameter 'intMode' = " << intMode_ << " --> ABORTING !!\n";
    assert(0);
  }

  //std::cout << "numDimensions = " << numDimensions_ << std::endl;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
  xl_[0] = -1.0*0.5*mH_; // vis1Px
  xu_[0] = +1.0*0.5*mH_;
  xl_[1] = -1.0*0.5*mH_; // vis1Py
  xu_[1] = +1.0*0.5*mH_;
  xl_[2] = -TMath::Power(sqrtS_, 1./3.); // vis1Pz^1/3
  xu_[2] = +TMath::Power(sqrtS_, 1./3.);
  xl_[3] = 0.;           // mVis^2
  xu_[3] = 2.*mH2_;
  //-----------------------------------------------------------------------------
  xl_[4] = -300.;        // jetPx
  xu_[4] = +300.;        
  xl_[5] = -300.;        // jetPx
  xu_[5] = +300.;  
  xl_[6] = -sqrtS_;      // jetPz
  xu_[6] = +sqrtS_;  
  xl_[idxLeg1_X] = 0.;
  xu_[idxLeg1_X] = 1.;
  xl_[idxLeg1_phi] = -TMath::Pi();
  xu_[idxLeg1_phi] = +TMath::Pi();
  if ( idxLeg1VisPtShift != -1 ) {
    xl_[idxLeg1VisPtShift] = 0.5;
    xu_[idxLeg1VisPtShift] = 1.5;
  }
  if ( idxLeg1_mNuNu != -1 ) {
    xl_[idxLeg1_mNuNu] = 0.;
    xu_[idxLeg1_mNuNu] = tauLeptonMass2;
  }
  xl_[idxLeg2_t] = -0.5*TMath::Pi();
  xu_[idxLeg2_t] = +0.5*TMath::Pi();
  xl_[idxLeg2_phi] = -TMath::Pi();
  xu_[idxLeg2_phi] = +TMath::Pi();
  if ( idxLeg2VisPtShift != -1 ) {
    xl_[idxLeg2VisPtShift] = 0.5;
    xu_[idxLeg2VisPtShift] = 1.5;
  }
  if ( idxLeg2_mNuNu != -1 ) {
    xl_[idxLeg2_mNuNu] = 0.;
    xu_[idxLeg2_mNuNu] = tauLeptonMass2;
  }
  if ( verbosity_ >= 2 ) { 
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::cout << " fitParameter #" << iDimension << ": xl = " << xl_[iDimension] << ", xu = " << xu_[iDimension] << std::endl;
    }
  }

  if ( intMode_ == kMarkovChain || intMode_ == kVEGAS ) { 
    intAlgo_->integrate(&g_C, xl_, xu_, numDimensions_, xSection_, xSectionErr_);
  } else if ( intMode_ == kVAMP ) {
    intAlgo_->integrate(&g_Fortran, xl_, xu_, numDimensions_, xSection_, xSectionErr_);
  } else assert(0);    

  if ( verbosity_ >= 1 ) {
    std::cout << "--> cross-section = " << xSection_ << " +/- " << xSectionErr_ << std::endl;
  }

  delete [] xl_;
  delete [] xu_;

  delete intAlgo_;

  //if ( verbosity_ >= 1 ) {
    clock_->Show("<HttXsectionWithTauDecays_nlo::integrate>");
  //}
}
