#include "TauAnalysis/SVfitMEM/interface/HttXsectionWithTauDecays.h"

#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include <TH1.h>

#include <algorithm>

using namespace svFitMEM;

namespace 
{
  double g(double* x, size_t dim, void* param)
  {    
    return HttXsectionIntegrandWithTauDecays::gHttXsectionIntegrandWithTauDecays->Eval(x);
  }
}

HttXsectionWithTauDecays::HttXsectionWithTauDecays(const std::string& madgraphFileName, double sqrtS, double mH, const std::string& pdfFileName, int verbosity) 
  : integrand_(0),
    sqrtS_(sqrtS),
    mH_(mH),
    vegasIntegrand_(0),
    vegasWorkspace_(0),
    vegasRnd_(0),
    numCallsGridOpt_(2000),
    numCallsIntEval_(8000),
    maxChi2_(2.),
    maxIntEvalIter_(5),
    precision_(1.e-5),
    numDimensions_(0),
    xl_(0),
    xu_(0),
    shiftVisPt_(false),
    lutVisPtResDM0_(0),
    lutVisPtResDM1_(0),
    lutVisPtResDM10_(0),
    verbosity_(verbosity)
{ 
  integrand_ = new HttXsectionIntegrandWithTauDecays(madgraphFileName, sqrtS_, pdfFileName, verbosity_);
  integrand_->setMtest(mH);

  clock_ = new TBenchmark();
}

HttXsectionWithTauDecays::~HttXsectionWithTauDecays() 
{
  delete integrand_;

  delete lutVisPtResDM0_;
  delete lutVisPtResDM1_;
  delete lutVisPtResDM10_;

  delete clock_;
}

namespace
{
  TH1* readHistogram(TFile* inputFile, const std::string& histogramName)
  {
    TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
    if ( !histogram ) {
      std::cerr << "<readHistogram>: Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
      assert(0);
    }
    return (TH1*)histogram->Clone();
  }
}

void 
HttXsectionWithTauDecays::shiftVisPt(bool value, TFile* inputFile)
{
  shiftVisPt_ = value;
  if ( shiftVisPt_ ) {
    delete lutVisPtResDM0_;
    lutVisPtResDM0_ = readHistogram(inputFile, "recTauPtDivGenTauPt_recDecayModeEq0");
    delete lutVisPtResDM1_;
    lutVisPtResDM1_ = readHistogram(inputFile, "recTauPtDivGenTauPt_recDecayModeEq1");
    delete lutVisPtResDM10_;
    lutVisPtResDM10_ = readHistogram(inputFile, "recTauPtDivGenTauPt_recDecayModeEq10");
  }
}

void
HttXsectionWithTauDecays::integrate(int tau1Type, int tau1DecayMode, double tau1Mass, int tau2Type, int tau2DecayMode, double tau2Mass, const TMatrixD& covMET)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<HttXsectionWithTauDecays::integrate>:" << std::endl;
    clock_->Start("<HttXsectionWithTauDecays::integrate>");
  }
  
//--- determine dimension of integration space 
  int idxLeg1_t = -1;
  int idxLeg1_phi = -1;
  int idxLeg1VisPtShift = -1;
  int idxLeg1_mNuNu = -1;
  const TH1* leg1lutVisPtRes = 0;

  int idxLeg2_X = -1;
  int idxLeg2_phi = -1;
  int idxLeg2VisPtShift = -1;
  int idxLeg2_mNuNu = -1;
  const TH1* leg2lutVisPtRes = 0;

  numDimensions_ = 8; // P, theta, phi of visTau1 and visTau2; Px and Py components of missing-ET
  
  idxLeg1_t = numDimensions_;
  numDimensions_ += 1;
  idxLeg1_phi = numDimensions_;
  numDimensions_ += 1;
  if ( tau1Type == MeasuredTauLepton::kTauToHadDecay ) { 
    if ( shiftVisPt_ ) {
      if ( tau1DecayMode == 0 ) {
	leg1lutVisPtRes = lutVisPtResDM0_;
      } else if ( tau1DecayMode == 1 || tau1DecayMode == 2 ) {
	leg1lutVisPtRes = lutVisPtResDM1_;
      } else if ( tau1DecayMode == 10 ) {
	leg1lutVisPtRes = lutVisPtResDM10_;
      } else {
	std::cerr << "Warning: shiftVisPt is enabled, but leg1 decay mode = " << tau1DecayMode << " is not supported" 
		  << " --> disabling shiftVisPt for this event !!" << std::endl;
      }
    }
    if ( leg1lutVisPtRes ) {
      idxLeg1VisPtShift = numDimensions_;
      ++numDimensions_;
    }
  } else {
    idxLeg1_mNuNu = numDimensions_;
    numDimensions_ += 1;
  }
  idxLeg2_X = numDimensions_;
  numDimensions_ += 1;
  idxLeg2_phi = numDimensions_;
  numDimensions_ += 1;
  if ( tau2Type == MeasuredTauLepton::kTauToHadDecay ) { 
    if ( shiftVisPt_ ) {
      if ( tau2DecayMode == 0 ) {
	leg2lutVisPtRes = lutVisPtResDM0_;
      } else if ( tau2DecayMode == 1 || tau2DecayMode == 2 ) {
	leg2lutVisPtRes = lutVisPtResDM1_;
      } else if ( tau2DecayMode == 10 ) {
	leg2lutVisPtRes = lutVisPtResDM10_;
      } else {
	std::cerr << "Warning: shiftVisPt is enabled, but leg2 decay mode = " << tau2DecayMode << " is not supported" 
		  << " --> disabling shiftVisPt for this event !!" << std::endl;
      }
    }
    if ( leg2lutVisPtRes ) {
      idxLeg2VisPtShift = numDimensions_;
      ++numDimensions_;
    }
  } else {
    idxLeg2_mNuNu = numDimensions_;
    numDimensions_ += 1;
  }

  integrand_->setInputs(tau1Type, tau1Mass, tau2Type, tau2Mass, covMET);
  integrand_->shiftVisPt(shiftVisPt_, leg1lutVisPtRes, leg2lutVisPtRes);
  integrand_->setIdxLeg1_t(idxLeg1_t);
  integrand_->setIdxLeg1_phi(idxLeg1_phi);
  integrand_->setIdxLeg1VisPtShift(idxLeg1VisPtShift);
  integrand_->setIdxLeg1_mNuNu(idxLeg1_mNuNu);
  integrand_->setIdxLeg2_X(idxLeg2_X);
  integrand_->setIdxLeg2_phi(idxLeg2_phi);
  integrand_->setIdxLeg2VisPtShift(idxLeg2VisPtShift);
  integrand_->setIdxLeg2_mNuNu(idxLeg2_mNuNu);
  HttXsectionIntegrandWithTauDecays::gHttXsectionIntegrandWithTauDecays = integrand_;

  vegasIntegrand_ = new gsl_monte_function;
  vegasIntegrand_->f = &g;
  vegasIntegrand_->dim = numDimensions_;
  vegasIntegrand_->params = new double[1];
  vegasWorkspace_ = gsl_monte_vegas_alloc(numDimensions_);
  gsl_rng_env_setup();
  vegasRnd_ = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(vegasRnd_, 12345); 
  
  //std::cout << "numDimensions = " << numDimensions_ << std::endl;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
  xl_[0] = tau1Mass;                    // vis1E
  xu_[0] = 1.e+3;
  xl_[1] = 0.;                          // vis1Theta
  xu_[1] = TMath::Pi();
  xl_[2] = -TMath::Pi();                // vis1Phi
  xu_[2] = +TMath::Pi();
  xl_[3] = square(tau1Mass + tau2Mass); // qVis2
  xu_[3] = square(2.*mH_);
  xl_[4] = 0.;                          // vis2Theta
  xu_[4] = TMath::Pi();
  xl_[5] = -TMath::Pi();                // vis2Phi
  xu_[5] = +TMath::Pi();
  xl_[6] = -100.;                       // MET: recPx - genPx
  xu_[6] = +100.;        
  xl_[7] = -100.;                       // MET: recPy - genPy
  xu_[7] = +100.;  
  xl_[idxLeg1_t] = -0.5*TMath::Pi();
  xu_[idxLeg1_t] = +0.5*TMath::Pi();
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
  xl_[idxLeg2_X] = 0.;
  xu_[idxLeg2_X] = 1.;
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

  gsl_monte_vegas_init(vegasWorkspace_);
  vegasWorkspace_->stage = 0;
  xSection_ = 0.;
  xSectionErr_ = 0.;
  gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsGridOpt_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &xSection_, &xSectionErr_);
  vegasWorkspace_->stage = 1;
  xSection_ = 0.;
  xSectionErr_ = 0.;

  // CV: repeat integration in case chi2 of estimated integral/uncertainty values
  //     indicates that result of integration cannot be trusted
  //    (up to maxIntEvalIter times in total)
  unsigned iteration = 0;
  double chi2 = -1.;
  do {
    gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsIntEval_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &xSection_, &xSectionErr_);
    vegasWorkspace_->stage = 3;
    ++iteration;
    chi2 = vegasWorkspace_->chisq;
  } while ( (chi2 > maxChi2_ || xSectionErr_ > (xSection_*precision_)) && iteration < maxIntEvalIter_ );	
  
  if ( verbosity_ >= 1 ) {
    std::cout << "--> cross-section = " << xSection_ << " +/- " << xSectionErr_ << std::endl;
  }

  delete [] xl_;
  delete [] xu_;

  delete [] (double*)vegasIntegrand_->params;
  delete vegasIntegrand_;
  gsl_monte_vegas_free(vegasWorkspace_);
  gsl_rng_free(vegasRnd_);

  if ( verbosity_ >= 1 ) {
    clock_->Show("<HttXsectionWithTauDecays::integrate>");
  }
}
