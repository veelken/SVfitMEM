#include "TauAnalysis/SVfitMEM/interface/SVfitMEM.h"

#include <TGraphErrors.h>
#include <TH1.h>

#include <algorithm>

using namespace svFitMEM;

namespace 
{
  double g(double* x, size_t dim, void* param)
  {    
    return SVfitIntegrand::gSVfitIntegrand->Eval(x);
  }
}

SVfitMEM::SVfitMEM(const std::string& madgraphFileName, double sqrtS, const std::string& pdfFileName, int verbosity) 
  : integrand_(0),
    sqrtS_(sqrtS),
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
  integrand_ = new SVfitIntegrand(madgraphFileName, sqrtS_, pdfFileName, verbosity_);

  clock_ = new TBenchmark();
}

SVfitMEM::~SVfitMEM() 
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
SVfitMEM::shiftVisPt(bool value, TFile* inputFile)
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

namespace
{
  struct sortMeasuredTauLeptons 
  {
    bool operator() (const MeasuredTauLepton& measuredTauLepton1, const MeasuredTauLepton& measuredTauLepton2)
    {
      if ( (measuredTauLepton1.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton1.type() == MeasuredTauLepton::kTauToMuDecay) && 
	   measuredTauLepton2.type() == MeasuredTauLepton::kTauToHadDecay  ) return true;
      if ( (measuredTauLepton2.type() == MeasuredTauLepton::kTauToElecDecay || measuredTauLepton2.type() == MeasuredTauLepton::kTauToMuDecay) && 
	   measuredTauLepton1.type() == MeasuredTauLepton::kTauToHadDecay ) return false;
      return ( measuredTauLepton1.pt() > measuredTauLepton2.pt() );
    }
  };
}

void
SVfitMEM::integrate(const std::vector<MeasuredTauLepton>& measuredTauLeptons, double measuredMETx, double measuredMETy, const TMatrixD& covMET, const std::string& likelihoodFileName)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<SVfitMEM::integrate>:" << std::endl;
    clock_->Start("<SVfitMEM::integrate>");
  }
  
  std::vector<MeasuredTauLepton> measuredTauLeptons_rounded;
  for ( std::vector<MeasuredTauLepton>::const_iterator measuredTauLepton = measuredTauLeptons.begin();
	measuredTauLepton != measuredTauLeptons.end(); ++measuredTauLepton ) {
    MeasuredTauLepton measuredTauLepton_rounded(
      measuredTauLepton->type(), 
      roundToNdigits(measuredTauLepton->pt()), 
      roundToNdigits(measuredTauLepton->eta()), 
      roundToNdigits(measuredTauLepton->phi()), 
      roundToNdigits(measuredTauLepton->mass()), 
      measuredTauLepton->decayMode());
    measuredTauLeptons_rounded.push_back(measuredTauLepton_rounded);
  }
  // for the VEGAS integration the order of MeasuredTauLeptons matters,
  // due to the choice of order in which the integration boundaries are defined.
  // The leptonic tau decay should go before the hadronic tau decays.
  // In case both taus decay to leptons or both taus decay hadronically,
  // the higher pT tau should go before the lower pT tau.
  std::sort(measuredTauLeptons_rounded.begin(), measuredTauLeptons_rounded.end(), sortMeasuredTauLeptons());
  measuredTauLeptons_ = measuredTauLeptons_rounded;
  if ( verbosity_ >= 1 ) {
    for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
      const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
      std::cout << "measuredTauLepton #" << idx << " (type = " << measuredTauLepton.type() << "): Pt = " << measuredTauLepton.pt() << "," 
		<< " eta = " << measuredTauLepton.eta() << " (theta = " << measuredTauLepton.p3().theta() << ")" << ", phi = " << measuredTauLepton.phi() << "," 
		<< " mass = " << measuredTauLepton.mass() << std::endl;
    }
  }
  double measuredMETx_rounded = roundToNdigits(measuredMETx);
  double measuredMETy_rounded = roundToNdigits(measuredMETy);
  TMatrixD covMET_rounded(2,2);
  covMET_rounded[0][0] = roundToNdigits(covMET[0][0]);
  covMET_rounded[1][0] = roundToNdigits(covMET[1][0]);
  covMET_rounded[0][1] = roundToNdigits(covMET[0][1]);
  covMET_rounded[1][1] = roundToNdigits(covMET[1][1]);

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

  numDimensions_ = 0; 
  
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
    if ( idx == 0 ) {      
      idxLeg1_t = numDimensions_;
      numDimensions_ += 1;
      idxLeg1_phi = numDimensions_;
      numDimensions_ += 1;
      if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay ) { 
	if ( shiftVisPt_ ) {
	  if ( measuredTauLepton.decayMode() == 0 ) {
	    leg1lutVisPtRes = lutVisPtResDM0_;
	  } else if ( measuredTauLepton.decayMode() == 1 || measuredTauLepton.decayMode() == 2 ) {
	    leg1lutVisPtRes = lutVisPtResDM1_;
	  } else if ( measuredTauLepton.decayMode() == 10 ) {
	    leg1lutVisPtRes = lutVisPtResDM10_;
	  } else {
	    std::cerr << "Warning: shiftVisPt is enabled, but leg1 decay mode = " << measuredTauLepton.decayMode() << " is not supported" 
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
    }
    if ( idx == 1 ) {
      idxLeg2_X = numDimensions_;
      numDimensions_ += 1;
      idxLeg2_phi = numDimensions_;
      numDimensions_ += 1;
      if ( measuredTauLepton.type() == MeasuredTauLepton::kTauToHadDecay ) { 
	if ( shiftVisPt_ ) {
	  if ( measuredTauLepton.decayMode() == 0 ) {
	    leg2lutVisPtRes = lutVisPtResDM0_;
	  } else if ( measuredTauLepton.decayMode() == 1 || measuredTauLepton.decayMode() == 2 ) {
	    leg2lutVisPtRes = lutVisPtResDM1_;
	  } else if ( measuredTauLepton.decayMode() == 10 ) {
	    leg2lutVisPtRes = lutVisPtResDM10_;
	  } else {
	    std::cerr << "Warning: shiftVisPt is enabled, but leg2 decay mode = " << measuredTauLepton.decayMode() << " is not supported" 
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
    }
  }

  integrand_->setInputs(measuredTauLeptons_rounded, measuredMETx_rounded, measuredMETy_rounded, covMET_rounded);
  integrand_->shiftVisPt(shiftVisPt_, leg1lutVisPtRes, leg2lutVisPtRes);
  integrand_->setIdxLeg1_t(idxLeg1_t);
  integrand_->setIdxLeg1_phi(idxLeg1_phi);
  integrand_->setIdxLeg1VisPtShift(idxLeg1VisPtShift);
  integrand_->setIdxLeg1_mNuNu(idxLeg1_mNuNu);
  integrand_->setIdxLeg2_X(idxLeg2_X);
  integrand_->setIdxLeg2_phi(idxLeg2_phi);
  integrand_->setIdxLeg2VisPtShift(idxLeg2VisPtShift);
  integrand_->setIdxLeg2_mNuNu(idxLeg2_mNuNu);
  SVfitIntegrand::gSVfitIntegrand = integrand_;

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

  assert(measuredTauLeptons_.size() == 2);
  double mVis = (measuredTauLeptons_[0].p4() + measuredTauLeptons_[1].p4()).mass();
  double mTest = 1.0125*mVis;

  std::vector<double> xGraph;
  std::vector<double> xErrGraph;
  std::vector<double> yGraph;
  std::vector<double> yErrGraph;

  double pMax = 0.;
  unsigned numMassParBelowThreshold = 0;
  bool skipHighMassTail = false;

//--- call VEGAS routine (part of GNU scientific library)
//    to perform actual integration
  double p    = 0.; 
  double pErr = 0.;
  while ( mTest < sqrtS_ && !skipHighMassTail ) {
    integrand_->setMtest(mTest);
    
    p = 0.;
    pErr = 0.;

    gsl_monte_vegas_init(vegasWorkspace_);
    vegasWorkspace_->stage = 0;
    gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsGridOpt_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &p, &pErr);
    vegasWorkspace_->stage = 1;

    // CV: repeat integration in case chi2 of estimated integral/uncertainty values
    //     indicates that result of integration cannot be trusted
    //    (up to maxIntEvalIter times in total)
    unsigned iteration = 0;
    double chi2 = -1.;
    do {
      gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_, numCallsIntEval_/vegasWorkspace_->iterations, vegasRnd_, vegasWorkspace_, &p, &pErr);
      vegasWorkspace_->stage = 3;
      ++iteration;
      chi2 = vegasWorkspace_->chisq;
    } while ( chi2 > maxChi2_ && iteration < maxIntEvalIter_ );	
    
    if ( verbosity_ >= 1 ) {
      std::cout << " M(test) = " << mTest << ": p = " << p << " +/- " << pErr << " (chi2 = " << chi2 << ")" << std::endl;
    }
 
    // CV: in order to reduce computing time, skip precise computation of integral
    //     if in high mass tail and probability negligible anyway
    if ( p > pMax ) pMax = p;
    if ( pMax > 1.e-20 && (p + 3.*TMath::Abs(pErr)) < (pMax*precision_) ) ++numMassParBelowThreshold;
    else numMassParBelowThreshold = 0;
    if ( numMassParBelowThreshold >= 5 ) {
      skipHighMassTail = true;
    }

    xGraph.push_back(mTest);
    xErrGraph.push_back(0.0125*mTest);
    yGraph.push_back(p);
    yErrGraph.push_back(pErr);

    mTest *= 1.025;
  }

  TGraphErrors* likelihoodGraph = makeGraph("svFitLikelihoodGraph", xGraph, xErrGraph, yGraph, yErrGraph);
  extractResult(likelihoodGraph, mass_, massErr_, Lmax_, verbosity_);
  if ( verbosity_ >= 1 ) {
    std::cout << "--> M = " << mass_ << " +/- " << massErr_ << " (Lmax = " << Lmax_ << ")" << std::endl;
  }
  
  if ( likelihoodFileName != "" ) {
    TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
    likelihoodGraph->Write();
    delete likelihoodFile;
  }

  delete likelihoodGraph;

  delete [] xl_;
  delete [] xu_;

  delete [] (double*)vegasIntegrand_->params;
  delete vegasIntegrand_;
  gsl_monte_vegas_free(vegasWorkspace_);
  gsl_rng_free(vegasRnd_);

  if ( verbosity_ >= 1 ) {
    clock_->Show("<SVfitMEM2::integrate>");
  }
}
