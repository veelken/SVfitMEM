#include "TauAnalysis/SVfitMEM/interface/SVfitMEM.h"

#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorMarkovChain.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVEGAS.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitIntegratorVAMP.h"

#include <TGraphErrors.h>
#include <TH1.h>

#include <algorithm>

using namespace svFitMEM;

namespace 
{
  double g_C(double* x, size_t dim, void* param)
  {    
    //std::cout << "<g_C>:" << std::endl;
    double retVal = SVfitIntegrand::gSVfitIntegrand->Eval(x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }

  double g_Fortran(double** x, size_t dim, void** param)
  {    
    //std::cout << "<g_Fortran>:" << std::endl;
    double retVal = SVfitIntegrand::gSVfitIntegrand->Eval(*x);
    //std::cout << " retVal = " <<  retVal << std::endl;
    return retVal;
  }
}

SVfitMEM::SVfitMEM(double sqrtS, const std::string& pdfName, int mode, const std::string& madgraphFileName, int verbosity) 
  : integrand_(0),
    sqrtS_(sqrtS),
    intMode_(kVAMP),
    intAlgo_(0),
    maxObjFunctionCalls_(20000),
    precision_(1.e-3),
    numDimensions_(0),
    xl_(0),
    xu_(0),
    graph_xSection_(0),
    graph_Acc_(0),
    minAcc_(1.e-2),
    addLogM_(true), 
    addLogM_power_(7.), // CV: best compatibility with "old" SVfitStandalone algorithm
    shiftVisPt_(false),
    lutVisPtResDM0_(0),
    lutVisPtResDM1_(0),
    lutVisPtResDM10_(0),
    clock_(0),
    verbosity_(verbosity)
{ 
  integrand_ = new SVfitIntegrand(sqrtS_, pdfName, mode, madgraphFileName, verbosity_);

  clock_ = new TBenchmark();
}

SVfitMEM::~SVfitMEM() 
{
  delete integrand_;

  // CV: do not delete graph_xSection and graph_Acc,
  //     as they will be deleted by calling code!!
  //delete graph_xSection_;
  //delete graph_Acc_;

  delete lutVisPtResDM0_;
  delete lutVisPtResDM1_;
  delete lutVisPtResDM10_;

  delete clock_;
}

void SVfitMEM::setCrossSection(const TGraphErrors* graph_xSection)
{
  graph_xSection_ = graph_xSection;
  graph_Acc_ = 0;
  minAcc_ = 1.;
}

void SVfitMEM::setCrossSection_and_Acc(const TGraphErrors* graph_xSection, const TGraphErrors* graph_Acc, double minAcc)
{
  graph_xSection_ = graph_xSection;
  graph_Acc_ = graph_Acc;
  minAcc_ = minAcc;
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

  double linearInterpolateY(double mTest, double x0, double x1, double y0, double y1)
  {
    double weight0 = (mTest - x0)/(x1 - x0);
    assert(weight0 >= 0 && weight0 <= 1);
    double weight1 = (x1 - mTest)/(x1 - x0);
    assert(weight1 >= 0 && weight1 <= 1);
    const double epsilon = 1.e-3;
    assert((weight0 + weight1) > (1 - epsilon) && (weight0 + weight1) < (1 + epsilon));
    return (weight0*y0 + weight1*y1);
  }

  double compCrossSection_or_Acc(const TGraphErrors* graph, double mTest, double& xSection_or_AccErr)
  {
    int firstPoint = 0;
    int lastPoint = graph->GetN() - 1;

    int idxPoint0 = firstPoint - 1;
    int idxPoint1 = lastPoint + 1;
    for ( int idxPoint = firstPoint; idxPoint <= lastPoint; ++idxPoint ) {
      double x, y;
      graph->GetPoint(idxPoint, x, y);
      //std::cout << "idxPoint #" << idxPoint << ": x = " << x << ", y = " << y << std::endl;
      if ( x < mTest && idxPoint > idxPoint0 ) idxPoint0 = idxPoint;
      if ( x > mTest && idxPoint < idxPoint1 ) idxPoint1 = idxPoint;
    }
    //std::cout << "idxPoint0 = " << idxPoint0 << ", idxPoint1 = " << idxPoint1 << std::endl;
    assert(idxPoint1 >= idxPoint0);
    double xSection_or_Acc = 0.;
    if ( idxPoint0 >= firstPoint && idxPoint1 <= lastPoint ) {
      double x0, y0;
      graph->GetPoint(idxPoint0, x0, y0);
      double x1, y1;
      graph->GetPoint(idxPoint1, x1, y1);
      xSection_or_Acc = linearInterpolateY(mTest, x0, x1, y0, y1);
      double yErr0 = graph->GetErrorY(idxPoint0);
      double yErr1 = graph->GetErrorY(idxPoint1);
      xSection_or_AccErr = linearInterpolateY(mTest, x0, x1, yErr0, yErr1);
    } else if ( idxPoint0 < firstPoint ) {
      double x, y;
      graph->GetPoint(firstPoint, x, y);
      xSection_or_Acc = y;
      double yErr = graph->GetErrorY(firstPoint);
      xSection_or_AccErr = yErr;
    } else if ( idxPoint1 > lastPoint ) {
      double x, y;
      graph->GetPoint(lastPoint, x, y);
      xSection_or_Acc = y;
      double yErr = graph->GetErrorY(lastPoint);
      xSection_or_AccErr = yErr;
    } else assert(0);
    return xSection_or_Acc;
  }

  struct sortGraphPoints
  {
    bool operator() (const GraphPoint& graphPoint1, const GraphPoint& graphPoint2)
    {
      return ( graphPoint1.x_ < graphPoint2.x_ );
    }
  };
}

void
SVfitMEM::integrate(const std::vector<MeasuredTauLepton>& measuredTauLeptons, double measuredMETx, double measuredMETy, const TMatrixD& covMET, const std::string& likelihoodFileName)
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<SVfitMEM::integrate>:" << std::endl;
    clock_->Reset();
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
  if ( verbosity_ >= 1 ) {
    std::cout << "MET: Px = " << measuredMETx_rounded << ", Py = " << measuredMETy_rounded << std::endl;
    std::cout << "covMET:" << std::endl;
    covMET_rounded.Print();
  }

//--- determine dimension of integration space 
  int idxLeg1_X = -1;
  int idxLeg1_phi = -1;
  int idxLeg1VisPtShift = -1;
  int idxLeg1_mNuNu = -1;
  const TH1* leg1lutVisPtRes = 0;

  int idxLeg2_t = -1;
  int idxLeg2_phi = -1;
  int idxLeg2VisPtShift = -1;
  int idxLeg2_mNuNu = -1;
  const TH1* leg2lutVisPtRes = 0;

  numDimensions_ = 0; 
  
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];
    if ( idx == 0 ) {      
      idxLeg1_X = numDimensions_;
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
      idxLeg2_t = numDimensions_;
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
  //integrand_->shiftVisPt(shiftVisPt_, leg1lutVisPtRes, leg2lutVisPtRes);
  integrand_->setIdxLeg1_X(idxLeg1_X);
  integrand_->setIdxLeg1_phi(idxLeg1_phi);
  integrand_->setIdxLeg1VisPtShift(idxLeg1VisPtShift);
  integrand_->setIdxLeg1_mNuNu(idxLeg1_mNuNu);
  integrand_->setIdxLeg2_t(idxLeg2_t);
  integrand_->setIdxLeg2_phi(idxLeg2_phi);
  integrand_->setIdxLeg2VisPtShift(idxLeg2VisPtShift);
  integrand_->setIdxLeg2_mNuNu(idxLeg2_mNuNu);
  SVfitIntegrand::gSVfitIntegrand = integrand_;

  if ( intMode_ == kMarkovChain ) {
    //unsigned numChains = TMath::Nint(maxObjFunctionCalls_/100000.);
    unsigned numChains = 1;
    unsigned numIterBurnin = TMath::Nint(0.10*maxObjFunctionCalls_/numChains);
    unsigned numIterSampling = TMath::Nint(0.90*maxObjFunctionCalls_/numChains);
    unsigned numIterSimAnnealingPhase1 = TMath::Nint(0.20*numIterBurnin);
    unsigned numIterSimAnnealingPhase2 = TMath::Nint(0.60*numIterBurnin);
    std::string treeFileName = "SVfitIntegratorMarkovChain_SVfitMEM.root";
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
    std::cerr << "<SVfitMEM::integrate>: Invalid Configuration Parameter 'intMode' = " << intMode_ << " --> ABORTING !!\n";
    assert(0);
  }

  //std::cout << "numDimensions = " << numDimensions_ << std::endl;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
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

  assert(measuredTauLeptons_.size() == 2);
  double mVis = (measuredTauLeptons_[0].p4() + measuredTauLeptons_[1].p4()).mass();

  std::vector<GraphPoint> graphPoints;
  std::set<int> mTest_computed;
  double pMax = 0.;

//--- call VEGAS routine (part of GNU scientific library)
//    to perform actual integration
  const int numIter = 3;
  for ( int idxIter = 0; idxIter < numIter; ++idxIter ) {
    double mTest_step;
    if      ( idxIter == 0 ) mTest_step = 1.1000;
    else if ( idxIter == 1 ) mTest_step = TMath::Power(1.1000, 1./4);
    else if ( idxIter == 2 ) mTest_step = TMath::Power(1.1000, 1./8);
    else assert(0);

    double mTest_min = mVis;
    double mTest_max = sqrtS_;
    if ( idxIter >= 1 ) {      
      double pMin;
      if      ( idxIter == 1 ) pMin = 1.e-2*pMax;
      else if ( idxIter == 2 ) pMin = 1.e-1*pMax;
      else assert(0);
      mTest_min = sqrtS_;
      mTest_max = 0.;
      for ( std::vector<GraphPoint>::const_iterator graphPoint = graphPoints.begin();
	    graphPoint != graphPoints.end(); ++graphPoint ) {
	if ( graphPoint->y_ > pMin ) {
	  if ( graphPoint->x_ < mTest_min ) mTest_min = graphPoint->x_;
	  if ( graphPoint->x_ > mTest_max ) mTest_max = graphPoint->x_;
	}
      }
    }
    if ( verbosity_ >= 1 ) {
      std::cout << "starting iteration #" << idxIter << ": mTest_min = " << mTest_min << ", mTest_max = " << mTest_max << std::endl;
    }

    unsigned numMassParBelowThreshold = 0;
    bool skipHighMassTail = false;

    double mTest = 1.0125*mVis;
    if ( mTest < mTest_min ) mTest = mTest_min;
    while ( mTest < mTest_max && !skipHighMassTail ) {
      if ( mTest_computed.find(TMath::Nint(mTest*1.e+2)) == mTest_computed.end() ) {
	integrand_->setMtest(mTest);
	
	double p = 0.; 
	double pErr = 0.;  
	if ( intMode_ == kMarkovChain || intMode_ == kVEGAS ) { 
	  intAlgo_->integrate(&g_C, xl_, xu_, numDimensions_, p, pErr);
	} else if ( intMode_ == kVAMP ) {
	  intAlgo_->integrate(&g_Fortran, xl_, xu_, numDimensions_, p, pErr);
	} else assert(0);    
	
	if ( graph_xSection_ ) {
	  double xSectionErr = 0.;
	  double xSection = compCrossSection_or_Acc(graph_xSection_, mTest, xSectionErr);
	  double xSection_times_Acc = xSection;
	  double xSection_times_AccErr = xSectionErr;
	  if ( graph_Acc_ ) {
	    double AccErr = 0.;
	    double Acc = compCrossSection_or_Acc(graph_Acc_, mTest, AccErr);
	    if ( Acc < minAcc_ ) Acc = minAcc_;
	    double xSection_times_AccRelErr2 = 0.;
	    if ( xSection > 0. ) xSection_times_AccRelErr2 += square(xSectionErr/xSection);
	    if ( Acc      > 0. ) xSection_times_AccRelErr2 += square(AccErr/Acc);
	    xSection_times_Acc = xSection*Acc;
	    xSection_times_AccErr = xSection_times_Acc*TMath::Sqrt(xSection_times_AccRelErr2);
	  }
	  if ( xSection_times_Acc > 0. ) {
	    p /= xSection_times_Acc;
	    if ( p > 0. && xSection_times_Acc > 0. ) pErr = p*TMath::Sqrt(square(pErr/p) + square(xSection_times_AccErr/xSection_times_Acc));
	    else pErr = 0.;
	  } else {
	    p = 0.;
	    pErr = 0.;
	  }
	}
	
	if ( addLogM_ ) {
	  double factor = TMath::Power(mTest/mVis, addLogM_power_);
	  p /= factor;
	  pErr /= factor;
	}
	
	if ( verbosity_ >= 1 ) {
	  std::cout << " M(test) = " << mTest << ": p = " << p << " +/- " << pErr << std::endl;
	}
      
	// CV: in order to reduce computing time, skip precise computation of integral
	//     if in high mass tail and probability negligible anyway
	if ( p > pMax ) pMax = p;
	if ( idxIter == 0 ) {
	  if ( pMax > 1.e-20 && (p + 3.*TMath::Abs(pErr)) < (pMax*precision_) ) ++numMassParBelowThreshold;
	  else numMassParBelowThreshold = 0;
	  if ( numMassParBelowThreshold >= 5 ) {
	    skipHighMassTail = true;
	  }
	} 
	
	GraphPoint graphPoint;
	graphPoint.x_ = mTest;
	graphPoint.xErr_ = 0.5*(mTest_step - 1.)*mTest;
	graphPoint.y_ = p;
	graphPoint.yErr_ = pErr;
	graphPoints.push_back(graphPoint);

	mTest_computed.insert(TMath::Nint(mTest*1.e+2));
      }
      
      mTest *= mTest_step;
    }    
  }

  std::sort(graphPoints.begin(), graphPoints.end(), sortGraphPoints());

  TGraphErrors* likelihoodGraph = makeGraph("svFitLikelihoodGraph", graphPoints);
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

  delete intAlgo_;

  if ( verbosity_ >= 1 ) {
    clock_->Show("<SVfitMEM::integrate>");
  }
}
