
/**
   \class testSVfitMEM testSVfitMEM.cc "TauAnalysis/SVfitStandalone/bin/testSVfitMEM.cc"
   \brief Basic example of the use of the standalone version of SVfit based on matrix elements
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/SVfitMEM.h"
#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TString.h"

using namespace svFitMEM;

namespace
{
  std::string findFile(const std::string& fileName)
  {
    edm::FileInPath inputFile(fileName);
    if ( !inputFile.isLocal() ) {
      std::cerr << "Error: Cannot find file = " << fileName << " !!" << std::endl;
      assert(0);
    }
    return inputFile.fullPath().data();
  }

  const TGraphErrors* readGraphErrors(TFile* inputFile, const std::string& graphName)
  {
    TGraphErrors* graph = dynamic_cast<TGraphErrors*>(inputFile->Get(graphName.data()));
    if ( !graph ) {
      std::cerr << "<readGraph>: Failed to load graph = " << graphName << " from file = " << inputFile->GetName() << " !!" << std::endl;
      assert(0);
    }
    return (TGraphErrors*)graph->Clone();
  }
}

void singleEvent()
{
  /* 
     This is a single event for testing purposes.
  */

  // define MET
  double measuredMETx =  18.24 - 20.00*TMath::Cos(-0.906);
  double measuredMETy = -23.07 - 20.00*TMath::Sin(-0.906); 

  // define MET covariance
  TMatrixD covMET(2,2);
  covMET[0][0] = 100.00;
  covMET[1][0] =   0.00;
  covMET[0][1] =   0.00;
  covMET[1][1] = 100.00;

  // define visible tau decay products
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  //measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay,  22.76, -0.566, -0.906, muonMass)); // tau -> muon decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, 42.76, -0.566, -0.906, chargedPionMass, 0)); // tau -> hadron decay (Pt, eta, phi, mass, tauDecayMode)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, 55.05, -1.543,  2.749, 0.712, 1)); // tau -> hadron decay (Pt, eta, phi, mass, tauDecayMode)
  /*
     tauDecayModes:  0 one-prong without neutral pions
                     1 one-prong with neutral pions
		    10 three-prong without neutral pions
  */

  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 13.e+3; 
  
  //int mode = SVfitIntegrand::kMadgraph;
  int mode = SVfitIntegrand::kLiterature;

  // CV: remove ".gz" suffix, as it is internally added by LHAPDF 
  std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/cteq65.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();

  std::string madgraphFileName = "TauAnalysis/SVfitMEM/data/param_card.dat";

  int verbosity = 1;
  SVfitMEM svFitAlgo(sqrtS, pdfFileName.data(), mode, findFile(madgraphFileName), verbosity);
  //-----------------------------------------------------------------------------
  // CV: enable the following lines to take experimental resolution on hadronic tau energy into account
  //std::string visPtResFileName = findFile("TauAnalysis/SVfitMEM/data/svFitVisMassAndPtResolutionPDF.root");
  //TH1::AddDirectory(false);  
  //TFile* visPtResFile = new TFile(visPtResFileName.data());
  //svFitAlgo.shiftVisPt(true, visPtResFile);
  //delete visPtResFile;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // CV: enable the following lines to take cross-section*signal acceptance/efficiency into account
  std::string xSection_times_AccFileName = "TauAnalysis/SVfitMEM/data/testHttXsectionWithTauDecays_hadhad.root";
  TFile* xSection_times_AccFile = new TFile(findFile(xSection_times_AccFileName).data());
  const TGraphErrors* xSection_times_AccGraph = readGraphErrors(xSection_times_AccFile, "graph_Xsection_wAcc");
  svFitAlgo.setCrossSection_times_Acc(xSection_times_AccGraph);
  delete xSection_times_AccFile;
  //-----------------------------------------------------------------------------
  svFitAlgo.setMaxObjFunctionCalls(2000);
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET, "testSVfitMEM.root");
  double mass = svFitAlgo.mass();
  double massErr = svFitAlgo.massErr();
  double Lmax = svFitAlgo.Lmax();
  std::cout << "mass = " << mass << " +/- " << massErr << ", Lmax = " << Lmax << " (expected values = 127.13, 23.47, 1.92e-18)" << std::endl;

  return;
}

int main(int argc, char* argv[]) 
{
  singleEvent();
  return 0;
}
