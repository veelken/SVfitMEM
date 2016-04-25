
/**
   \class testSVfitMEM testSVfitMEM.cc "TauAnalysis/SVfitMEM/bin/testSVfitMEM.cc"
   \brief Basic example of the use of the standalone version of the SVfit algorithm, based on the matrix element method
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/SVfitMEM_lo.h"
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
    if ( inputFile.fullPath() == "" ) {
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

int main(int argc, char* argv[]) 
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
  
  //int mode = SVfitIntegrand_lo::kMadgraph;
  int mode = SVfitIntegrand_lo::kLiterature;

  //std::string pdfName = "cteq66";
  std::string pdfName = "MSTW2008lo68cl";

  std::string madgraphFileName = "TauAnalysis/SVfitMEM/data/param_card.dat";

  int verbosity = 1;
  SVfitMEM_lo svFitAlgo(sqrtS, pdfName.data(), mode, findFile(madgraphFileName), verbosity);
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
  std::string xSection_times_AccFileName = "TauAnalysis/SVfitMEM/data/svFitMEM_xSection_and_AccCorr_hadhad.root";
  TFile* xSection_times_AccFile = new TFile(findFile(xSection_times_AccFileName).data());
  const TGraphErrors* graph_xSection = readGraphErrors(xSection_times_AccFile, "graph_Xsection_woAcc_hadhad_vamp");
  //const TGraphErrors* graph_Acc = readGraphErrors(xSection_times_AccFile, "graph_Acc_hadhad_vamp");
  //svFitAlgo.setCrossSection_and_Acc(graph_xSection, graph_Acc, 1.e-2);
  svFitAlgo.setCrossSection(graph_xSection);
  delete xSection_times_AccFile;
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // CV: enable the following line to add an additional log(mTauTau) term to the nll to suppress high mass tail in mTauTau distribution (default is false)
  svFitAlgo.addLogM(true, 1.e+1);
  //-----------------------------------------------------------------------------
  svFitAlgo.setMaxObjFunctionCalls(20000);
  svFitAlgo.setIntMode(SVfitMEM_lo::kVAMP);
  //svFitAlgo.setIntMode(SVfitMEM_lo::kVEGAS);
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET, "testSVfitMEM.root");
  double mass = svFitAlgo.mass();
  double massErr = svFitAlgo.massErr();
  double Lmax = svFitAlgo.Lmax();
  std::cout << "mass = " << mass << " +/- " << massErr << ", Lmax = " << Lmax << " (expected values = 121.505 +/- 15.281, Lmax = 1.11915e-10)" << std::endl;

  return 0;
}
