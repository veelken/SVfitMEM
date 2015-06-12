
/**
   \class testSVfitMEM testSVfitMEM.cc "TauAnalysis/SVfitStandalone/bin/testSVfitMEM.cc"
   \brief Basic example of the use of the standalone version of SVfit based on matrix elements
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/SVfitMEM.h"
#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include "TFile.h"
#include "TH1.h"
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
}

void singleEvent()
{
  /* 
     This is a single event for testing purposes.
  */

  // define MET
  double measuredMETx =  18.24;
  double measuredMETy = -23.07; 

  // define MET covariance
  TMatrixD covMET(2,2);
  covMET[0][0] = 100.00;
  covMET[1][0] =   0.00;
  covMET[0][1] =   0.00;
  covMET[1][1] = 100.00;

  // define visible tau decay products
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay,  22.76, -0.566, -0.906, muonMass)); // tau -> electron decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, 55.05, -1.543,  2.749, 0.712, 1)); // tau -> hadron decay (Pt, eta, phi, mass, tauDecayMode)
  /*
     tauDecayModes:  0 one-prong without neutral pions
                     1 one-prong with neutral pions
		    10 three-prong without neutral pions
  */

  std::string madgraphFileName = findFile("TauAnalysis/SVfitMEM/data/param_card.dat");
  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 13.e+3; 
  // CV: remove ".gz" suffix, as it is internally added by LHAPDF 
  std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/cteq65.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();
  int verbosity = 1;
  SVfitMEM svFitAlgo(madgraphFileName, sqrtS, pdfFileName.data(), verbosity);
  //-----------------------------------------------------------------------------
  // CV: enable the following lines to take experimental resolution on hadronic tau energy into account
  std::string visPtResFileName = findFile("TauAnalysis/SVfitMEM/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  TFile* visPtResFile = new TFile(visPtResFileName.data());
  svFitAlgo.shiftVisPt(true, visPtResFile);
  delete visPtResFile;
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
