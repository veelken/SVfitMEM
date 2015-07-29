
/**
   \class testHttXsectionStableHiggs testHttXsectionStableHiggs.cc "TauAnalysis/SVfitStandalone/bin/testHttXsectionStableHiggs.cc"
   \brief Compute leading order gg -> Higgs -> tautau cross-section, assuming taus are stable particles
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/HttXsectionStableHiggs.h"
#include "TauAnalysis/SVfitMEM/interface/HttXsectionIntegrandStableHiggs.h"

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
    //std::cout << "<findFile>: returning '" << inputFile.fullPath() << "'" << std::endl; 
    return inputFile.fullPath().data();
  }
}

int main(int argc, char* argv[]) 
{
  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 14.e+3; 

  bool applyNWA = true;
  //bool applyNWA = false;

  //int mode = HttXsectionStableHiggs::kMadgraph;
  int mode = HttXsectionIntegrandStableHiggs::kLiterature;

  // define Higgs mass points 
  std::vector<double> mH;
  //mH.push_back(90.);
  mH.push_back(100.);
  mH.push_back(125.);
  mH.push_back(150.);
  //mH.push_back(160.);
  mH.push_back(200.);
  mH.push_back(250.);
  mH.push_back(300.);
  mH.push_back(350.);
  mH.push_back(400.);
  mH.push_back(450.);
  mH.push_back(500.);

  // CV: compare cross-sections computed using MadGraph 
  //     with values computed by LHC XS working-group
  //    - cross-sections for 8 TeV center-of-mass energy from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV ,
  //      scaled by parton luminosity ratios from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
  //    - branching ratios from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
  //    
  // NOTE: The cross-sections computed using MadGraph are expected to be lower by about a factor 2,
  //       as they are leading order only.
  //       The k-factors are computed by taking the ratio of the cross-sections computed by the LHC XS working-group
  //       and the leading order cross-sections given by Eq. 38 (with Eqs. 30, 31 and 36) of http://www.itp.phys.ethz.ch/education/fs10/aft/Thesis_MB.pdf 
  //
  std::map<double, double> branchingRatios;
  branchingRatios[90.]  = 8.33e-02;
  branchingRatios[100.] = 8.28e-02;
  branchingRatios[125.] = 6.32e-02;
  branchingRatios[150.] = 1.78e-02;
  branchingRatios[160.] = 3.96e-03;
  branchingRatios[200.] = 2.87e-04;
  branchingRatios[250.] = 1.27e-04;
  branchingRatios[300.] = 7.34e-05;
  branchingRatios[350.] = 4.76e-05;
  branchingRatios[400.] = 2.84e-05;
  branchingRatios[450.] = 1.99e-04;
  branchingRatios[500.] = 1.53e-05;

  std::map<double, double> xSection_targets_13TeV;
  // format: parton luminosity ratio * NLO cross-section [pb] / k factor
  xSection_targets_13TeV[90.]  = 2.137*36.23/1.89;	
  xSection_targets_13TeV[100.] = 2.185*29.68/1.89;
  xSection_targets_13TeV[125.] = 2.296*19.27/1.90;
  xSection_targets_13TeV[150.] = 2.399*13.55/1.91;
  xSection_targets_13TeV[160.] = 2.439*11.96/1.92;
  xSection_targets_13TeV[200.] = 2.592*7.081/1.94;
  xSection_targets_13TeV[250.] = 2.776*4.783/1.96;
  xSection_targets_13TeV[300.] = 2.956*3.594/1.99;

  std::map<double, double> xSection_targets_14TeV;
  // format: parton luminosity ratio * NLO cross-section [pb] / k factor
  xSection_targets_14TeV[90.]  = 2.385*36.23/1.89;	
  xSection_targets_14TeV[100.] = 2.446*29.68/(2.446*29.68/25.759);
  xSection_targets_14TeV[125.] = 2.587*19.27/(2.587*19.27/17.317);
  xSection_targets_14TeV[150.] = 2.720*13.55/(2.720*13.55/12.448);
  xSection_targets_14TeV[160.] = 2.771*11.96/1.92;
  xSection_targets_14TeV[200.] = 2.969*7.081/(2.969*7.081/7.366);
  xSection_targets_14TeV[250.] = 3.208*4.783/(3.208*4.783/4.993);
  xSection_targets_14TeV[300.] = 3.444*3.594/(3.444*3.594/3.829);
  xSection_targets_14TeV[350.] = 3.208*4.783/(3.208*4.783/3.532);
  xSection_targets_14TeV[400.] = 3.444*3.594/(3.444*3.594/3.786);
  xSection_targets_14TeV[450.] = 3.208*4.783/(3.208*4.783/2.846);
  xSection_targets_14TeV[500.] = 3.444*3.594/(3.444*3.594/1.960);

  std::map<double, std::string> madgraphFileNames;
  madgraphFileNames[90.]  = "TauAnalysis/SVfitMEM/data/param_card_mH90.dat";
  madgraphFileNames[100.] = "TauAnalysis/SVfitMEM/data/param_card_mH100.dat";
  madgraphFileNames[125.] = "TauAnalysis/SVfitMEM/data/param_card_mH125.dat";
  madgraphFileNames[150.] = "TauAnalysis/SVfitMEM/data/param_card_mH150.dat";
  madgraphFileNames[160.] = "TauAnalysis/SVfitMEM/data/param_card_mH160.dat";
  madgraphFileNames[200.] = "TauAnalysis/SVfitMEM/data/param_card_mH200.dat";
  madgraphFileNames[250.] = "TauAnalysis/SVfitMEM/data/param_card_mH250.dat";
  madgraphFileNames[300.] = "TauAnalysis/SVfitMEM/data/param_card_mH300.dat";
  madgraphFileNames[350.] = "TauAnalysis/SVfitMEM/data/param_card_mH350.dat";
  madgraphFileNames[400.] = "TauAnalysis/SVfitMEM/data/param_card_mH400.dat";
  madgraphFileNames[450.] = "TauAnalysis/SVfitMEM/data/param_card_mH450.dat";
  madgraphFileNames[500.] = "TauAnalysis/SVfitMEM/data/param_card_mH500.dat";

  // CV: remove ".gz" suffix, as it is internally added by LHAPDF 
  std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/cteq65.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();
  //std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/MSTW2008lo68cl.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();

  int verbosity = 0;

  for ( std::vector<double>::const_iterator mH_i = mH.begin();
	mH_i != mH.end(); ++mH_i ) {
    HttXsectionStableHiggs* HttXsection = 0;
    if ( mode == HttXsectionIntegrandStableHiggs::kMadgraph ) {
      HttXsection = new HttXsectionStableHiggs(sqrtS, *mH_i, pdfFileName.data(), applyNWA, HttXsectionIntegrandStableHiggs::kMadgraph, findFile(madgraphFileNames[*mH_i]), verbosity);
    } else if ( mode == HttXsectionIntegrandStableHiggs::kLiterature ) { 
      HttXsection = new HttXsectionStableHiggs(sqrtS, *mH_i, pdfFileName.data(), applyNWA, HttXsectionIntegrandStableHiggs::kLiterature, "", verbosity);
    } else {
      assert(0);
    }
    HttXsection->setBR(branchingRatios[*mH_i]);
    //HttXsection->setMaxObjFunctionCalls(500000); 
    HttXsection->setMaxObjFunctionCalls(10000); 
    HttXsection->integrate();
    double xSection = HttXsection->xSection();
    double xSectionErr = HttXsection->xSectionErr();
    std::cout << "mH = " << (*mH_i) << " (GammaH = " << HttXsection->GammaH() << "):" 
	      << " cross-section = " << xSection << " +/- " << xSectionErr << " pb" 
	      << " (expected = " << xSection_targets_14TeV[*mH_i] << " pb, ratio = " << xSection/xSection_targets_14TeV[*mH_i] << ")" << std::endl;
    delete HttXsection;
  }

  return 0;
}
