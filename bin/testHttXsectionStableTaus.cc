
/**
   \class testHttXsectionStableTaus testHttXsectionStableTaus.cc "TauAnalysis/SVfitStandalone/bin/testHttXsectionStableTaus.cc"
   \brief Compute leading order gg -> Higgs -> tautau cross-section, assuming taus are stable particles
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/HttXsectionStableTaus.h"

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

int main(int argc, char* argv[]) 
{
  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 13.e+3; 

  // define Higgs mass points 
  std::vector<double> mH;
  mH.push_back(90.);
  //mH.push_back(100.);
  //mH.push_back(105.);
  //mH.push_back(110.);
  //mH.push_back(115.);
  //mH.push_back(120.);
  mH.push_back(125.);
  //mH.push_back(130.);
  //mH.push_back(135.);
  //mH.push_back(140.);
  //mH.push_back(145.);
  //mH.push_back(150.);
  mH.push_back(160.);
  //mH.push_back(200.);
  //mH.push_back(250.);
  //mH.push_back(300.);

  // CV: compare cross-sections computed using MadGraph 
  //     with values computed by LHC XS working-group
  //    - cross-sections for 8 TeV center-of-mass energy from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV ,
  //      scaled by parton luminosity ratios from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
  //    - branching ratios from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
  //    
  // NOTE: The cross-sections computed using MadGraph are expected to be lower by about a factor 2,
  //       as they are leading order only.
  //       The k-factors are taken from http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/HIGGS/higgs-xsec/cross.pdf .
  //
  std::map<double, double> xSection_targets;
  // format: parton luminosity ratio * NLO cross-section [pb] * branching ratio / k factor
  xSection_targets[90.]  = 2.137*36.23*8.33e-02/1.89;	
  xSection_targets[100.] = 2.185*29.68*8.28e-02/1.89;
  xSection_targets[105.] = 2.208*27.01*8.17e-02/1.89;
  xSection_targets[110.] = 2.230*24.70*7.95e-02/1.90;
  xSection_targets[115.] = 2.253*22.66*7.58e-02/1.90;
  xSection_targets[120.] = 2.274*20.86*7.04e-02/1.90;
  xSection_targets[125.] = 2.296*19.27*6.32e-02/1.90;
  xSection_targets[130.] = 2.317*17.85*5.45e-02/1.90;
  xSection_targets[135.] = 2.338*16.57*4.49e-02/1.91;
  xSection_targets[140.] = 2.358*15.42*3.52e-02/1.91;
  xSection_targets[145.] = 2.379*14.46*2.61e-02/1.91;
  xSection_targets[150.] = 2.399*13.55*1.78e-02/1.91;
  xSection_targets[160.] = 2.439*11.96*3.96e-03/1.92;
  xSection_targets[200.] = 2.592*7.081*2.87e-04/1.94;
  xSection_targets[250.] = 2.776*4.783*1.27e-04/1.94;
  xSection_targets[300.] = 2.956*3.594*7.34e-05/1.94;

  std::map<double, std::string> madgraphFileNames;
  madgraphFileNames[90.]  = "TauAnalysis/SVfitMEM/data/param_card_mH90.dat";
  madgraphFileNames[125.] = "TauAnalysis/SVfitMEM/data/param_card_mH125.dat";
  madgraphFileNames[160.] = "TauAnalysis/SVfitMEM/data/param_card_mH160.dat";
  madgraphFileNames[200.] = "TauAnalysis/SVfitMEM/data/param_card_mH200.dat";
  madgraphFileNames[250.] = "TauAnalysis/SVfitMEM/data/param_card_mH250.dat";
  madgraphFileNames[300.] = "TauAnalysis/SVfitMEM/data/param_card_mH300.dat";

  // CV: remove ".gz" suffix, as it is internally added by LHAPDF 
  std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/cteq65.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();

  int verbosity = 0;

  for ( std::vector<double>::const_iterator mH_i = mH.begin();
	mH_i != mH.end(); ++mH_i ) {
    HttXsectionStableTaus HttXsection(findFile(madgraphFileNames[*mH_i]), sqrtS, *mH_i, pdfFileName.data(), verbosity);
    HttXsection.setMaxObjFunctionCalls(500000); 
    HttXsection.integrate();
    double xSection = HttXsection.xSection();
    double xSectionErr = HttXsection.xSectionErr();
    std::cout << "mH = " << (*mH_i) << " (GammaH = " << HttXsection.GammaH() << "):" 
	      << " cross-section = " << xSection << " +/- " << xSectionErr << " pb" 
	      << " (expected = " << xSection_targets[*mH_i] << " pb, ratio = " << xSection/xSection_targets[*mH_i] << ")" << std::endl;
  }

  return 0;
}
