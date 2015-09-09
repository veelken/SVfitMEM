
/**
   \class testHttXsectionStableTaus testHttXsectionStableTaus.cc "TauAnalysis/SVfitStandalone/bin/testHttXsectionStableTaus.cc"
   \brief Compute leading order gg -> Higgs -> tautau cross-section, assuming taus are stable particles
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/HttXsectionWithTauDecays.h"
#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"

#include <TFile.h>
#include <TGraphAsymmErrors.h>

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
}

class acceptanceTypeHadHad : public acceptanceBaseType
{
 public:
  double operator()(const svFitMEM::LorentzVector& visTau1P4, const svFitMEM::LorentzVector& visTau2P4, double metPx, double metPy) const
  {
    //std::cout << "<acceptanceHadHad>:" << std::endl;
    //std::cout << " visTau1: Pt = " << visTau1P4.pt() << ", eta = " << visTau1P4.eta() << ", phi = " << visTau1P4.phi() << ", mass = " << visTau1P4.mass() << std::endl;
    //std::cout << " visTau2: Pt = " << visTau2P4.pt() << ", eta = " << visTau2P4.eta() << ", phi = " << visTau2P4.phi() << ", mass = " << visTau2P4.mass() << std::endl;
    double acceptance = 0.;
    if ( (visTau1P4.pt() > 45. && TMath::Abs(visTau1P4.eta()) < 2.1 &&
	  visTau2P4.pt() > 45. && TMath::Abs(visTau2P4.eta()) < 2.1) ) {
      acceptance = 1.; // CV: acceptance only, efficiency not taken into account yet
    } else {
      acceptance = 0.;
    }
    //std::cout << "--> returning acceptance = " << acceptance << std::endl;
    return acceptance;
  }
};

double getBR(int decayMode)
{
  if ( decayMode == MeasuredTauLepton::kTauToElecDecay ) {
    return 0.178;
  } else if ( decayMode == MeasuredTauLepton::kTauToMuDecay ) {
    return 0.174;
  } else if ( decayMode == MeasuredTauLepton::kTauToHadDecay ) {
    return 0.648;
  } else {
    std::cerr << "Error: Invalid decayMode = " << decayMode << " !!" << std::endl;
    assert(0);
  }
}

int main(int argc, char* argv[]) 
{
  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 14.e+3; 

  //int mode = HttXsectionIntegrandWithTauDecays::kMadgraph;
  int mode = HttXsectionIntegrandWithTauDecays::kLiterature;

  //int leg1decayMode = MeasuredTauLepton::kTauToMuDecay;
  int leg1decayMode = MeasuredTauLepton::kTauToHadDecay;
  double leg1Mass = 0.;
  if      ( leg1decayMode == MeasuredTauLepton::kTauToElecDecay ) leg1Mass = svFitMEM::electronMass;
  else if ( leg1decayMode == MeasuredTauLepton::kTauToMuDecay   ) leg1Mass = svFitMEM::muonMass;
  else if ( leg1decayMode == MeasuredTauLepton::kTauToHadDecay  ) leg1Mass = svFitMEM::rhoMesonMass;
  else assert(0);

  //int leg2decayMode = MeasuredTauLepton::kTauToMuDecay;
  int leg2decayMode = MeasuredTauLepton::kTauToHadDecay;
  double leg2Mass = 0.;
  if      ( leg2decayMode == MeasuredTauLepton::kTauToElecDecay ) leg2Mass = svFitMEM::electronMass;
  else if ( leg2decayMode == MeasuredTauLepton::kTauToMuDecay   ) leg2Mass = svFitMEM::muonMass;
  else if ( leg2decayMode == MeasuredTauLepton::kTauToHadDecay  ) leg2Mass = svFitMEM::rhoMesonMass;
  else assert(0);
  
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
  branchingRatios[450.] = 1.99e-05;
  branchingRatios[500.] = 1.53e-05;
  for ( std::map<double, double>::iterator branchingRatio = branchingRatios.begin();
	branchingRatio != branchingRatios.end(); ++branchingRatio ) {
    branchingRatio->second *= (getBR(leg1decayMode)*getBR(leg2decayMode));
  }

  std::map<double, double> xSection_targets_13TeV;
  // format: parton luminosity ratio * NLO cross-section [pb] / k factor
  xSection_targets_13TeV[90.]  = 2.137*36.23*branchingRatios[90.]/1.89;	
  xSection_targets_13TeV[100.] = 2.185*29.68*branchingRatios[100.]/1.89;
  xSection_targets_13TeV[150.] = 2.399*13.55*branchingRatios[150.]/1.91;
  xSection_targets_13TeV[160.] = 2.439*11.96*branchingRatios[160.]/1.92;
  xSection_targets_13TeV[200.] = 2.592*7.081*branchingRatios[200.]/1.94;
  xSection_targets_13TeV[250.] = 2.776*4.783*branchingRatios[250.]/1.96;
  xSection_targets_13TeV[300.] = 2.956*3.594*branchingRatios[300.]/1.99;

  std::map<double, double> xSection_targets_14TeV;
  // format: parton luminosity ratio * NLO cross-section [pb] / k factor
  xSection_targets_14TeV[90.]  = 2.385*36.23*branchingRatios[90.]/1.89;	
  xSection_targets_14TeV[100.] = 2.446*29.68*branchingRatios[100.]/(2.446*29.68/25.759);
  xSection_targets_14TeV[125.] = 2.587*19.27*branchingRatios[125.]/(2.587*19.27/17.317);
  xSection_targets_14TeV[150.] = 2.720*13.55*branchingRatios[150.]/(2.720*13.55/12.448);
  xSection_targets_14TeV[160.] = 2.771*11.96*branchingRatios[160.]/1.92;
  xSection_targets_14TeV[200.] = 2.969*7.081*branchingRatios[200.]/(2.969*7.081/7.366);
  xSection_targets_14TeV[250.] = 3.208*4.783*branchingRatios[250.]/(3.208*4.783/4.993);
  xSection_targets_14TeV[300.] = 3.444*3.594*branchingRatios[300.]/(3.444*3.594/3.829);
  xSection_targets_14TeV[350.] = 3.208*4.783*branchingRatios[350.]/(3.208*4.783/3.532);
  xSection_targets_14TeV[400.] = 3.444*3.594*branchingRatios[400.]/(3.444*3.594/3.786);
  xSection_targets_14TeV[450.] = 3.208*4.783*branchingRatios[450.]/(3.208*4.783/2.846);
  xSection_targets_14TeV[500.] = 3.444*3.594*branchingRatios[500.]/(3.444*3.594/1.960);

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

  // define MET covariance
  TMatrixD covMET(2,2);
  covMET[0][0] = 100.00;
  covMET[1][0] =   0.00;
  covMET[0][1] =   0.00;
  covMET[1][1] = 100.00;

  // CV: remove ".gz" suffix, as it is internally added by LHAPDF 
  //std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/cteq65.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();
  std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/MSTW2008lo68cl.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();

  int verbosity = 0;

  TGraphAsymmErrors* graph_Xsection = new TGraphAsymmErrors(mH.size());
  TGraphAsymmErrors* graph_Acc = new TGraphAsymmErrors(mH.size());

  int idxPoint = 0;
  for ( std::vector<double>::const_iterator mH_i = mH.begin();
	mH_i != mH.end(); ++mH_i ) {
    std::cout << "computing cross-sections for mH = " << (*mH_i) << "..." << std::endl;
    std::cout << "without acceptance cuts:" << std::endl;    
    HttXsectionWithTauDecays HttXsection_woAcc(sqrtS, *mH_i, pdfFileName.data(), mode, findFile(madgraphFileNames[*mH_i]), verbosity);
    HttXsection_woAcc.setBR(branchingRatios[*mH_i]);
    HttXsection_woAcc.disableAcceptanceCuts();
    HttXsection_woAcc.setMaxObjFunctionCalls(5000000); 
    //HttXsection_woAcc.setMaxObjFunctionCalls(100000); 
    HttXsection_woAcc.integrate(leg1decayMode, -1, leg1Mass, leg2decayMode, -1, leg2Mass, covMET);
    double xSection = HttXsection_woAcc.xSection();
    double xSectionErr = HttXsection_woAcc.xSectionErr();
    std::cout << "mH = " << (*mH_i) << " (GammaH = " << 1.e-2*(*mH_i) << "):" 
	      << " cross-section = " << xSection << " +/- " << xSectionErr << " pb" 
	      << " (expected = " << xSection_targets_14TeV[*mH_i] << " pb, ratio = " << xSection/xSection_targets_14TeV[*mH_i] << ")" << std::endl;
    //continue;
    graph_Xsection->SetPoint(idxPoint, *mH_i, xSection);
    graph_Xsection->SetPointError(idxPoint, 0., 0., xSectionErr, xSectionErr);
/*
    std::cout << "with acceptance cuts:" << std::endl;
    HttXsectionWithTauDecays HttXsection_wAcc(sqrtS, *mH_i, pdfFileName.data(), mode, findFile(madgraphFileNames[*mH_i]), verbosity);
    HttXsection_wAcc.setBR(branchingRatios[*mH_i]);
    acceptanceTypeHadHad acceptanceHadHad;
    HttXsection_wAcc.enableAcceptanceCuts(acceptanceHadHad);
    //HttXsection_wAcc.setMaxObjFunctionCalls(500000); 
    HttXsection_wAcc.setMaxObjFunctionCalls(10000); 
    HttXsection_wAcc.integrate(leg1decayMode, -1, leg1Mass, leg2decayMode, -1, leg2Mass, covMET);
    double xSection_times_Acc = HttXsection_wAcc.xSection();
    double xSection_times_AccErr = HttXsection_wAcc.xSectionErr();
    std::cout << "cross-section*acceptance = " << xSection_times_Acc << " +/- " << xSection_times_AccErr << " pb" << std::endl;
    double Acc = xSection_times_Acc/xSection;
    double AccErr = Acc*TMath::Sqrt(square(xSection_times_AccErr/xSection_times_Acc) + square(xSectionErr/xSection));
    graph_Acc->SetPoint(idxPoint, *mH_i, Acc);
    std::cout << "acceptance = " << Acc << " +/- " << AccErr << std::endl;    
 */
    ++idxPoint;
  }

  TFile* outputFile = new TFile("testHttXsectionWithTauDecays_hadhad.root", "RECREATE");
  graph_Xsection->Write();
  graph_Acc->Write();
  delete outputFile;

  return 0;
}
