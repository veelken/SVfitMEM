
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
    if ( !inputFile.isLocal() ) {
      std::cerr << "Error: Cannot find file = " << fileName << " !!" << std::endl;
      assert(0);
    }
    return inputFile.fullPath().data();
  }
}

const double pb = 1.e-12; // conversion factor for cross-section to picobarn

double acceptanceTauTau(const svFitMEM::LorentzVector& visTau1P4, const svFitMEM::LorentzVector& visTau2P4, double metPx, double metPy)
{
  //std::cout << "<acceptanceTauTau>:" << std::endl;
  //std::cout << " visTau1: Pt = " << visTau1P4.pt() << ", eta = " << visTau1P4.eta() << ", phi = " << visTau1P4.phi() << ", mass = " << visTau1P4.mass() << std::endl;
  //std::cout << " visTau1: Pt = " << visTau2P4.pt() << ", eta = " << visTau2P4.eta() << ", phi = " << visTau2P4.phi() << ", mass = " << visTau2P4.mass() << std::endl;
  double acceptance = 0.;
  if ( (visTau1P4.pt() > 45. && TMath::Abs(visTau1P4.eta()) < 2.3 &&
	visTau2P4.pt() > 45. && TMath::Abs(visTau2P4.eta()) < 2.3) ) {
    acceptance = 1.; // CV: acceptance only, efficiency not taken into account yet
  } else {
    acceptance = 0.;
  }
  //std::cout << "--> returning acceptance = " << acceptance << std::endl;
  return acceptance;
}

int main(int argc, char* argv[]) 
{
  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 13.e+3; 

  // define Higgs mass points 
  std::vector<double> mH;
  mH.push_back(90.);
  mH.push_back(100.);
  mH.push_back(105.);
  mH.push_back(110.);
  mH.push_back(115.);
  mH.push_back(120.);
  mH.push_back(125.);
  mH.push_back(130.);
  mH.push_back(135.);
  mH.push_back(140.);
  mH.push_back(145.);
  mH.push_back(150.);
  mH.push_back(160.);

  // define MET covariance
  TMatrixD covMET(2,2);
  covMET[0][0] = 100.00;
  covMET[1][0] =   0.00;
  covMET[0][1] =   0.00;
  covMET[1][1] = 100.00;

  std::string madgraphFileName = findFile("TauAnalysis/SVfitMEM/data/param_card.dat");

  // CV: remove ".gz" suffix, as it is internally added by LHAPDF 
  std::string pdfFileName = TString(findFile("TauAnalysis/SVfitMEM/data/cteq65.LHgrid.gz").data()).ReplaceAll(".gz", "").Data();

  int verbosity = 0;

  TGraphAsymmErrors* graph_Xsection = new TGraphAsymmErrors(mH.size());
  TGraphAsymmErrors* graph_Acc = new TGraphAsymmErrors(mH.size());

  int idxPoint = 0;
  for ( std::vector<double>::const_iterator mH_i = mH.begin();
	mH_i != mH.end(); ++mH_i ) {
    std::cout << "computing cross-sections for mH = " << (*mH_i) << "..." << std::endl;
    std::cout << "without acceptance cuts:" << std::endl;
    HttXsectionWithTauDecays HttXsection(madgraphFileName, sqrtS, *mH_i, pdfFileName.data(), verbosity);
    HttXsection.disableAcceptanceCuts();
    HttXsection.integrate(MeasuredTauLepton::kTauToHadDecay, -1, svFitMEM::chargedPionMass, MeasuredTauLepton::kTauToHadDecay, -1, svFitMEM::rhoMesonMass, covMET);
    double xSection = HttXsection.xSection()/pb;
    double xSectionErr = HttXsection.xSectionErr()/pb;
    graph_Xsection->SetPoint(idxPoint, *mH_i, xSection);
    graph_Xsection->SetPointError(idxPoint, 0., 0., xSectionErr, xSectionErr);
    std::cout << "with acceptance cuts:" << std::endl;
    HttXsection.enableAcceptanceCuts(&acceptanceTauTau);
    HttXsection.integrate(MeasuredTauLepton::kTauToHadDecay, -1, svFitMEM::chargedPionMass, MeasuredTauLepton::kTauToHadDecay, -1, svFitMEM::rhoMesonMass, covMET);
    double xSection_times_Acc = HttXsection.xSection()/pb;
    double xSection_times_AccErr = HttXsection.xSectionErr()/pb;
    double Acc = xSection_times_Acc/xSection;
    double AccErr = Acc*TMath::Sqrt(square(xSection_times_AccErr/xSection_times_Acc) + square(xSectionErr/xSection));
    graph_Acc->SetPoint(idxPoint, *mH_i, Acc);
    std::cout << "cross-section = " << xSection << " +/- " << xSectionErr << " pb, acceptance = " << Acc << " +/- " << AccErr << std::endl;    
    ++idxPoint;
  }

  TFile* outputFile = new TFile("testHttXsectionWithTauDecays.root", "RECREATE");
  graph_Xsection->Write();
  graph_Acc->Write();
  delete outputFile;

  return 0;
}
