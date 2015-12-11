
/**
   \class compSVfitMEM_xSection_times_Acc compSVfitMEM_xSection_times_Acc.cc "TauAnalysis/SVfitStandalone/bin/compSVfitMEM_xSection_times_Acc.cc"
   \brief Compute leading order gg -> Higgs -> tautau cross-section with and without acceptance cuts
*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitMEM/interface/HttXsectionWithTauDecays.h"
#include "TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h"
#include "TauAnalysis/SVfitMEM/interface/SVfitMEM.h"

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>

#include <string>
#include <vector>
#include <iostream>

using namespace svFitMEM;

typedef std::vector<double> vdouble;

namespace
{
  std::string format_vdouble(const std::vector<double>& vd)
  {
    std::ostringstream os;
    
    os << "{ ";
    
    unsigned numEntries = vd.size();
    for ( unsigned iEntry = 0; iEntry < numEntries; ++iEntry ) {
      os << vd[iEntry];
      if ( iEntry < (numEntries - 1) ) os << ", ";
    }
    
    os << " }";
    
    return os.str();
  }
}

class acceptanceType : public acceptanceBaseType
{
 public:
  acceptanceType(double minLeg1Pt_1, double minLeg1Pt_2, double maxLeg1AbsEta, double minLeg2Pt_1, double minLeg2Pt_2, double maxLeg2AbsEta)
    : minLeg1Pt_1_(minLeg1Pt_1),
      minLeg1Pt_2_(minLeg1Pt_2),
      maxLeg1AbsEta_(maxLeg1AbsEta),
      minLeg2Pt_1_(minLeg2Pt_1),
      minLeg2Pt_2_(minLeg2Pt_2),
      maxLeg2AbsEta_(maxLeg2AbsEta)
  {}
  ~acceptanceType() {}

  double operator()(const svFitMEM::LorentzVector& visTau1P4, const svFitMEM::LorentzVector& visTau2P4, double metPx, double metPy) const
  {
    //std::cout << "<acceptance>:" << std::endl;
    //std::cout << " visTau1: Pt = " << visTau1P4.pt() << ", eta = " << visTau1P4.eta() << ", phi = " << visTau1P4.phi() << ", mass = " << visTau1P4.mass() << std::endl;
    //std::cout << " visTau2: Pt = " << visTau2P4.pt() << ", eta = " << visTau2P4.eta() << ", phi = " << visTau2P4.phi() << ", mass = " << visTau2P4.mass() << std::endl;
    double acceptance = 0.;
    if ( ((visTau1P4.pt() > minLeg1Pt_1_ && visTau2P4.pt() > minLeg2Pt_2_) ||
	  (visTau1P4.pt() > minLeg1Pt_2_ && visTau2P4.pt() > minLeg2Pt_1_)) &&
	 TMath::Abs(visTau1P4.eta()) < maxLeg1AbsEta_ && TMath::Abs(visTau2P4.eta()) < maxLeg2AbsEta_ ) {
      acceptance = 1.; // CV: acceptance only, efficiency not taken into account yet
    } else {
      acceptance = 0.;
    }
    //std::cout << "--> returning acceptance = " << acceptance << std::endl;
    return acceptance;
  }
  
  double minLeg1Pt_1_;
  double minLeg1Pt_2_;
  double maxLeg1AbsEta_;
  double minLeg2Pt_1_; 
  double minLeg2Pt_2_; 
  double maxLeg2AbsEta_;
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
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<compSVfitMEM_xSection_times_Acc>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("compSVfitMEM_xSection_times_Acc");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("compSVfitMEM_xSection_times_Acc") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgCompSVfitMEM_xSection_times_Acc = cfg.getParameter<edm::ParameterSet>("compSVfitMEM_xSection_times_Acc");
  
  double sqrtS = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("sqrtS");

  int mode = HttXsectionIntegrandWithTauDecays::kLiterature;
  
  std::string leg1decayMode_string = cfgCompSVfitMEM_xSection_times_Acc.getParameter<std::string>("leg1decayMode");
  int leg1decayMode = -1;
  double leg1Mass = 0.;
  if ( leg1decayMode_string == "tauToElecDecay" ) {
    leg1decayMode = MeasuredTauLepton::kTauToElecDecay;
    leg1Mass = svFitMEM::electronMass;
  } else if ( leg1decayMode_string == "tauToMuDecay" ) {
    leg1decayMode = MeasuredTauLepton::kTauToMuDecay;
    leg1Mass = svFitMEM::muonMass;
  } else if ( leg1decayMode_string == "tauToHadDecay" ) {
    leg1decayMode = MeasuredTauLepton::kTauToHadDecay;
    leg1Mass = svFitMEM::chargedPionMass;
  } else throw cms::Exception("compSVfitMEM_xSection_times_Acc") 
    << "Invalid Configuration parameter 'leg1decayMode' = " << leg1decayMode_string << " !!\n";
  double minLeg1Pt_1 = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("minLeg1Pt_1");
  double minLeg1Pt_2 = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("minLeg1Pt_2");
  double maxLeg1AbsEta = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("maxLeg1AbsEta");

  std::string leg2decayMode_string = cfgCompSVfitMEM_xSection_times_Acc.getParameter<std::string>("leg2decayMode");
  int leg2decayMode = -1;
  double leg2Mass = 0.;
  if ( leg2decayMode_string == "tauToElecDecay" ) {
    leg2decayMode = MeasuredTauLepton::kTauToElecDecay;
    leg2Mass = svFitMEM::electronMass;
  } else if ( leg2decayMode_string == "tauToMuDecay" ) {
    leg2decayMode = MeasuredTauLepton::kTauToMuDecay;
    leg2Mass = svFitMEM::muonMass;
  } else if ( leg2decayMode_string == "tauToHadDecay" ) {
    leg2decayMode = MeasuredTauLepton::kTauToHadDecay;
    leg2Mass = svFitMEM::chargedPionMass;
  } else throw cms::Exception("compSVfitMEM_xSection_times_Acc") 
    << "Invalid Configuration parameter 'leg2decayMode' = " << leg2decayMode_string << " !!\n";
  double minLeg2Pt_1 = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("minLeg2Pt_1");
  double minLeg2Pt_2 = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("minLeg2Pt_2");
  double maxLeg2AbsEta = cfgCompSVfitMEM_xSection_times_Acc.getParameter<double>("maxLeg2AbsEta");

  // define Higgs mass points 
  std::vector<double> mH = cfgCompSVfitMEM_xSection_times_Acc.getParameter<vdouble>("mH");
  if ( !(mH.size() >= 1) ) throw cms::Exception("compSVfitMEM_xSection_times_Acc") 
    << "Invalid Configuration parameter 'mH' = " << format_vdouble(mH) << " !!\n";

  double mH_min = mH[0];
  double mH_max = mH[mH.size() - 1];

  double branchingRatio = 1.e-1; // set Higgs -> tautau decay branching fraction to value used by SVfitIntegrand

  // define MET covariance
  TMatrixD covMET(2,2);
  covMET[0][0] = 100.00;
  covMET[1][0] =   0.00;
  covMET[0][1] =   0.00;
  covMET[1][1] = 100.00;

  //std::string pdfName = "cteq66";
  std::string pdfName = "MSTW2008lo68cl";

  unsigned maxObjFunctionCalls = cfgCompSVfitMEM_xSection_times_Acc.getParameter<unsigned>("maxObjFunctionCalls");

  std::string intMode_string = cfgCompSVfitMEM_xSection_times_Acc.getParameter<std::string>("intMode");
  int intMode = 1;
  if      ( intMode_string == "vamp"  ) intMode = SVfitMEM::kVAMP;
  else if ( intMode_string == "vegas" ) intMode = SVfitMEM::kVEGAS;
  else throw cms::Exception("compSVfitMEM_xSection_times_Acc") 
    << "Invalid Configuration parameter 'intMode' = " << intMode_string << " !!\n";

  std::string outputFileName = cfgCompSVfitMEM_xSection_times_Acc.getParameter<std::string>("outputFileName");
  std::cout << " outputFileName = " << outputFileName << std::endl;

  int verbosity = cfgCompSVfitMEM_xSection_times_Acc.getParameter<int>("verbosity");

  TTree* tree = new TTree("tree", "tree");
  Float_t mH_value;
  tree->Branch("mH", &mH_value, "mH/F");
  Int_t leg1decayMode_value, leg2decayMode_value;
  tree->Branch("leg1decayMode", &leg1decayMode_value, "leg1decayMode/I");
  tree->Branch("leg2decayMode", &leg2decayMode_value, "leg2decayMode/I");
  Float_t xSection, xSectionErr;
  tree->Branch("xSection", &xSection, "xSection/F");
  tree->Branch("xSectionErr", &xSectionErr, "xSectionErr/F");
  Float_t xSection_times_BR, xSection_times_BRerr;
  tree->Branch("xSection_times_BR", &xSection_times_BR, "xSection_times_BR/F");
  tree->Branch("xSection_times_BRerr", &xSection_times_BRerr, "xSection_times_BRerr/F");
  Float_t xSection_times_Acc, xSection_times_AccErr;
  tree->Branch("xSection_times_Acc", &xSection_times_Acc, "xSection_times_Acc/F");
  tree->Branch("xSection_times_AccErr", &xSection_times_AccErr, "xSection_times_AccErr/F");
  Float_t xSection_times_BR_times_Acc, xSection_times_BR_times_AccErr;
  tree->Branch("xSection_times_BR_times_Acc", &xSection_times_BR_times_Acc, "xSection_times_BR_times_Acc/F");
  tree->Branch("xSection_times_BR_times_AccErr", &xSection_times_BR_times_AccErr, "xSection_times_BR_times_AccErr/F");
  Float_t Acc, AccErr;
  tree->Branch("Acc", &Acc, "Acc/F");
  tree->Branch("AccErr", &AccErr, "AccErr/F");
  
  TGraphErrors* graph_Xsection_woAcc = new TGraphErrors(mH.size()); 
  std::string graphName_Xsection_woAcc = Form("graph_Xsection_woAcc_mH%f1.0to%f1.0_numCalls%u_intMode%i", mH_min, mH_max, maxObjFunctionCalls, intMode);
  graph_Xsection_woAcc->SetName(graphName_Xsection_woAcc.data());
  TGraphErrors* graph_Xsection_times_BR_woAcc = new TGraphErrors(mH.size()); 
  std::string graphName_Xsection_times_BR_woAcc = Form("graph_Xsection_times_BR_woAcc_mH%f1.0to%f1.0_numCalls%u_intMode%i", mH_min, mH_max, maxObjFunctionCalls, intMode);
  graph_Xsection_times_BR_woAcc->SetName(graphName_Xsection_times_BR_woAcc.data());

  TGraphErrors* graph_Xsection_wAcc = new TGraphErrors(mH.size()); 
  std::string graphName_Xsection_wAcc = Form("graph_Xsection_wAcc_mH%f1.0to%f1.0_numCalls%u_intMode%i", mH_min, mH_max, maxObjFunctionCalls, intMode);
  graph_Xsection_wAcc->SetName(graphName_Xsection_wAcc.data());
  TGraphErrors* graph_Xsection_times_BR_wAcc = new TGraphErrors(mH.size()); 
  std::string graphName_Xsection_times_BR_wAcc = Form("graph_Xsection_times_BR_wAcc_mH%f1.0to%f1.0_numCalls%u_intMode%i", mH_min, mH_max, maxObjFunctionCalls, intMode);
  graph_Xsection_times_BR_wAcc->SetName(graphName_Xsection_times_BR_wAcc.data());

  TGraphErrors* graph_Acc = new TGraphErrors(mH.size()); 
  std::string graphName_Acc = Form("graph_Acc_mH%f1.0to%f1.0_numCalls%u_intMode%i", mH_min, mH_max, maxObjFunctionCalls, intMode);
  graph_Acc->SetName(graphName_Acc.data());

  int idxPoint = 0;
  for ( std::vector<double>::const_iterator mH_i = mH.begin();
	mH_i != mH.end(); ++mH_i ) {
    std::cout << "computing cross-sections for mH = " << (*mH_i) << ", GammaH = " << 1.e-2*(*mH_i) << " (intMode = " << intMode << ", numCalls = " << maxObjFunctionCalls << "):" << std::endl;
    
    mH_value = (*mH_i);

    leg1decayMode_value = leg1decayMode;
    leg2decayMode_value = leg2decayMode;

    std::cout << "without acceptance cuts:" << std::endl;    
    HttXsectionWithTauDecays HttXsection_woAcc(sqrtS, *mH_i, pdfName.data(), mode, "", verbosity);
    HttXsection_woAcc.setBR(branchingRatio);
    HttXsection_woAcc.disableAcceptanceCuts();
    HttXsection_woAcc.setMaxObjFunctionCalls(maxObjFunctionCalls);
    HttXsection_woAcc.setIntMode(intMode);
    HttXsection_woAcc.integrate(leg1decayMode, -1, leg1Mass, leg2decayMode, -1, leg2Mass, covMET);
    xSection_times_BR = HttXsection_woAcc.xSection();
    xSection_times_BRerr = HttXsection_woAcc.xSectionErr();    
    std::cout << " cross-section = " << xSection_times_BR << " +/- " << xSection_times_BRerr << " pb" << std::endl;
    xSection = xSection_times_BR/branchingRatio;
    xSectionErr = xSection_times_BRerr/branchingRatio;
    graph_Xsection_woAcc->SetPoint(idxPoint, *mH_i, xSection);
    graph_Xsection_woAcc->SetPointError(idxPoint, 0., xSectionErr);
    graph_Xsection_times_BR_woAcc->SetPoint(idxPoint, *mH_i, xSection_times_BR);
    graph_Xsection_times_BR_woAcc->SetPointError(idxPoint, 0., xSection_times_BRerr);
    
    std::cout << "with acceptance cuts:" << std::endl;
    HttXsectionWithTauDecays HttXsection_wAcc(sqrtS, *mH_i, pdfName.data(), mode, "", verbosity);
    HttXsection_wAcc.setBR(branchingRatio);
    acceptanceType acceptance(minLeg1Pt_1, minLeg1Pt_2, maxLeg1AbsEta, minLeg2Pt_1, minLeg2Pt_2, maxLeg2AbsEta);
    HttXsection_wAcc.enableAcceptanceCuts(acceptance);
    HttXsection_wAcc.setMaxObjFunctionCalls(maxObjFunctionCalls);
    HttXsection_wAcc.setIntMode(intMode);
    HttXsection_wAcc.integrate(leg1decayMode, -1, leg1Mass, leg2decayMode, -1, leg2Mass, covMET);
    xSection_times_BR_times_Acc = HttXsection_wAcc.xSection();
    xSection_times_BR_times_AccErr = HttXsection_wAcc.xSectionErr();
    std::cout << " cross-section*acceptance = " << xSection_times_BR_times_Acc << " +/- " << xSection_times_BR_times_AccErr << " pb" << std::endl;
    xSection_times_Acc = xSection_times_BR_times_Acc/branchingRatio;
    xSection_times_AccErr = xSection_times_BR_times_AccErr/branchingRatio;
    graph_Xsection_wAcc->SetPoint(idxPoint, *mH_i, xSection_times_Acc);
    graph_Xsection_wAcc->SetPointError(idxPoint, 0., xSection_times_AccErr);
    graph_Xsection_times_BR_wAcc->SetPoint(idxPoint, *mH_i, xSection_times_BR_times_Acc);
    graph_Xsection_times_BR_wAcc->SetPointError(idxPoint, 0., xSection_times_BR_times_AccErr);
    
    Acc = xSection_times_BR_times_Acc/xSection_times_BR;
    AccErr = Acc*TMath::Sqrt(square(xSection_times_BR_times_AccErr/xSection_times_BR_times_Acc) + square(xSection_times_BRerr/xSection_times_BR));
    graph_Acc->SetPoint(idxPoint, *mH_i, Acc);
    graph_Acc->SetPointError(idxPoint, 0., AccErr);
    std::cout << "acceptance = " << Acc << " +/- " << AccErr << std::endl;    

    tree->Fill();

    ++idxPoint;
  }

  TFile* outputFile = new TFile(outputFileName.data(), "RECREATE");
  tree->Write();
  graph_Xsection_woAcc->Write();
  graph_Xsection_times_BR_woAcc->Write();
  graph_Xsection_wAcc->Write();
  graph_Xsection_times_BR_wAcc->Write();
  graph_Acc->Write();
  delete outputFile;

  return 0;
}
