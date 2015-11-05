
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TROOT.h>

#include <vector>
#include <iostream>
#include <iomanip>

#include "assert.h"

TTree* loadTree(const std::string& inputFileName, std::vector<TFile*>& inputFilesToClose, const std::string& treeName = "tree")
{
  TFile* inputFile = new TFile(inputFileName.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName << "' !!" << std::endl;
    assert(0);
  }
  
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get(treeName.data()));
  if ( !tree ) {
    std::cerr << "Failed to load tree = '" << treeName << "' from file = '" << inputFileName.data() << "' !!" << std::endl;
    inputFile->ls();
    assert(0);
  }
  
  inputFilesToClose.push_back(inputFile);

  return tree;
}

struct entryType
{
  double mH_;
  double xSection_;
  double xSectionErr_;
  double xSection_times_BR_;
  double xSection_times_BRerr_;
  double xSection_times_Acc_;
  double xSection_times_AccErr_;
  double xSection_times_BR_times_Acc_;
  double xSection_times_BR_times_AccErr_;
  double Acc_;
  double AccErr_;
};

bool isHigherMass(const entryType& entry1, const entryType& entry2)
{
  return (entry1.mH_ < entry2.mH_);
}

std::vector<TGraphErrors*> processTree(TTree* tree, int selLeg1decayMode, int selLeg2decayMode, const std::string& sqrtS)
{
  Float_t mH;
  tree->SetBranchAddress("mH", &mH);
  Int_t leg1decayMode, leg2decayMode;
  tree->SetBranchAddress("leg1decayMode", &leg1decayMode);
  tree->SetBranchAddress("leg2decayMode", &leg2decayMode);
  Float_t xSection, xSectionErr;
  tree->SetBranchAddress("xSection", &xSection);
  tree->SetBranchAddress("xSectionErr", &xSectionErr);
  Float_t xSection_times_BR, xSection_times_BRerr;
  tree->SetBranchAddress("xSection_times_BR", &xSection_times_BR);
  tree->SetBranchAddress("xSection_times_BRerr", &xSection_times_BRerr);
  Float_t xSection_times_Acc, xSection_times_AccErr;
  tree->SetBranchAddress("xSection_times_Acc", &xSection_times_Acc);
  tree->SetBranchAddress("xSection_times_AccErr", &xSection_times_AccErr);
  Float_t xSection_times_BR_times_Acc, xSection_times_BR_times_AccErr;
  tree->SetBranchAddress("xSection_times_BR_times_Acc", &xSection_times_BR_times_Acc);
  tree->SetBranchAddress("xSection_times_BR_times_AccErr", &xSection_times_BR_times_AccErr);
  Float_t Acc, AccErr;
  tree->SetBranchAddress("Acc", &Acc);
  tree->SetBranchAddress("AccErr", &AccErr);

  std::vector<entryType> entries;

  int numEntries = tree->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
    if ( iEntry > 0 && (iEntry % 1000) == 0 ) {
      std::cout << "processing Entry " << iEntry << std::endl;
    }

    tree->GetEntry(iEntry);

    if ( !(leg1decayMode == selLeg1decayMode && leg2decayMode == selLeg2decayMode) ) continue;

    entryType entry;
    entry.mH_ = mH;
    entry.xSection_ = xSection;
    entry.xSectionErr_ = xSectionErr;
    entry.xSection_times_BR_ = xSection_times_BR;
    entry.xSection_times_BRerr_ = xSection_times_BRerr;
    entry.xSection_times_Acc_ = xSection_times_Acc;
    entry.xSection_times_AccErr_ = xSection_times_AccErr;
    entry.xSection_times_BR_times_Acc_ = xSection_times_BR_times_Acc;
    entry.xSection_times_BR_times_AccErr_ = xSection_times_BR_times_AccErr;
    entry.Acc_ = Acc;
    entry.AccErr_ = AccErr;

    entries.push_back(entry);
  }

  std::sort(entries.begin(), entries.end(), isHigherMass);

  // WARNING: mapping of selLeg1decayMode and selLeg2decayMode to decayMode_string needs to match definition in
  //          TauAnalysis/SVfitMEM/interface/MeasuredTauLepton.h !!
  std::string decayMode_string = "";
  if      ( selLeg1decayMode == 1 && selLeg2decayMode == 1 ) decayMode_string = "hadhad";
  else if ( selLeg1decayMode != 1 && selLeg2decayMode == 1 ) decayMode_string = "lephad";
  else if ( selLeg1decayMode == 1 && selLeg2decayMode != 1 ) decayMode_string = "hadlep";
  else if ( selLeg1decayMode != 1 && selLeg2decayMode != 1 ) decayMode_string = "leplep";
  else assert(0);

  // WARNING: mapping of selIntMode to intMode_string needs to match definition in
  //          TauAnalysis/SVfitMEM/interface/SVfitMEM.h !!
  //std::string intMode_string;
  //if      ( selIntMode == 1 ) intMode_string = "vegas";
  //else if ( selIntMode == 2 ) intMode_string = "vamp";
  //else assert(0);
  std::string intMode_string = "vamp";

  TGraphErrors* graph_Xsection_woAcc = new TGraphErrors(entries.size()); 
  std::string graphName_Xsection_woAcc = Form("graph_Xsection_woAcc_%s_%s_%s", sqrtS.data(), decayMode_string.data(), intMode_string.data());
  graph_Xsection_woAcc->SetName(graphName_Xsection_woAcc.data());
  TGraphErrors* graph_Xsection_times_BR_woAcc = new TGraphErrors(entries.size()); 
  std::string graphName_Xsection_times_BR_woAcc = Form("graph_Xsection_times_BR_woAcc_%s_%s_%s", sqrtS.data(), decayMode_string.data(), intMode_string.data());
  graph_Xsection_times_BR_woAcc->SetName(graphName_Xsection_times_BR_woAcc.data());
  
  TGraphErrors* graph_Xsection_wAcc = new TGraphErrors(entries.size()); 
  std::string graphName_Xsection_wAcc = Form("graph_Xsection_wAcc_%s_%s_%s", sqrtS.data(), decayMode_string.data(), intMode_string.data());
  graph_Xsection_wAcc->SetName(graphName_Xsection_wAcc.data());
  TGraphErrors* graph_Xsection_times_BR_wAcc = new TGraphErrors(entries.size()); 
  std::string graphName_Xsection_times_BR_wAcc = Form("graph_Xsection_times_BR_wAcc_%s_%s_%s", sqrtS.data(), decayMode_string.data(), intMode_string.data());
  graph_Xsection_times_BR_wAcc->SetName(graphName_Xsection_times_BR_wAcc.data());
  TGraphErrors* graph_Acc = new TGraphErrors(entries.size()); 
  std::string graphName_Acc = Form("graph_Acc_%s_%s_%s", sqrtS.data(), decayMode_string.data(), intMode_string.data());
  graph_Acc->SetName(graphName_Acc.data());

  int idxPoint = 0;
  for ( std::vector<entryType>::const_iterator entry = entries.begin();
	entry != entries.end(); ++entry ) {
    graph_Xsection_woAcc->SetPoint(idxPoint, entry->mH_, entry->xSection_);
    graph_Xsection_woAcc->SetPointError(idxPoint, 0., entry->xSectionErr_);
    graph_Xsection_times_BR_woAcc->SetPoint(idxPoint, entry->mH_, entry->xSection_times_BR_);
    graph_Xsection_times_BR_woAcc->SetPointError(idxPoint, 0., entry->xSection_times_BRerr_);
    
    graph_Xsection_wAcc->SetPoint(idxPoint, entry->mH_, entry->xSection_times_Acc_);
    graph_Xsection_wAcc->SetPointError(idxPoint, 0., entry->xSection_times_AccErr_);
    graph_Xsection_times_BR_wAcc->SetPoint(idxPoint, entry->mH_, entry->xSection_times_BR_times_Acc_);
    graph_Xsection_times_BR_wAcc->SetPointError(idxPoint, 0., entry->xSection_times_BR_times_AccErr_);
	
    graph_Acc->SetPoint(idxPoint, entry->mH_, entry->Acc_);
    graph_Acc->SetPointError(idxPoint, 0., entry->AccErr_);
    
    ++idxPoint;
  }

  std::vector<TGraphErrors*> graphs;
  graphs.push_back(graph_Xsection_woAcc);
  graphs.push_back(graph_Xsection_times_BR_woAcc);
  graphs.push_back(graph_Xsection_wAcc);
  graphs.push_back(graph_Xsection_times_BR_wAcc);
  graphs.push_back(graph_Acc);
  
  return graphs;
}

void makeSVfitMEM_xSection_times_AccGraphs()
{
  gROOT->SetBatch(true);

  std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/SVfitMEM_with_vamp/CMSSW_7_4_6/src/TauAnalysis/SVfitMEM/macros/data/";

  std::map<int, std::map<int, std::map<std::string, std::string> > > inputFileNames; // key = selLeg1decayMode, selLeg2decayMode, sqrtS
  inputFileNames[1][1]["8TeV"] = "compSVfitMEM_xSection_times_Acc_8TeV_hadhad_vamp_all.root";
  inputFileNames[1][2]["8TeV"] = "compSVfitMEM_xSection_times_Acc_8TeV_hadlep_vamp_all.root";
  inputFileNames[2][1]["8TeV"] = "compSVfitMEM_xSection_times_Acc_8TeV_lephad_vamp_all.root";
  inputFileNames[2][2]["8TeV"] = "compSVfitMEM_xSection_times_Acc_8TeV_leplep_vamp_all.root";
  inputFileNames[1][1]["13TeV"] = "compSVfitMEM_xSection_times_Acc_13TeV_hadhad_vamp_all.root";
  inputFileNames[1][2]["13TeV"] = "compSVfitMEM_xSection_times_Acc_13TeV_hadlep_vamp_all.root";
  inputFileNames[2][1]["13TeV"] = "compSVfitMEM_xSection_times_Acc_13TeV_lephad_vamp_all.root";
  inputFileNames[2][2]["13TeV"] = "compSVfitMEM_xSection_times_Acc_13TeV_leplep_vamp_all.root";

  std::vector<std::string> center_of_mass_energies;
  center_of_mass_energies.push_back("8TeV");
  center_of_mass_energies.push_back("13TeV");

  for ( std::vector<std::string>::const_iterator sqrtS = center_of_mass_energies.begin(); sqrtS != center_of_mass_energies.end(); ++sqrtS ) {
    
    std::vector<TFile*> inputFilesToClose;
    
    std::vector<TGraphErrors*> graphs;
    
    for ( int selLeg1decayMode = 1; selLeg1decayMode <= 2; ++selLeg1decayMode ) {
      for ( int selLeg2decayMode = 1; selLeg2decayMode <= 2; ++selLeg2decayMode ) {
	assert(inputFileNames[selLeg1decayMode][selLeg2decayMode][*sqrtS] != "");
	std::string inputFileName_job = std::string(inputFilePath).append(inputFileNames[selLeg1decayMode][selLeg2decayMode][*sqrtS]);
	TTree* tree = loadTree(inputFileName_job, inputFilesToClose);
	std::vector<TGraphErrors*> graphs_job = processTree(tree, selLeg1decayMode, selLeg2decayMode, *sqrtS);
	graphs.insert(graphs.end(), graphs_job.begin(), graphs_job.end());
      }
    }
    
    std::string outputFileName = Form("makeSVfitMEM_xSection_times_AccGraphs_%s.root", sqrtS->data());
    TFile* outputFile = new TFile(outputFileName.data(), "RECREATE");
    for ( std::vector<TGraphErrors*>::iterator graph = graphs.begin();
	  graph != graphs.end(); ++graph ) {
      (*graph)->Write();
    }
    delete outputFile;
    
    for ( std::vector<TFile*>::iterator inputFile = inputFilesToClose.begin();
	  inputFile != inputFilesToClose.end(); ++inputFile ) {
      delete (*inputFile);
    }
    
    for ( std::vector<TGraphErrors*>::iterator graph = graphs.begin();
	  graph != graphs.end(); ++graph ) {
      delete (*graph);
    }
  }
}
