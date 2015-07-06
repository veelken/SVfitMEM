
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>


TH1* loadHistogram(TFile* inputFile, const std::string& dqmDirectory, const std::string& histogramName)
{
  TString histogramName_full = dqmDirectory.data();
  if ( histogramName_full.Length() > 0 && !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(histogramName.data());

  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName_full.Data()));
  std::cout << "histogram = " << histogram << std::endl;
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName_full << " from file = " << inputFile->GetName() << " !!" << std::endl;
    delete inputFile;
    assert(0);
  }
  // if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  //if ( histogram->Integral() > 0. ) histogram->Scale(1./histogram->Integral());
  return histogram;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram, 
		    double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.15);
  canvas->SetBottomMargin(0.15);
  canvas->SetLogy(useLogScale);

  histogram->SetTitle("");
  histogram->SetStats(true);
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);
  histogram->SetLineColor(1);
  histogram->SetLineWidth(2);
  histogram->SetMarkerColor(1);
  histogram->SetMarkerStyle(20);
  histogram->SetMarkerSize(1.5);
  histogram->Draw("hist");

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetRangeUser(xMin, xMax);
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.060);
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetLabelSize(0.050);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.060);
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetLabelSize(0.050);
  yAxis->SetNdivisions(505);

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  //canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete canvas;  
}

void makeTFNormalizationPlots()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_5_3_14/src/TauAnalysis/SVfitMEM/";
  std::string inputFileName = "testTFNormalization.root";
  TString inputFileName_full = inputFilePath;
  if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
  inputFileName_full.Append(inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.Data());

  TH1* histogram_norm_met    = loadHistogram(inputFile, "", "norm_met");
  TH1* histogram_norm_lepTau = loadHistogram(inputFile, "", "norm_lepTau");
  TH1* histogram_norm_hadTau = loadHistogram(inputFile, "", "norm_hadTau");

  showHistograms(800, 600, 
		 histogram_norm_met,
		 0.99, 1.01, "Normalization", 1.15,
		 true, 1.e0, 2.9e+3, "Toys", 1.20,
		 "plots/norm_met.pdf");
  showHistograms(800, 600, 
		 histogram_norm_lepTau,
		 0.99, 1.01, "Normalization", 1.15,
		 true, 1.e0, 2.9e+3, "Toys", 1.20,
		 "plots/norm_lepTau.pdf");
  showHistograms(800, 600, 
		 histogram_norm_hadTau,
		 0.99, 1.01, "Normalization", 1.15,
		 true, 1.e0, 2.9e+3, "Toys", 1.20,
		 "plots/norm_hadTau.pdf");

  delete inputFile;
}


