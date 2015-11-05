
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TGraph.h>
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

std::vector<double> extractVarFromTree(TTree* tree, const std::string& branchName, double mH)
{
  Float_t genLeg1Pt, genLeg2Pt;
  tree->SetBranchAddress("l1GenPt", &genLeg1Pt);
  tree->SetBranchAddress("l2GenPt", &genLeg2Pt);
  Float_t genHiggsMass;
  tree->SetBranchAddress("genHiggsMass", &genHiggsMass);
  Float_t var;
  tree->SetBranchAddress(branchName.data(), &var);

  std::vector<double> values;

  int numEntries = tree->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
    if ( iEntry > 0 && (iEntry % 1000) == 0 ) {
      std::cout << "processing Event " << iEntry << std::endl;
    }

    tree->GetEntry(iEntry);

    if ( !(genLeg1Pt > 20. && genLeg2Pt > 20.) ) continue; 
    if ( !(genHiggsMass > (0.95*mH) && genHiggsMass < (1.05*mH)) ) continue;

    values.push_back(var);
  }

  return values;
}

//-------------------------------------------------------------------------------
// CV: functions to compute tau pair mass using collinear approximation 
//     copied from TauAnalysis/CandidateTools/src/candidateAuxFunctions.cc

void compX1X2byCollinearApprox(double& x1, double& x2, double pxLeg1, double pyLeg1, double pxLeg2, double pyLeg2, double pxMEt, double pyMEt)
{
  double x1_numerator = pxLeg1*pyLeg2 - pxLeg2*pyLeg1;
  double x1_denominator = pyLeg2*(pxLeg1 + pxMEt) - pxLeg2*(pyLeg1 + pyMEt);
  x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
  //std::cout << "x1 = " << x1 << std::endl;

  double x2_numerator = x1_numerator;
  double x2_denominator = pxLeg1*(pyLeg2 + pyMEt) - pyLeg1*(pxLeg2 + pxMEt);
  x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
  //std::cout << "x2 = " << x2 << std::endl;
}

double getPhysX(double x, bool& isWithinPhysRange)
{
  double physX = x;

  isWithinPhysRange = true;

  if ( x < 0. ) {
    physX = 0.;
    isWithinPhysRange = false;
  }

  if ( x > 1. ) {
    physX = 1.;
    isWithinPhysRange = false;
  }

  return physX;
}
//-------------------------------------------------------------------------------

std::vector<double> extractCollinearApproxMassFromTree(TTree* tree, double mH, int optUnphysicalSolution = 0)
{
  Float_t genLeg1Pt, genLeg2Pt;
  tree->SetBranchAddress("l1GenPt", &genLeg1Pt);
  tree->SetBranchAddress("l2GenPt", &genLeg2Pt);
  Float_t genHiggsMass;
  tree->SetBranchAddress("genHiggsMass", &genHiggsMass);
  Float_t leg1Px, leg1Py;
  tree->SetBranchAddress("l1Px", &leg1Px);
  tree->SetBranchAddress("l1Py", &leg1Py);
  Float_t leg2Px, leg2Py;
  tree->SetBranchAddress("l2Px", &leg2Px);
  tree->SetBranchAddress("l2Py", &leg2Py);
  Float_t metPx, metPy;
  tree->SetBranchAddress("mex", &metPx);
  tree->SetBranchAddress("mey", &metPy);
  Float_t visMass;
  tree->SetBranchAddress("visMass", &visMass);
  
  std::vector<double> values;

  int numEntries = tree->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
    if ( iEntry > 0 && (iEntry % 1000) == 0 ) {
      std::cout << "processing Event " << iEntry << std::endl;
    }

    tree->GetEntry(iEntry);

    if ( !(genLeg1Pt > 20. && genLeg2Pt > 20.) ) continue; 
    if ( !(genHiggsMass > (0.95*mH) && genHiggsMass < (1.05*mH)) ) continue;

    double x1, x2;
    compX1X2byCollinearApprox(x1, x2, leg1Px, leg1Py, leg2Px, leg2Py, metPx, metPy);

    bool isX1withinPhysRange, isX2withinPhysRange;
    double x1phys = getPhysX(x1, isX1withinPhysRange);
    double x2phys = getPhysX(x2, isX2withinPhysRange);
    
    if ( isX1withinPhysRange && x1phys != 0. && 
	 isX2withinPhysRange && x2phys != 0. ) {
      values.push_back(visMass/TMath::Sqrt(x1phys*x2phys));
    } else {
      if ( optUnphysicalSolution == 1 ) values.push_back(1.e-3);
    }
  }

  return values;
}

void fillHistogram1d(TH1* histogram, const std::vector<double>& values)
{
  TAxis* xAxis = histogram->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
  size_t numValues = values.size();
  for ( size_t idxValue = 0; idxValue < numValues; ++idxValue ) {
    double value = values[idxValue];
    if ( value < (1.01*xMin) ) value = 1.01*xMin;
    if ( value > (0.99*xMax) ) value = 0.99*xMax;
    histogram->Fill(value);
  }
}

void fillHistogram2d(TH2* histogram, const std::vector<double>& valuesX, const std::vector<double>& valuesY)
{
  assert(valuesX.size() == valuesY.size());
  TAxis* xAxis = histogram->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
  TAxis* yAxis = histogram->GetYaxis();
  double yMin = yAxis->GetXmin();
  double yMax = yAxis->GetXmax();
  size_t numValues = valuesX.size();
  for ( size_t idxValue = 0; idxValue < numValues; ++idxValue ) {
    double valueX = valuesX[idxValue];
    if ( valueX < (1.01*xMin) ) valueX = 1.01*xMin;
    if ( valueX > (0.99*xMax) ) valueX = 0.99*xMax;
    double valueY = valuesY[idxValue];
    if ( valueY < (1.01*yMin) ) valueY = 1.01*yMin;
    if ( valueY > (0.99*yMax) ) valueY = 0.99*yMax;
    histogram->Fill(valueX, valueY);
  }
}

void fillHistogram2dFromTree(TH2* histogram, TTree* tree1, const std::string& branchName1, TTree* tree2, const std::string& branchName2, double mH)
{
  ULong64_t run_1, event_1, lumi_1;
  tree1->SetBranchAddress("run", &run_1);
  tree1->SetBranchAddress("event", &event_1);
  tree1->SetBranchAddress("lumi", &lumi_1);
  Float_t genLeg1Pt_1, genLeg2Pt_1;
  tree1->SetBranchAddress("l1GenPt", &genLeg1Pt_1);
  tree1->SetBranchAddress("l2GenPt", &genLeg2Pt_1);
  Float_t genHiggsMass_1;
  tree1->SetBranchAddress("genHiggsMass", &genHiggsMass_1);
  Float_t var_1;
  tree1->SetBranchAddress(branchName1.data(), &var_1);

  ULong64_t run_2, event_2, lumi_2;
  tree2->SetBranchAddress("run", &run_2);
  tree2->SetBranchAddress("event", &event_2);
  tree2->SetBranchAddress("lumi", &lumi_2);
  Float_t genLeg1Pt_2, genLeg2Pt_2;
  tree2->SetBranchAddress("l1GenPt", &genLeg1Pt_2);
  tree2->SetBranchAddress("l2GenPt", &genLeg2Pt_2);
  Float_t genHiggsMass_2;
  tree2->SetBranchAddress("genHiggsMass", &genHiggsMass_2);
  Float_t var_2;
  tree2->SetBranchAddress(branchName2.data(), &var_2);

  int numEntries1 = tree1->GetEntries();
  int numEntries2 = tree2->GetEntries();
  for ( int iEntry1 = 0; iEntry1 < numEntries1; ++iEntry1 ) {
    if ( iEntry1 > 0 && (iEntry1 % 100) == 0 ) {
      std::cout << "processing Event " << iEntry1 << std::endl;
    }

    tree1->GetEntry(iEntry1);

    if ( genLeg1Pt_1 > 20. && genLeg2Pt_1 > 20. && genHiggsMass_1 > (0.95*mH) && genHiggsMass_1 < (1.05*mH) ) {

      // CV: find entry that matches by run, ls and event number in tree2
      for ( int iEntry2 = 0; iEntry2 < numEntries2; ++iEntry2 ) {

	tree2->GetEntry(iEntry2);

	if ( run_1 == run_2 && lumi_1 == lumi_2 && event_1 == event_2 ) {
	  if ( genLeg1Pt_2 > 20. && genLeg2Pt_2 > 20. && genHiggsMass_2 > (0.95*mH) && genHiggsMass_2 < (1.05*mH) ) {
	    histogram->Fill(var_1, var_2);
	  }
	  break;
	}
      }
    }
  }
}

void normalizeHistogram(TH1* histogram)
{
  double integral = 0.;
  int numBins = histogram->GetNbinsX();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    integral += histogram->GetBinContent(idxBin);
  }
  if ( integral > 0. ) {
    histogram->Sumw2();
    histogram->Scale(1./integral);
  }
}

TH1* compRatioHistogram(const std::string& ratioHistogramName, const TH1* numerator, const TH1* denominator)
{
  TH1* histogramRatio = 0;
  
  if ( numerator->GetDimension() == denominator->GetDimension() &&
       numerator->GetNbinsX() == denominator->GetNbinsX() ) {
    histogramRatio = (TH1*)numerator->Clone(ratioHistogramName.data());
    histogramRatio->Divide(denominator);
    
    int nBins = histogramRatio->GetNbinsX();
    for ( int iBin = 1; iBin <= nBins; ++iBin ){
      double binContent = histogramRatio->GetBinContent(iBin);
      histogramRatio->SetBinContent(iBin, binContent);
    }
    
    histogramRatio->SetLineColor(numerator->GetLineColor());
    histogramRatio->SetLineWidth(numerator->GetLineWidth());
    histogramRatio->SetMarkerColor(numerator->GetMarkerColor());
    histogramRatio->SetMarkerStyle(numerator->GetMarkerStyle());
    histogramRatio->SetMarkerSize(numerator->GetMarkerSize());
  }

  return histogramRatio;
}

void showHistograms1d(double canvasSizeX, double canvasSizeY,
		      TH1* histogram_ref, const std::string& legendEntry_ref, 
		      TH1* histogram2, const std::string& legendEntry2, 
		      TH1* histogram3, const std::string& legendEntry3, 
		      TH1* histogram4, const std::string& legendEntry4, 
		      TH1* histogram5, const std::string& legendEntry5, 
		      TH1* histogram6, const std::string& legendEntry6, 
		      const std::string& xAxisTitle, double xAxisOffset,
		      bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		      double legendX0, double legendY0, 
		      const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.27, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.04);
  topPad->SetLeftMargin(0.15);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.05);
  topPad->SetLogy(useLogScale);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.27);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  int colors[6] = { 1, 2, 3, 4, 6, 7 };
  int markerStyles[6] = { 24, 25, 20, 21, 22, 23 };
  int markerSizes[6] = { 1, 1, 1, 1, 1, 1 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.70, legendY0 + 0.25, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);

  histogram_ref->SetTitle("");
  histogram_ref->SetStats(false);
  histogram_ref->SetMinimum(yMin);
  histogram_ref->SetMaximum(yMax);
  histogram_ref->SetLineColor(colors[0]);
  histogram_ref->SetLineWidth(1);
  histogram_ref->SetMarkerColor(colors[0]);
  histogram_ref->SetMarkerStyle(markerStyles[0]);
  histogram_ref->SetMarkerSize(markerSizes[0]);
  histogram_ref->Draw("e1p");
  legend->AddEntry(histogram_ref, Form("%s: mean=%1.1f, rms=%1.1f", legendEntry_ref.data(), histogram_ref->GetMean(), histogram_ref->GetRMS()), "p");

  TAxis* xAxis_top = histogram_ref->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogram_ref->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(1);
    histogram2->SetMarkerColor(colors[1]);
    histogram2->SetMarkerStyle(markerStyles[1]);
    histogram2->SetMarkerSize(markerSizes[1]);
    histogram2->Draw("e1psame");
    legend->AddEntry(histogram2, Form("%s: mean=%1.1f, rms=%1.1f", legendEntry2.data(), histogram2->GetMean(), histogram2->GetRMS()), "p");
  }

  if ( histogram3 ) {
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineWidth(1);
    histogram3->SetMarkerColor(colors[2]);
    histogram3->SetMarkerStyle(markerStyles[2]);
    histogram3->SetMarkerSize(markerSizes[2]);
    histogram3->Draw("e1psame");
    legend->AddEntry(histogram3, Form("%s: mean=%1.1f, rms=%1.1f", legendEntry3.data(), histogram3->GetMean(), histogram3->GetRMS()), "p");
  }

  if ( histogram4 ) {
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineWidth(1);
    histogram4->SetMarkerColor(colors[3]);
    histogram4->SetMarkerStyle(markerStyles[3]);
    histogram4->SetMarkerSize(markerSizes[3]);
    histogram4->Draw("e1psame");
    legend->AddEntry(histogram4, Form("%s: mean=%1.1f, rms=%1.1f", legendEntry4.data(), histogram4->GetMean(), histogram4->GetRMS()), "p");
  }

  if ( histogram5 ) {
    histogram5->SetLineColor(colors[4]);
    histogram5->SetLineWidth(1);
    histogram5->SetMarkerColor(colors[4]);
    histogram5->SetMarkerStyle(markerStyles[4]);
    histogram5->SetMarkerSize(markerSizes[4]);
    histogram5->Draw("e1psame");
    legend->AddEntry(histogram5, Form("%s: mean=%1.1f, rms=%1.1f", legendEntry5.data(), histogram5->GetMean(), histogram5->GetRMS()), "p");
  }
  
  if ( histogram6 ) {
    histogram6->SetLineColor(colors[5]);
    histogram6->SetLineWidth(1);
    histogram6->SetMarkerColor(colors[5]);
    histogram6->SetMarkerStyle(markerStyles[5]);
    histogram6->SetMarkerSize(markerSizes[5]);
    histogram6->Draw("e1psame");
    legend->AddEntry(histogram6, Form("%s: mean=%1.1f, rms=%1.1f", legendEntry6.data(), histogram6->GetMean(), histogram6->GetRMS()), "p");
  }

  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogram2_div_ref = 0;
  if ( histogram2 ) {
    std::string histogramName2_div_ref = std::string(histogram2->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram2_div_ref = compRatioHistogram(histogramName2_div_ref, histogram2, histogram_ref);
    if ( histogram2_div_ref ) {
      histogram2_div_ref->SetTitle("");
      histogram2_div_ref->SetStats(false);
      histogram2_div_ref->SetMinimum(-0.10);
      histogram2_div_ref->SetMaximum(+1.90);
      
      TAxis* xAxis_bottom = histogram2_div_ref->GetXaxis();
      xAxis_bottom->SetTitle(xAxis_top->GetTitle());
      xAxis_bottom->SetLabelColor(1);
      xAxis_bottom->SetTitleColor(1);
      xAxis_bottom->SetTitleOffset(1.20);
      xAxis_bottom->SetTitleSize(0.08);
      xAxis_bottom->SetLabelOffset(0.02);
      xAxis_bottom->SetLabelSize(0.08);
      xAxis_bottom->SetTickLength(0.055);
      
      TAxis* yAxis_bottom = histogram2_div_ref->GetYaxis();
      yAxis_bottom->SetTitle("Ratio");
      yAxis_bottom->SetTitleOffset(0.85);
      yAxis_bottom->SetNdivisions(505);
      yAxis_bottom->CenterTitle();
      yAxis_bottom->SetTitleSize(0.08);
      yAxis_bottom->SetLabelSize(0.08);
      yAxis_bottom->SetTickLength(0.04);  
      
      histogram2_div_ref->Draw("e1p");
    }
  }

  TH1* histogram3_div_ref = 0;
  if ( histogram3 ) {
    std::string histogramName3_div_ref = std::string(histogram3->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram3_div_ref = compRatioHistogram(histogramName3_div_ref, histogram3, histogram_ref);
    if ( histogram3_div_ref ) {
      histogram3_div_ref->SetTitle("");
      histogram3_div_ref->SetStats(false);
      histogram3_div_ref->SetMinimum(-0.50);
      histogram3_div_ref->SetMaximum(+0.50);
      
      TAxis* xAxis_bottom = histogram3_div_ref->GetXaxis();
      xAxis_bottom->SetTitle(xAxis_top->GetTitle());
      xAxis_bottom->SetLabelColor(1);
      xAxis_bottom->SetTitleColor(1);
      xAxis_bottom->SetTitleOffset(1.20);
      xAxis_bottom->SetTitleSize(0.08);
      xAxis_bottom->SetLabelOffset(0.02);
      xAxis_bottom->SetLabelSize(0.08);
      xAxis_bottom->SetTickLength(0.055);
      
      TAxis* yAxis_bottom = histogram3_div_ref->GetYaxis();
      yAxis_bottom->SetTitle("Ratio");
      yAxis_bottom->SetTitleOffset(0.85);
      yAxis_bottom->SetNdivisions(505);
      yAxis_bottom->CenterTitle();
      yAxis_bottom->SetTitleSize(0.08);
      yAxis_bottom->SetLabelSize(0.08);
      yAxis_bottom->SetTickLength(0.04);  
      
      if ( histogram2 ) histogram3_div_ref->Draw("e1psame");
      else histogram3_div_ref->Draw("e1p");
    }
  }

  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_top->GetXmin(), 0.);
  graph_line->SetPoint(1, xAxis_top->GetXmax(), 0.);
  graph_line->SetLineColor(8);
  graph_line->SetLineWidth(1);
  graph_line->Draw("L");

  if ( histogram2_div_ref ) histogram2_div_ref->Draw("e1psame");
  if ( histogram3_div_ref ) histogram3_div_ref->Draw("e1psame");

  TH1* histogram4_div_ref = 0;
  if ( histogram4 ) {
    std::string histogramName4_div_ref = std::string(histogram4->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram4_div_ref = compRatioHistogram(histogramName4_div_ref, histogram4, histogram_ref);
    if ( histogram4_div_ref ) {
      histogram4_div_ref->Draw("e1psame");
    }
  }

  TH1* histogram5_div_ref = 0;
  if ( histogram5 ) {
    std::string histogramName5_div_ref = std::string(histogram5->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram5_div_ref = compRatioHistogram(histogramName5_div_ref, histogram5, histogram_ref);
    if ( histogram5_div_ref ) {
      histogram5_div_ref->Draw("e1psame");
    }
  }
  
  TH1* histogram6_div_ref = 0;
  if ( histogram6 ) {
    std::string histogramName6_div_ref = std::string(histogram6->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram6_div_ref = compRatioHistogram(histogramName6_div_ref, histogram6, histogram_ref);
    if ( histogram6_div_ref ) {
      histogram6_div_ref->Draw("e1psame");
    }
  }

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete legend;
  delete histogram2_div_ref;
  delete histogram3_div_ref;
  delete histogram4_div_ref;
  delete histogram5_div_ref;
  delete histogram6_div_ref;
  delete topPad;
  delete bottomPad;
  delete canvas;  
}

void showHistogram2d(double canvasSizeX, double canvasSizeY,
		     TH2* histogram, 
		     const std::string& xAxisTitle, double xAxisOffset,
		     const std::string& yAxisTitle, double yAxisOffset,
		     const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  int colors[1] = { 1 };
  int markerStyles[1] = { 8 };
  int markerSizes[1] = { 1 };

  histogram->SetTitle("");
  histogram->SetStats(false);
  //histogram->SetLineColor(colors[0]);
  //histogram->SetLineWidth(1);
  //histogram->SetMarkerColor(colors[0]);
  //histogram->SetMarkerStyle(markerStyles[0]);
  //histogram->SetMarkerSize(markerSizes[0]);
  histogram->Draw("scat");
      
  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  
  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);

  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis->GetXmin(), yAxis->GetXmin());
  graph_line->SetPoint(1, xAxis->GetXmax(), yAxis->GetXmax());
  graph_line->SetLineColor(8);
  graph_line->SetLineWidth(1);
  graph_line->Draw("L");

  histogram->Draw("scatsame");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete canvas;  
}

void makeSVfitMEM_PerformancePlots()
{
  gROOT->SetBatch(true);

  std::string inputFilePath = "/data2/veelken/CMSSW_5_3_x/Ntuples/AHtoTauTau/2014Jun09/v3_15_addSVfitMEM/";
  //std::string inputFilePath = "/tmp/veelken/";

  std::vector<std::string> samplesToAnalyze;
  samplesToAnalyze.push_back("HiggsSUSYGluGlu120");
  samplesToAnalyze.push_back("HiggsSUSYGluGlu200");
  samplesToAnalyze.push_back("HiggsSUSYGluGlu300");
  samplesToAnalyze.push_back("HiggsSUSYGluGlu500");
  samplesToAnalyze.push_back("HiggsSUSYGluGlu700");
  samplesToAnalyze.push_back("HiggsSUSYGluGlu1000");

  std::vector<std::string> intModes;
  //intModes.push_back("vegas");
  intModes.push_back("vamp");

  std::vector<std::string> xSection_times_Acc;
  //xSection_times_Acc.push_back("none");
  xSection_times_Acc.push_back("xSection");
  //xSection_times_Acc.push_back("xSection_times_Acc");

  std::map<std::string, std::string> branchName_svFitMEM; // key = xSection_times_Acc
  branchName_svFitMEM["none"]               = "svfitMEMuncorrMass";
  branchName_svFitMEM["xSection"]           = "svfitMEMxSectionCorrMass";
  branchName_svFitMEM["xSection_times_Acc"] = "svfitMEMxSection_times_AccCorrMass";
  
  std::vector<int> logM;
  logM.push_back(0);
  logM.push_back(1);
  logM.push_back(2);
  logM.push_back(3);
  logM.push_back(5);
  logM.push_back(6);
  logM.push_back(7);
  logM.push_back(10);
  logM.push_back(15);
  logM.push_back(20);
  logM.push_back(50);
  logM.push_back(100);
  
  std::vector<int> numCalls;
  numCalls.push_back(10000); 
  numCalls.push_back(20000); 
  numCalls.push_back(50000); 
  numCalls.push_back(100000); 

  std::string treeName = "tauTauNtupleProducer/H2TauTauTreeProducerTauTau";

  std::map<std::string, double> mH; // key = sample
  mH["HiggsSUSYGluGlu120"]  =  120.;
  mH["HiggsSUSYGluGlu200"]  =  200.;
  mH["HiggsSUSYGluGlu300"]  =  300.;
  mH["HiggsSUSYGluGlu500"]  =  500.;
  mH["HiggsSUSYGluGlu700"]  =  700.;
  mH["HiggsSUSYGluGlu1000"] = 1000.;

  //-----------------------------------------------------------------------------
  // make 2d correlation plots between SVfit and SVfitMEM mass

  for ( std::vector<std::string>::const_iterator sample = samplesToAnalyze.begin();
	sample != samplesToAnalyze.end(); ++sample ) {

    assert(mH.find(*sample) != mH.end());
    double mH_i = mH[*sample];

    for ( std::vector<std::string>::const_iterator xSection_times_Acc_i = xSection_times_Acc.begin();
	  xSection_times_Acc_i != xSection_times_Acc.end(); ++xSection_times_Acc_i ) {

      assert(branchName_svFitMEM.find(*xSection_times_Acc_i) != branchName_svFitMEM.end());
      std::string branchName_svFitMEM_i = branchName_svFitMEM[*xSection_times_Acc_i];

      for ( std::vector<int>::const_iterator logM_i = logM.begin();
	    logM_i != logM.end(); ++logM_i ) {
	
	if ( (*xSection_times_Acc_i) == "none" && (*logM_i) != 0 ) continue;

	for ( std::vector<int>::const_iterator numCalls_i = numCalls.begin();
	      numCalls_i != numCalls.end(); ++numCalls_i ) {
	  
	  if ( (*numCalls_i) != 100000 && (*logM_i) != 5 && !((*numCalls_i) == 20000 && (*logM_i) == 6) ) continue;
	  if ( (*numCalls_i) !=  20000 && (*logM_i) == 6 ) continue;

	  std::vector<TFile*> inputFilesToClose2d;

	  std::cout << "processing sample = " << (*sample) << ", xSection_times_Acc = " << (*xSection_times_Acc_i) << ", logM = " << (*logM_i) << ", numCalls = " << (*numCalls_i) << "..." << std::endl;
	  
	  std::string inputFileName_vamp;
	  if ( (*logM_i) == 6 ) {
	    inputFileName_vamp = std::string(inputFilePath).append(Form(
              "%s/addLogM%i/vamp/maxObjFunctionCalls%ik/nom/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
	      xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000, sample->data()));
	  } else {
	    inputFileName_vamp = std::string(inputFilePath).append(Form(
              "%s/addLogM%i/vamp/maxObjFunctionCalls%ik/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
	      xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000, sample->data()));
	  }
	  TTree* tree_vamp = loadTree(inputFileName_vamp, inputFilesToClose2d, treeName);	  
	  std::vector<double> svFit = extractVarFromTree(tree_vamp, "svfitMass", mH_i);
	  delete tree_vamp;
	  tree_vamp = loadTree(inputFileName_vamp, inputFilesToClose2d, treeName);
	  std::vector<double> svFitMEM_vamp = extractVarFromTree(tree_vamp, branchName_svFitMEM_i, mH_i);
	  delete tree_vamp;
	  tree_vamp = loadTree(inputFileName_vamp, inputFilesToClose2d, treeName);
	  std::vector<double> mTtotal = extractVarFromTree(tree_vamp, "mTtotal", mH_i);
	  delete tree_vamp;
	  tree_vamp = loadTree(inputFileName_vamp, inputFilesToClose2d, treeName);
	  std::vector<double> visMass = extractVarFromTree(tree_vamp, "visMass", mH_i);
	  delete tree_vamp;
	  tree_vamp = loadTree(inputFileName_vamp, inputFilesToClose2d, treeName);
	  std::vector<double> caMass = extractCollinearApproxMassFromTree(tree_vamp, mH_i, 1);

	  std::string histogramName_svFitMEM_vamp_vs_svFit = Form(
            "histogram_svFitMEM_vamp_vs_svFit_%s_%s_addLogM%i_numCalls%ik", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  TH2* histogram_svFitMEM_vamp_vs_svFit = new TH2D(
            histogramName_svFitMEM_vamp_vs_svFit.data(), 
	    histogramName_svFitMEM_vamp_vs_svFit.data(), 100, 0., 2.*mH_i, 100, 0., 2.*mH_i);
	  fillHistogram2d(histogram_svFitMEM_vamp_vs_svFit, svFit, svFitMEM_vamp);
	  std::string outputFileName_svFitMEM_vamp_vs_svFit = Form(
            "plots/makeSVfitMEM_PerformancePlots_svFitMEM_vamp_vs_svFit_%s_%s_addLogM%i_numCalls%ik.png", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  showHistogram2d(
	    800, 800, histogram_svFitMEM_vamp_vs_svFit, 
	    "m_{#tau#tau}^{svFit} [GeV]", 1.3, "m_{#tau#tau}^{svFitMEM} [GeV]", 1.65, 
	    outputFileName_svFitMEM_vamp_vs_svFit);
	  delete histogram_svFitMEM_vamp_vs_svFit;

	  std::string histogramName_svFitMEM_vamp_vs_mTtotal = Form(
            "histogram_svFitMEM_vamp_vs_mTtotal_%s_%s_addLogM%i_numCalls%ik", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  TH2* histogram_svFitMEM_vamp_vs_mTtotal = new TH2D(
            histogramName_svFitMEM_vamp_vs_mTtotal.data(), 
	    histogramName_svFitMEM_vamp_vs_mTtotal.data(), 100, 0., 2.*mH_i, 100, 0., 2.*mH_i);
	  fillHistogram2d(histogram_svFitMEM_vamp_vs_mTtotal, mTtotal, svFitMEM_vamp);
	  std::string outputFileName_svFitMEM_vamp_vs_mTtotal = Form(
            "plots/makeSVfitMEM_PerformancePlots_svFitMEM_vamp_vs_mTtotal_%s_%s_addLogM%i_numCalls%ik.png", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  showHistogram2d(
	    800, 800, histogram_svFitMEM_vamp_vs_mTtotal, 
	    "m_{T}^{total} [GeV]", 1.3, "m_{#tau#tau}^{svFitMEM} [GeV]", 1.65, 
	    outputFileName_svFitMEM_vamp_vs_mTtotal);
	  delete histogram_svFitMEM_vamp_vs_mTtotal;

	  std::string histogramName_svFitMEM_vamp_vs_visMass = Form(
            "histogram_svFitMEM_vamp_vs_visMass_%s_%s_addLogM%i_numCalls%ik", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  TH2* histogram_svFitMEM_vamp_vs_visMass = new TH2D(
            histogramName_svFitMEM_vamp_vs_visMass.data(), 
	    histogramName_svFitMEM_vamp_vs_visMass.data(), 100, 0., 2.*mH_i, 100, 0., 2.*mH_i);
	  fillHistogram2d(histogram_svFitMEM_vamp_vs_visMass, visMass, svFitMEM_vamp);
	  std::string outputFileName_svFitMEM_vamp_vs_visMass = Form(
            "plots/makeSVfitMEM_PerformancePlots_svFitMEM_vamp_vs_visMass_%s_%s_addLogM%i_numCalls%ik.png", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  showHistogram2d(
	    800, 800, histogram_svFitMEM_vamp_vs_visMass, 
	    "m_{vis} [GeV]", 1.3, "m_{#tau#tau}^{svFitMEM} [GeV]", 1.65, 
	    outputFileName_svFitMEM_vamp_vs_visMass);
	  delete histogram_svFitMEM_vamp_vs_visMass;

	  std::string histogramName_svFitMEM_vamp_vs_caMass = Form(
            "histogram_svFitMEM_vamp_vs_caMass_%s_%s_addLogM%i_numCalls%ik", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  TH2* histogram_svFitMEM_vamp_vs_caMass = new TH2D(
            histogramName_svFitMEM_vamp_vs_caMass.data(), 
	    histogramName_svFitMEM_vamp_vs_caMass.data(), 100, 0., 2.*mH_i, 100, 0., 2.*mH_i);
	  fillHistogram2d(histogram_svFitMEM_vamp_vs_caMass, caMass, svFitMEM_vamp);
	  std::string outputFileName_svFitMEM_vamp_vs_caMass = Form(
            "plots/makeSVfitMEM_PerformancePlots_svFitMEM_vamp_vs_caMass_%s_%s_addLogM%i_numCalls%ik.png", 
	    sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	  showHistogram2d(
	    800, 800, histogram_svFitMEM_vamp_vs_caMass, 
	    "m_{#tau#tau}^{CA} [GeV]", 1.3, "m_{#tau#tau}^{svFitMEM} [GeV]", 1.65, 
	    outputFileName_svFitMEM_vamp_vs_caMass);
	  delete histogram_svFitMEM_vamp_vs_caMass;
/*
	  if ( (*numCalls_i) != 100000 ) {
	    std::string inputFileName_vamp_100k = std::string(inputFilePath).append(Form(
              "%s/addLogM%i/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
	      xSection_times_Acc_i->data(), *logM_i, sample->data()));
	    TTree* tree_vamp_100k = loadTree(inputFileName_vamp_100k, inputFilesToClose2d, treeName);

	    std::string histogramName_svFitMEM_vamp_vs_svFitMEM_vamp_100k = Form(
              "histogram_svFitMEM_vamp_vs_svFitMEM_vamp_100k_%s_%s_addLogM%i_numCalls%ik", 
	      sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
	    TH2* histogram_svFitMEM_vamp_vs_svFitMEM_vamp_100k = new TH2D(
              histogramName_svFitMEM_vamp_vs_svFitMEM_vamp_100k.data(), 
	      histogramName_svFitMEM_vamp_vs_svFitMEM_vamp_100k.data(), 100, 0., 2.*mH_i, 100, 0., 2.*mH_i);
	    fillHistogram2dFromTree(histogram_svFitMEM_vamp_vs_svFitMEM_vamp_100k, tree_vamp_100k, branchName_svFitMEM_i, tree_vamp, branchName_svFitMEM_i, mH_i);
	    std::string outputFileName_svFitMEM_vamp_vs_svFitMEM_vamp_100k = Form(
              "plots/makeSVfitMEM_PerformancePlots_svFitMEM_vamp_vs_svFitMEM_vamp_100k_%s_%s_addLogM%i_numCalls%ik.png", 
	      sample->data(), xSection_times_Acc_i->data(), *logM_i, (*numCalls_i) / 1000);
  	    showHistogram2d(
	      800, 800, histogram_svFitMEM_vamp_vs_svFitMEM_vamp_100k, 
	      "m_{#tau#tau}^{VAMP,100k} [GeV]", 1.3, Form("m_{#tau#tau}^{VAMP,%i} [GeV]", (*numCalls_i) / 1000), 1.65, 
	      outputFileName_svFitMEM_vamp_vs_svFitMEM_vamp_100k);
	    delete histogram_svFitMEM_vamp_vs_svFitMEM_vamp_100k;
	  }
 */
	  for ( std::vector<TFile*>::iterator inputFile = inputFilesToClose2d.begin();
		inputFile != inputFilesToClose2d.end(); ++inputFile ) {
	    std::cout << "closing inputFile = " << (*inputFile)->GetName() << std::endl;
	    delete (*inputFile);
	  }
	}
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // make 1d plots of SVfit and SVfitMEM mass distributions

  for ( std::vector<std::string>::const_iterator sample = samplesToAnalyze.begin();
	sample != samplesToAnalyze.end(); ++sample ) {
    
    assert(mH.find(*sample) != mH.end());
    double mH_i = mH[*sample];

    for ( std::vector<std::string>::const_iterator xSection_times_Acc_i = xSection_times_Acc.begin();
	  xSection_times_Acc_i != xSection_times_Acc.end(); ++xSection_times_Acc_i ) {

      std::cout << "processing sample = " << (*sample) << ", xSection_times_Acc = " << (*xSection_times_Acc_i) << " for different values of numCalls..." << std::endl;

      std::vector<TFile*> inputFilesToClose1d;

      assert(branchName_svFitMEM.find(*xSection_times_Acc_i) != branchName_svFitMEM.end());
      std::string branchName_svFitMEM_i = branchName_svFitMEM[*xSection_times_Acc_i];

      std::string inputFileName_vamp_10k = std::string(inputFilePath).append(Form(
        "%s/addLogM5/vamp/maxObjFunctionCalls10k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_10k = loadTree(inputFileName_vamp_10k, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_10k = extractVarFromTree(tree_vamp_10k, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_10k = Form("histogram_svFitMEM_vamp_10k_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_10k = new TH1D(histogramName_svFitMEM_vamp_10k.data(), histogramName_svFitMEM_vamp_10k.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_10k, svFitMEM_vamp_10k);
      normalizeHistogram(histogram_svFitMEM_vamp_10k);

      std::string inputFileName_vamp_20k = std::string(inputFilePath).append(Form(
        "%s/addLogM5/vamp/maxObjFunctionCalls20k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_20k = loadTree(inputFileName_vamp_20k, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_20k = extractVarFromTree(tree_vamp_20k, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_20k = Form("histogram_svFitMEM_vamp_20k_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_20k = new TH1D(histogramName_svFitMEM_vamp_20k.data(), histogramName_svFitMEM_vamp_20k.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_20k, svFitMEM_vamp_20k);
      normalizeHistogram(histogram_svFitMEM_vamp_20k);

      std::string inputFileName_vamp_50k = std::string(inputFilePath).append(Form(
        "%s/addLogM5/vamp/maxObjFunctionCalls50k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_50k = loadTree(inputFileName_vamp_50k, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_50k = extractVarFromTree(tree_vamp_50k, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_50k = Form("histogram_svFitMEM_vamp_50k_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_50k = new TH1D(histogramName_svFitMEM_vamp_50k.data(), histogramName_svFitMEM_vamp_50k.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_50k, svFitMEM_vamp_50k);
      normalizeHistogram(histogram_svFitMEM_vamp_50k);

      std::string inputFileName_vamp_100k = std::string(inputFilePath).append(Form(
        "%s/addLogM5/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_100k = loadTree(inputFileName_vamp_100k, inputFilesToClose1d, treeName);      
      std::vector<double> svFitMEM_vamp_100k = extractVarFromTree(tree_vamp_100k, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_100k = Form("histogram_svFitMEM_vamp_100k_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_100k = new TH1D(histogramName_svFitMEM_vamp_100k.data(), histogramName_svFitMEM_vamp_100k.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_100k, svFitMEM_vamp_100k);
      normalizeHistogram(histogram_svFitMEM_vamp_100k);

      delete tree_vamp_100k;
      tree_vamp_100k = loadTree(inputFileName_vamp_100k, inputFilesToClose1d, treeName);
      std::vector<double> svFit = extractVarFromTree(tree_vamp_100k, "svfitMass", mH_i);
      std::string histogramName_svFit = Form("histogram_svFit_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFit = new TH1D(histogramName_svFit.data(), histogramName_svFit.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFit, svFit);
      normalizeHistogram(histogram_svFit);

      std::string outputFileName_numCalls = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_addLogM5_vamp_%s_%s_numCalls.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_svFitMEM_vamp_10k, "VAMP,10k",
	histogram_svFitMEM_vamp_20k, "VAMP,20k",
	histogram_svFitMEM_vamp_50k, "VAMP,50k",
	histogram_svFitMEM_vamp_100k, "VAMP,100k",
	0, "",
	"m_{#tau#tau} [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm_{#tau#tau} [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_numCalls);

      delete histogram_svFitMEM_vamp_10k;
      delete histogram_svFitMEM_vamp_20k;
      delete histogram_svFitMEM_vamp_50k;
      delete histogram_svFitMEM_vamp_100k;

      std::cout << "processing sample = " << (*sample) << ", xSection_times_Acc = " << (*xSection_times_Acc_i) << " for different values of logM..." << std::endl;

      std::string inputFileName_vamp_logM0 = std::string(inputFilePath).append(Form(
        "%s/addLogM0/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM0 = loadTree(inputFileName_vamp_logM0, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM0 = extractVarFromTree(tree_vamp_logM0, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM0 = Form("histogram_svFitMEM_vamp_logM0_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM0 = new TH1D(histogramName_svFitMEM_vamp_logM0.data(), histogramName_svFitMEM_vamp_logM0.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM0, svFitMEM_vamp_logM0);
      normalizeHistogram(histogram_svFitMEM_vamp_logM0);

      std::string inputFileName_vamp_logM1 = std::string(inputFilePath).append(Form(
        "%s/addLogM1/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM1 = loadTree(inputFileName_vamp_logM1, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM1 = extractVarFromTree(tree_vamp_logM1, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM1 = Form("histogram_svFitMEM_vamp_logM1_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM1 = new TH1D(histogramName_svFitMEM_vamp_logM1.data(), histogramName_svFitMEM_vamp_logM1.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM1, svFitMEM_vamp_logM1);
      normalizeHistogram(histogram_svFitMEM_vamp_logM1);

      std::string inputFileName_vamp_logM2 = std::string(inputFilePath).append(Form(
        "%s/addLogM2/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM2 = loadTree(inputFileName_vamp_logM2, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM2 = extractVarFromTree(tree_vamp_logM2, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM2 = Form("histogram_svFitMEM_vamp_logM2_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM2 = new TH1D(histogramName_svFitMEM_vamp_logM2.data(), histogramName_svFitMEM_vamp_logM2.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM2, svFitMEM_vamp_logM2);
      normalizeHistogram(histogram_svFitMEM_vamp_logM2);

      std::string inputFileName_vamp_logM3 = std::string(inputFilePath).append(Form(
        "%s/addLogM3/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM3 = loadTree(inputFileName_vamp_logM3, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM3 = extractVarFromTree(tree_vamp_logM3, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM3 = Form("histogram_svFitMEM_vamp_logM3_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM3 = new TH1D(histogramName_svFitMEM_vamp_logM3.data(), histogramName_svFitMEM_vamp_logM3.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM3, svFitMEM_vamp_logM3);
      normalizeHistogram(histogram_svFitMEM_vamp_logM3);

      std::string inputFileName_vamp_logM5 = std::string(inputFilePath).append(Form(
        "%s/addLogM5/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM5 = loadTree(inputFileName_vamp_logM5, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM5 = extractVarFromTree(tree_vamp_logM5, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM5 = Form("histogram_svFitMEM_vamp_logM5_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM5 = new TH1D(histogramName_svFitMEM_vamp_logM5.data(), histogramName_svFitMEM_vamp_logM5.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM5, svFitMEM_vamp_logM5);
      normalizeHistogram(histogram_svFitMEM_vamp_logM5);

      std::string inputFileName_vamp_logM6 = std::string(inputFilePath).append(Form(
        "%s/addLogM6/vamp/maxObjFunctionCalls20k/nom/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM6 = loadTree(inputFileName_vamp_logM6, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM6 = extractVarFromTree(tree_vamp_logM6, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM6 = Form("histogram_svFitMEM_vamp_logM6_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM6 = new TH1D(histogramName_svFitMEM_vamp_logM6.data(), histogramName_svFitMEM_vamp_logM6.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM6, svFitMEM_vamp_logM6);
      normalizeHistogram(histogram_svFitMEM_vamp_logM6);

      std::string inputFileName_vamp_logM7 = std::string(inputFilePath).append(Form(
        "%s/addLogM7/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM7 = loadTree(inputFileName_vamp_logM7, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM7 = extractVarFromTree(tree_vamp_logM7, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM7 = Form("histogram_svFitMEM_vamp_logM7_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM7 = new TH1D(histogramName_svFitMEM_vamp_logM7.data(), histogramName_svFitMEM_vamp_logM7.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM7, svFitMEM_vamp_logM7);
      normalizeHistogram(histogram_svFitMEM_vamp_logM7);

      std::string inputFileName_vamp_logM10 = std::string(inputFilePath).append(Form(
        "%s/addLogM10/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM10 = loadTree(inputFileName_vamp_logM10, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM10 = extractVarFromTree(tree_vamp_logM10, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM10 = Form("histogram_svFitMEM_vamp_logM10_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM10 = new TH1D(histogramName_svFitMEM_vamp_logM10.data(), histogramName_svFitMEM_vamp_logM10.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM10, svFitMEM_vamp_logM10);
      normalizeHistogram(histogram_svFitMEM_vamp_logM10);

      std::string inputFileName_vamp_logM15 = std::string(inputFilePath).append(Form(
        "%s/addLogM15/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM15 = loadTree(inputFileName_vamp_logM15, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM15 = extractVarFromTree(tree_vamp_logM15, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM15 = Form("histogram_svFitMEM_vamp_logM15_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM15 = new TH1D(histogramName_svFitMEM_vamp_logM15.data(), histogramName_svFitMEM_vamp_logM15.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM15, svFitMEM_vamp_logM15);
      normalizeHistogram(histogram_svFitMEM_vamp_logM15);

      std::string inputFileName_vamp_logM20 = std::string(inputFilePath).append(Form(
        "%s/addLogM20/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM20 = loadTree(inputFileName_vamp_logM20, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM20 = extractVarFromTree(tree_vamp_logM20, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM20 = Form("histogram_svFitMEM_vamp_logM20_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM20 = new TH1D(histogramName_svFitMEM_vamp_logM20.data(), histogramName_svFitMEM_vamp_logM20.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM20, svFitMEM_vamp_logM20);
      normalizeHistogram(histogram_svFitMEM_vamp_logM20);

      std::string inputFileName_vamp_logM50 = std::string(inputFilePath).append(Form(
        "%s/addLogM50/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM50 = loadTree(inputFileName_vamp_logM50, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM50 = extractVarFromTree(tree_vamp_logM50, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM50 = Form("histogram_svFitMEM_vamp_logM50_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM50 = new TH1D(histogramName_svFitMEM_vamp_logM50.data(), histogramName_svFitMEM_vamp_logM50.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM50, svFitMEM_vamp_logM50);
      normalizeHistogram(histogram_svFitMEM_vamp_logM50);

      std::string inputFileName_vamp_logM100 = std::string(inputFilePath).append(Form(
        "%s/addLogM100/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
        xSection_times_Acc_i->data(), sample->data()));
      TTree* tree_vamp_logM100 = loadTree(inputFileName_vamp_logM100, inputFilesToClose1d, treeName);
      std::vector<double> svFitMEM_vamp_logM100 = extractVarFromTree(tree_vamp_logM100, branchName_svFitMEM_i, mH_i);
      std::string histogramName_svFitMEM_vamp_logM100 = Form("histogram_svFitMEM_vamp_logM100_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_svFitMEM_vamp_logM100 = new TH1D(histogramName_svFitMEM_vamp_logM100.data(), histogramName_svFitMEM_vamp_logM100.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_svFitMEM_vamp_logM100, svFitMEM_vamp_logM100);
      normalizeHistogram(histogram_svFitMEM_vamp_logM100);

      delete tree_vamp_logM100;
      tree_vamp_logM100 = loadTree(inputFileName_vamp_logM100, inputFilesToClose1d, treeName);
      std::vector<double> mTtotal = extractVarFromTree(tree_vamp_logM100, "mTtotal", mH_i);
      std::string histogramName_mTtotal = Form("histogram_mTtotal_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_mTtotal = new TH1D(histogramName_mTtotal.data(), histogramName_mTtotal.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_mTtotal, mTtotal);
      normalizeHistogram(histogram_mTtotal);

      delete tree_vamp_logM100;
      tree_vamp_logM100 = loadTree(inputFileName_vamp_logM100, inputFilesToClose1d, treeName);
      std::vector<double> visMass = extractVarFromTree(tree_vamp_logM100, "visMass", mH_i);
      std::string histogramName_visMass = Form("histogram_visMass_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_visMass = new TH1D(histogramName_visMass.data(), histogramName_visMass.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_visMass, visMass);
      normalizeHistogram(histogram_visMass);

      delete tree_vamp_logM100;
      tree_vamp_logM100 = loadTree(inputFileName_vamp_logM100, inputFilesToClose1d, treeName);
      std::vector<double> caMass = extractCollinearApproxMassFromTree(tree_vamp_logM100, mH_i, 0);
      std::string histogramName_caMass = Form("histogram_caMass_%s_%s", sample->data(), xSection_times_Acc_i->data());
      TH1* histogram_caMass = new TH1D(histogramName_caMass.data(), histogramName_caMass.data(), 60, 0., 3.*mH_i);
      fillHistogram1d(histogram_caMass, caMass);
      normalizeHistogram(histogram_caMass);

      std::string outputFileName_logM_1 = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_vamp_%s_%s_logM_1.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_caMass, "m_{CA}",
	histogram_mTtotal, "m_{T}^{total}",
	histogram_visMass, "m_{vis}",
	histogram_svFitMEM_vamp_logM0, "VAMP,0*logM",
	histogram_svFitMEM_vamp_logM1, "VAMP,1*logM",
	"m [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_logM_1);

      std::string outputFileName_logM_2 = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_vamp_%s_%s_logM_2.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_caMass, "m_{CA}",
	histogram_mTtotal, "m_{T}^{total}",
	histogram_visMass, "m_{vis}",
	histogram_svFitMEM_vamp_logM2, "VAMP,2*logM",
	histogram_svFitMEM_vamp_logM3, "VAMP,3*logM",
	"m [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_logM_2);

      std::string outputFileName_logM_3 = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_vamp_%s_%s_logM_3.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_caMass, "m_{CA}",
	histogram_mTtotal, "m_{T}^{total}",
	histogram_visMass, "m_{vis}",
	histogram_svFitMEM_vamp_logM5, "VAMP,5*logM",
	histogram_svFitMEM_vamp_logM6, "VAMP,6*logM",
	"m [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_logM_3);

      std::string outputFileName_logM_4 = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_vamp_%s_%s_logM_4.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_caMass, "m_{CA}",
	histogram_mTtotal, "m_{T}^{total}",
	histogram_visMass, "m_{vis}",
	histogram_svFitMEM_vamp_logM7, "VAMP,7*logM",
	histogram_svFitMEM_vamp_logM10, "VAMP,10*logM",
	"m [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_logM_4);
      
      std::string outputFileName_logM_5 = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_vamp_%s_%s_logM_5.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_caMass, "m_{CA}",
	histogram_mTtotal, "m_{T}^{total}",
	histogram_visMass, "m_{vis}",
	histogram_svFitMEM_vamp_logM15, "VAMP,15*logM",
	histogram_svFitMEM_vamp_logM20, "VAMP,20*logM",
	"m [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_logM_5);

      std::string outputFileName_logM_6 = Form(
	"plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_vamp_%s_%s_logM_6.png", 
	sample->data(), xSection_times_Acc_i->data());
      showHistograms1d(
        800, 900,
	histogram_svFit, "SVfit",
	histogram_caMass, "m_{CA}",
	histogram_mTtotal, "m_{T}^{total}",
	histogram_visMass, "m_{vis}",
	histogram_svFitMEM_vamp_logM50, "VAMP,50*logM",
	histogram_svFitMEM_vamp_logM100, "VAMP,100*logM",
	"m [GeV]", 1.10,
	true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
	0.24, 0.69,
	outputFileName_logM_6);

      delete histogram_svFit;
      delete histogram_mTtotal;
      delete histogram_visMass;
      delete histogram_caMass;
      delete histogram_svFitMEM_vamp_logM0;
      delete histogram_svFitMEM_vamp_logM1;
      delete histogram_svFitMEM_vamp_logM2;
      delete histogram_svFitMEM_vamp_logM3;
      delete histogram_svFitMEM_vamp_logM5;
      delete histogram_svFitMEM_vamp_logM7;
      delete histogram_svFitMEM_vamp_logM10;
      delete histogram_svFitMEM_vamp_logM15;
      delete histogram_svFitMEM_vamp_logM20;
      delete histogram_svFitMEM_vamp_logM50;
      delete histogram_svFitMEM_vamp_logM100;

      for ( std::vector<TFile*>::iterator inputFile = inputFilesToClose1d.begin();
	    inputFile != inputFilesToClose1d.end(); ++inputFile ) {
	std::cout << "closing inputFile = " << (*inputFile)->GetName() << std::endl;
	delete (*inputFile);
      }
    }
  }
  
  for ( std::vector<std::string>::const_iterator sample = samplesToAnalyze.begin();
	sample != samplesToAnalyze.end(); ++sample ) {
    
    assert(mH.find(*sample) != mH.end());
    double mH_i = mH[*sample];
    
    std::cout << "processing sample = " << (*sample) << " for different values of xSection_times_Acc..." << std::endl;
    
    std::vector<TFile*> inputFilesToClose1d;
    
    assert(branchName_svFitMEM.find("xSection") != branchName_svFitMEM.end());
    std::string branchName_svFitMEM_xSection = branchName_svFitMEM["xSection"];

    std::string inputFileName_xSection_vamp_logM0 = std::string(inputFilePath).append(Form(
      "%s/addLogM0/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
      "xSection", sample->data()));
    TTree* tree_xSection_vamp_logM0 = loadTree(inputFileName_xSection_vamp_logM0, inputFilesToClose1d, treeName);
    std::vector<double> svFitMEM_xSection_vamp_logM0 = extractVarFromTree(tree_xSection_vamp_logM0, branchName_svFitMEM_xSection, mH_i);
    std::string histogramName_svFitMEM_xSection_vamp_logM0 = Form("histogram_svFitMEM_xSection_vamp_logM0_%s_%s", sample->data(), "xSection");
    TH1* histogram_svFitMEM_xSection_vamp_logM0 = new TH1D(histogramName_svFitMEM_xSection_vamp_logM0.data(), histogramName_svFitMEM_xSection_vamp_logM0.data(), 60, 0., 3.*mH_i);
    fillHistogram1d(histogram_svFitMEM_xSection_vamp_logM0, svFitMEM_xSection_vamp_logM0);
    normalizeHistogram(histogram_svFitMEM_xSection_vamp_logM0);

    assert(branchName_svFitMEM.find("none") != branchName_svFitMEM.end());
    std::string branchName_svFitMEM_none = branchName_svFitMEM["none"];

    std::string inputFileName_none_vamp_logM0 = std::string(inputFilePath).append(Form(
      "%s/addLogM0/vamp/maxObjFunctionCalls100k/%s/H2TauTauTreeProducerTauTau_addSVfitMEM_all.root", 
      "none", sample->data()));
    TTree* tree_none_vamp_logM0 = loadTree(inputFileName_none_vamp_logM0, inputFilesToClose1d, treeName);
    std::vector<double> svFitMEM_none_vamp_logM0 = extractVarFromTree(tree_none_vamp_logM0, branchName_svFitMEM_none, mH_i);
    std::string histogramName_svFitMEM_none_vamp_logM0 = Form("histogram_svFitMEM_none_vamp_logM0_%s_%s", sample->data(), "none");
    TH1* histogram_svFitMEM_none_vamp_logM0 = new TH1D(histogramName_svFitMEM_none_vamp_logM0.data(), histogramName_svFitMEM_none_vamp_logM0.data(), 60, 0., 3.*mH_i);
    fillHistogram1d(histogram_svFitMEM_none_vamp_logM0, svFitMEM_none_vamp_logM0);
    normalizeHistogram(histogram_svFitMEM_none_vamp_logM0);

    delete tree_xSection_vamp_logM0;
    tree_xSection_vamp_logM0 = loadTree(inputFileName_xSection_vamp_logM0, inputFilesToClose1d, treeName);
    std::vector<double> svFit = extractVarFromTree(tree_xSection_vamp_logM0, "svfitMass", mH_i);
    std::string histogramName_svFit = Form("histogram_svFit_%s_%s", sample->data(), "xSection");
    TH1* histogram_svFit = new TH1D(histogramName_svFit.data(), histogramName_svFit.data(), 60, 0., 3.*mH_i);
    fillHistogram1d(histogram_svFit, svFit);
    normalizeHistogram(histogram_svFit);

    delete tree_xSection_vamp_logM0;
    tree_xSection_vamp_logM0 = loadTree(inputFileName_xSection_vamp_logM0, inputFilesToClose1d, treeName);
    std::vector<double> mTtotal = extractVarFromTree(tree_xSection_vamp_logM0, "mTtotal", mH_i);
    std::string histogramName_mTtotal = Form("histogram_mTtotal_%s_%s", sample->data(), "xSection");
    TH1* histogram_mTtotal = new TH1D(histogramName_mTtotal.data(), histogramName_mTtotal.data(), 60, 0., 3.*mH_i);
    fillHistogram1d(histogram_mTtotal, mTtotal);
    normalizeHistogram(histogram_mTtotal);
    
    delete tree_xSection_vamp_logM0;
    tree_xSection_vamp_logM0 = loadTree(inputFileName_xSection_vamp_logM0, inputFilesToClose1d, treeName);
    std::vector<double> visMass = extractVarFromTree(tree_xSection_vamp_logM0, "visMass", mH_i);
    std::string histogramName_visMass = Form("histogram_visMass_%s_%s", sample->data(), "xSection");
    TH1* histogram_visMass = new TH1D(histogramName_visMass.data(), histogramName_visMass.data(), 60, 0., 3.*mH_i);
    fillHistogram1d(histogram_visMass, visMass);
    normalizeHistogram(histogram_visMass);
    
    delete tree_xSection_vamp_logM0;
    tree_xSection_vamp_logM0 = loadTree(inputFileName_xSection_vamp_logM0, inputFilesToClose1d, treeName);
    std::vector<double> caMass = extractCollinearApproxMassFromTree(tree_xSection_vamp_logM0, mH_i, 0);
    std::string histogramName_caMass = Form("histogram_caMass_%s_%s", sample->data(), "xSection");
    TH1* histogram_caMass = new TH1D(histogramName_caMass.data(), histogramName_caMass.data(), 60, 0., 3.*mH_i);
    fillHistogram1d(histogram_caMass, caMass);
    normalizeHistogram(histogram_caMass);

    std::string outputFileName_xSection_times_Acc = Form(
      "plots/makeSVfitMEM_PerformancePlots_svFitMEM_numCalls100k_xSection_vs_none_vamp_%s.png", 
      sample->data());
    showHistograms1d(
      800, 900,
      histogram_svFit, "SVfit",
      histogram_caMass, "m_{CA}",
      histogram_mTtotal, "m_{T}^{total}",
      histogram_visMass, "m_{vis}",
      histogram_svFitMEM_xSection_vamp_logM0, "'xSection',0*logM",
      histogram_svFitMEM_none_vamp_logM0, "'none',0*logM",
      "m [GeV]", 1.10,
      true, 1.e-5, 19.9e0, "1/dm [1/GeV]", 1.30,
      0.24, 0.69,
      outputFileName_xSection_times_Acc);
    
    delete histogram_svFit;
    delete histogram_mTtotal;
    delete histogram_visMass;
    delete histogram_caMass;
    delete histogram_svFitMEM_xSection_vamp_logM0;
    delete histogram_svFitMEM_none_vamp_logM0;

    for ( std::vector<TFile*>::iterator inputFile = inputFilesToClose1d.begin();
	  inputFile != inputFilesToClose1d.end(); ++inputFile ) {
      std::cout << "closing inputFile = " << (*inputFile)->GetName() << std::endl;
      delete (*inputFile);
    }
  }
  //-----------------------------------------------------------------------------
}
