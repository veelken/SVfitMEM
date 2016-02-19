

#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>

#include "assert.h"

//TGraphErrors* makeFittedGraph(const TGraphErrors* graph, const TF1* fitFunction1, double fitFunction1_xMax, const TF1* fitFunction2, double fitFunction2_xMax, const TF1* fitFunction3)
//{
//  std::string graphName_fitted = Form("%s_fitted", graph->GetName());
//  TGraphErrors* graph_fitted = (TGraphErrors*)graph->Clone(graphName_fitted.data());
//  graph_fitted->SetMarkerColor(8);
//  graph_fitted->SetLineColor(8);
//  graph_fitted->SetLineStyle(7);
//  graph_fitted->SetLineWidth(1);
//  int numPoints = graph->GetN();
//  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
//    double x, y;
//    graph->GetPoint(idxPoint, x, y);
//    double y_fitted = 0.;
//    if      ( x < fitFunction1_xMax ) y_fitted = fitFunction1->Eval(x);
//    else if ( x < fitFunction2_xMax ) y_fitted = fitFunction2->Eval(x);
//    else                              y_fitted = fitFunction3->Eval(x);
//    graph_fitted->SetPoint(idxPoint, x, y_fitted);
//    graph_fitted->SetPointError(idxPoint, 0.5, 0.);
//  }
//  return graph_fitted;
//}

double findIntersectionX(const TGraphErrors* graph, const TF1* fitFunction1, const TF1* fitFunction2)
{
  double fitFunction1_xMin, fitFunction1_xMax;
  fitFunction1->GetRange(fitFunction1_xMin, fitFunction1_xMax);
  double fitFunction2_xMin, fitFunction2_xMax;
  fitFunction2->GetRange(fitFunction2_xMin, fitFunction2_xMax);
  
  double x0 = -1.;
  double min_dy_rel = -1.;
  bool isFirst = true;

  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    
    if ( x >= fitFunction1_xMin && x <= fitFunction1_xMax &&
	 x >= fitFunction2_xMin && x <= fitFunction2_xMax ) {
      double fitFunction1_y = fitFunction1->Eval(x);
      double fitFunction2_y = fitFunction2->Eval(x);
      double dy_rel = TMath::Abs(fitFunction1_y - fitFunction2_y)/(fitFunction1_y + fitFunction2_y);

      if ( isFirst || dy_rel < min_dy_rel ) {
	x0 = x;
	min_dy_rel = dy_rel;
	isFirst = false;
      }
    }
  }

  return x0;
}

TGraphErrors* makeFittedGraph(const TGraphErrors* graph, const TF1* fitFunction1, const TF1* fitFunction2, const TF1* fitFunction3)
{
  std::cout << "<makeFittedGraph>:" << std::endl;

  TString graphName_fitted = Form("%s_fitted", graph->GetName());
  graphName_fitted = graphName_fitted.ReplaceAll("13TeV", "13TeV_lo");
  TGraphErrors* graph_fitted = (TGraphErrors*)graph->Clone(graphName_fitted.Data());
  graph_fitted->SetMarkerColor(8);
  graph_fitted->SetLineColor(8);
  graph_fitted->SetLineStyle(7);
  graph_fitted->SetLineWidth(1);

  double fitFunction1_xMax = -1.;
  if ( fitFunction1 ) {
    fitFunction1_xMax = fitFunction1->GetXmax();
    if ( fitFunction2 ) fitFunction1_xMax = findIntersectionX(graph, fitFunction1, fitFunction2);
  }
  std::cout << "fitFunction1_xMax = " << fitFunction1_xMax << ":";
  if ( fitFunction1 ) std::cout << " fitFunction1_y = " << fitFunction1->Eval(fitFunction1_xMax);
  if ( fitFunction2 ) std::cout << ", fitFunction2_y = " << fitFunction2->Eval(fitFunction1_xMax);
  std::cout << std::endl;
  double fitFunction2_xMax = -1.;
  if ( fitFunction2 ) {
    fitFunction2_xMax = fitFunction2->GetXmax();
    if ( fitFunction3 ) fitFunction2_xMax = findIntersectionX(graph, fitFunction2, fitFunction3);
  }
  std::cout << "fitFunction2_xMax = " << fitFunction2_xMax << ":";
  if ( fitFunction2 ) std::cout << " fitFunction2_y = " << fitFunction2->Eval(fitFunction2_xMax);
  if ( fitFunction3 ) std::cout << ", fitFunction3_y = " << fitFunction3->Eval(fitFunction2_xMax);
  std::cout << std::endl;

  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    double y_fitted = 0.;
    if ( x <= fitFunction1_xMax ) {
      if ( fitFunction1 ) y_fitted = fitFunction1->Eval(x);
    } else if ( x <= fitFunction2_xMax ) {
      if ( fitFunction2 ) y_fitted = fitFunction2->Eval(x);
    } else {
      if ( fitFunction3 ) y_fitted = fitFunction3->Eval(x);
    }
    //std::cout << "point #" << idxPoint << ": x = " << x << ", y = " << y << " +/- " << yErr << " (fitted = " << y_fitted << ")" << std::endl;
    graph_fitted->SetPoint(idxPoint, x, y_fitted);
    graph_fitted->SetPointError(idxPoint, 0.5, 0.);
  }
  return graph_fitted;
}

void showGraph(double canvasSizeX, double canvasSizeY,
	       TGraph* graph, 
	       TF1* fitFunction1, TF1* fitFunction2, TF1* fitFunction3, 
	       TGraph* graph_fitted, 
	       bool useLogScaleX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
	       bool useLogScaleY, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
	       const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.19);
  canvas->SetBottomMargin(0.19);
  canvas->SetRightMargin(0.05);
  canvas->SetLogx(useLogScaleX);
  canvas->SetLogy(useLogScaleY);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 10, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);
  dummyHistogram->Draw("axis");

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(0.065);
  xAxis->SetLabelSize(0.055);
  xAxis->SetLabelOffset(0.01);
  xAxis->SetTickLength(0.055);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(0.070);
  yAxis->SetLabelSize(0.055);
  yAxis->SetLabelOffset(0.01);
  yAxis->SetTickLength(0.055);
  yAxis->SetNdivisions(505);

  graph->Draw("Psame");

  if ( fitFunction1 ) fitFunction1->Draw("same");
  if ( fitFunction2 ) fitFunction2->Draw("same");
  if ( fitFunction3 ) fitFunction3->Draw("same");

  graph_fitted->Draw("Lsame");

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScaleY ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete dummyHistogram;
  delete canvas;  
}

void fitSVfitMEM_xSection_lo()
{
  gROOT->SetBatch(true);

  std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/SVfitMEM_with_vamp/CMSSW_7_4_6/src/TauAnalysis/SVfitMEM/data/";
  std::string inputFileName = "svFitMEM_xSection_and_AccCorr_13TeV_lo.root";
  std::string inputFileName_full = Form("%s%s", inputFilePath.data(), inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }
  
  std::vector<std::string> channels;
  channels.push_back("hadhad");
  channels.push_back("lephad");
  //channels.push_back("hadlep");
  channels.push_back("leplep");

  std::map<std::string, std::string> graphNames_xSection_woAcc; // key = channel
  graphNames_xSection_woAcc["hadhad"] = "graph_Xsection_woAcc_13TeV_hadhad_vamp";
  graphNames_xSection_woAcc["lephad"] = "graph_Xsection_woAcc_13TeV_lephad_vamp";
  graphNames_xSection_woAcc["hadlep"] = "graph_Xsection_woAcc_13TeV_hadlep_vamp";
  graphNames_xSection_woAcc["leplep"] = "graph_Xsection_woAcc_13TeV_leplep_vamp";

  std::map<std::string, std::string> graphNames_Acc; // key = channel
  graphNames_Acc["hadhad"] = "graph_Acc_13TeV_hadhad_vamp";
  graphNames_Acc["lephad"] = "graph_Acc_13TeV_lephad_vamp";
  graphNames_Acc["hadlep"] = "graph_Acc_13TeV_hadlep_vamp";
  graphNames_Acc["leplep"] = "graph_Acc_13TeV_leplep_vamp";
  
  double fitFunction1_xSection_woAcc_xMin = 5.e+1;
  double fitFunction1_xSection_woAcc_xMax = 3.e+2;
  //double fitFunction2_xSection_woAcc_xMin = 3.e+2;
  double fitFunction2_xSection_woAcc_xMin = 2.8e+2;
  double fitFunction2_xSection_woAcc_xMax = 4.6e+2;
  //double fitFunction3_xSection_woAcc_xMin = 4.5e+2;
  double fitFunction3_xSection_woAcc_xMin = 4.4e+2;
  double fitFunction3_xSection_woAcc_xMax = 5.e+3;

  std::vector<TGraph*> graphs_xSection_woAcc_fitted;
/*
  for ( std::vector<std::string>::const_iterator channel = channels.begin();
	channel != channels.end(); ++channel ) {
    if ( graphNames_xSection_woAcc.find(*channel) == graphNames_xSection_woAcc.end() ) {
      std::cerr << "No graph defined for channel = " << (*channel) << " !!" << std::endl;
      assert(0);
    }
    std::string graphName_xSection_woAcc = graphNames_xSection_woAcc[*channel];
    assert(graphName_xSection_woAcc != "");
    TGraphErrors* graph_xSection_woAcc = dynamic_cast<TGraphErrors*>(inputFile->Get(graphName_xSection_woAcc.data()));
    if ( !graph_xSection_woAcc ) {
      std::cerr << "Failed to load graph = '" << graphName_xSection_woAcc << "' from file = '" << inputFileName.data() << "' !!" << std::endl;
      inputFile->ls();
      assert(0);
    }

    std::string graphName_xSection_woAcc_cloned = Form("%s_cloned", graph_xSection_woAcc->GetName());
    TGraphErrors* graph_xSection_woAcc_cloned = (TGraphErrors*)graph_xSection_woAcc->Clone(graphName_xSection_woAcc_cloned.data());

    std::string fitFunction1_xSection_woAcc_formula = "[0]*TMath::Power(x, [1]) + [2]*TMath::Power(x, [3]) + [4]*TMath::Exp([5]*(TMath::Power(x, [6]) + TMath::Power(x, [7])))";
    TF1* fitFunction1_xSection_woAcc = new TF1("fitFunction1_xSection_woAcc", fitFunction1_xSection_woAcc_formula.data(), fitFunction1_xSection_woAcc_xMin, fitFunction1_xSection_woAcc_xMax);
    fitFunction1_xSection_woAcc->SetLineColor(2);
    fitFunction1_xSection_woAcc->SetLineWidth(1);

    std::string fitFunction2_xSection_woAcc_formula = "[1] + (x - [0])*([2] + (x - [0])*([3] + (x - [0])*([4] + (x - [0])*([5] + (x - [0])*([6] + (x - [0])*([7]";
    fitFunction2_xSection_woAcc_formula.append(" + (x - [0])*([8] + (x - [0])*([9] + (x - [0])*([10] + (x - [0])*([11] + (x - [0])*([12] + (x - [0])*([13] + (x - [0])*[14]))))))))))))");
    TF1* fitFunction2_xSection_woAcc = new TF1("fitFunction2_xSection_woAcc", fitFunction2_xSection_woAcc_formula.data(), fitFunction2_xSection_woAcc_xMin, fitFunction2_xSection_woAcc_xMax);
    fitFunction2_xSection_woAcc->SetLineColor(6);
    fitFunction2_xSection_woAcc->SetLineWidth(1);

    std::string fitFunction3_xSection_woAcc_formula = "[0]/(x*x*x) + [1]/(x*x) + [2]/x + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]*x*x*x*x + [8]*x*x*x*x*x + [9]*x*x*x*x*x*x";
    TF1* fitFunction3_xSection_woAcc = new TF1("fitFunction3_xSection_woAcc", fitFunction3_xSection_woAcc_formula.data(), fitFunction3_xSection_woAcc_xMin, fitFunction3_xSection_woAcc_xMax);
    fitFunction3_xSection_woAcc->SetLineColor(4);
    fitFunction3_xSection_woAcc->SetLineWidth(1);

    for ( int idxPass = 0; idxPass <= 1; ++idxPass ) { // CV: add second pass to reject outliers (points on graph with large distance to fit)
      int numPoints = graph_xSection_woAcc->GetN();
      for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
	double x, y;
	graph_xSection_woAcc->GetPoint(idxPoint, x, y);
	if ( x >= 5.e+2 ) {
	  double xErr = graph_xSection_woAcc->GetErrorX(idxPoint);
	  double yErr = graph_xSection_woAcc->GetErrorY(idxPoint);
	  if ( idxPass == 0 ) {
	    yErr = TMath::Max(yErr, 1.e-4);
	    yErr = TMath::Max(yErr, 0.1*y);
	  } else {
	    double y_fitted = 0.;
	    if      ( x <= fitFunction1_xSection_woAcc_xMax ) y_fitted = fitFunction1_xSection_woAcc->Eval(x);
	    else if ( x <= fitFunction2_xSection_woAcc_xMax ) y_fitted = fitFunction2_xSection_woAcc->Eval(x);
	    else                                              y_fitted = fitFunction3_xSection_woAcc->Eval(x);
	    yErr = TMath::Max(yErr, 0.1*TMath::Abs(y - y_fitted));
	    yErr = TMath::Max(yErr, 0.1*y_fitted);	    
	  }
	  graph_xSection_woAcc_cloned->SetPointError(idxPoint, xErr, yErr);
	}
      }
  
      if ( idxPass == 0 ) {
	fitFunction1_xSection_woAcc->SetParameter(0,  1.e+2);
	fitFunction1_xSection_woAcc->SetParameter(1, -1.);
	fitFunction1_xSection_woAcc->SetParameter(2,  1.e+1);
	fitFunction1_xSection_woAcc->SetParameter(3, -1.e-2);
	fitFunction1_xSection_woAcc->SetParameter(4,  1.e-4);
	fitFunction1_xSection_woAcc->SetParameter(5, -1.e-3);
	fitFunction1_xSection_woAcc->SetParameter(6,  1.);
	fitFunction1_xSection_woAcc->SetParameter(7,  2.);
      }
      graph_xSection_woAcc_cloned->Fit(fitFunction1_xSection_woAcc, "", "", fitFunction1_xSection_woAcc_xMin, fitFunction1_xSection_woAcc_xMax);
      
      if ( idxPass == 0 ) {
	fitFunction3_xSection_woAcc->SetParameter(0,  1.);
	fitFunction3_xSection_woAcc->SetParameter(1,  1.e-3);
	fitFunction3_xSection_woAcc->SetParameter(2,  1.e-6);
	fitFunction3_xSection_woAcc->SetParameter(3,  1.e-9);
	fitFunction3_xSection_woAcc->SetParameter(4,  1.e-12);
	fitFunction3_xSection_woAcc->SetParameter(5,  1.e-15);
	fitFunction3_xSection_woAcc->SetParameter(6,  1.e-18);
	fitFunction3_xSection_woAcc->SetParameter(7,  1.e-21);
	fitFunction3_xSection_woAcc->SetParameter(8,  1.e-24);
	fitFunction3_xSection_woAcc->SetParameter(9,  1.e-27);
      }
      graph_xSection_woAcc_cloned->Fit(fitFunction3_xSection_woAcc, "", "", fitFunction3_xSection_woAcc_xMin, fitFunction3_xSection_woAcc_xMax);

      if ( idxPass == 0 ) {
	fitFunction2_xSection_woAcc->FixParameter(0,  3.5e+2);
	fitFunction2_xSection_woAcc->SetParameter(1,  1.e+6);
	fitFunction2_xSection_woAcc->SetParameter(2,  1.e+3);
	fitFunction2_xSection_woAcc->SetParameter(3,  1.);
	fitFunction2_xSection_woAcc->SetParameter(4,  1.e-3);
	fitFunction2_xSection_woAcc->SetParameter(5,  1.e-6);
	fitFunction2_xSection_woAcc->SetParameter(6,  1.e-9);
	fitFunction2_xSection_woAcc->SetParameter(7,  1.e-12);
	fitFunction2_xSection_woAcc->SetParameter(8,  1.e-15);
	fitFunction2_xSection_woAcc->SetParameter(9,  1.e-18);
	fitFunction2_xSection_woAcc->SetParameter(10, 1.e-21);
	fitFunction2_xSection_woAcc->SetParameter(11, 1.e-24);
	fitFunction2_xSection_woAcc->SetParameter(12, 1.e-27);
	fitFunction2_xSection_woAcc->SetParameter(13, 1.e-30);
	fitFunction2_xSection_woAcc->SetParameter(14, 1.e-33);
      }
      graph_xSection_woAcc_cloned->Fit(fitFunction2_xSection_woAcc, "", "", fitFunction2_xSection_woAcc_xMin, fitFunction2_xSection_woAcc_xMax);
    }

    TGraphErrors* graph_xSection_woAcc_fitted = makeFittedGraph(graph_xSection_woAcc, fitFunction1_xSection_woAcc, fitFunction2_xSection_woAcc, fitFunction3_xSection_woAcc);

    std::string outputFileName_xSection_woAcc = Form("./plots/fitSVfitMEM_lo_xSection_woAcc_%s.pdf", channel->data());
    showGraph(800, 650,
	      graph_xSection_woAcc, 
	      fitFunction1_xSection_woAcc, fitFunction2_xSection_woAcc, fitFunction3_xSection_woAcc,
	      graph_xSection_woAcc_fitted, 
	      false, 50., 5000., "m_{H} [GeV]", 1.30,
	      true, 5.e-7, 1.e+2, "#sigma(gg #rightarrow H) [pb]", 1.35,
	      outputFileName_xSection_woAcc);
    
    graphs_xSection_woAcc_fitted.push_back(graph_xSection_woAcc_fitted);
  }
 */
  double fitFunction1_Acc_xMin = 5.e+1;
  double fitFunction1_Acc_xMax = 1.15e+3;
  double fitFunction2_Acc_xMin = 0.85e+3;
  double fitFunction2_Acc_xMax = 5.e+3;
  double fitFunction3_Acc_xMin = 2.5e+3;
  double fitFunction3_Acc_xMax = 5.e+3;

  std::vector<TGraph*> graphs_Acc_fitted;
  
  for ( std::vector<std::string>::const_iterator channel = channels.begin();
	channel != channels.end(); ++channel ) {
    if ( graphNames_Acc.find(*channel) == graphNames_Acc.end() ) {
      std::cerr << "No graph defined for channel = " << (*channel) << " !!" << std::endl;
      assert(0);
    }
    std::string graphName_Acc = graphNames_Acc[*channel];
    assert(graphName_Acc != "");
    TGraphErrors* graph_Acc = dynamic_cast<TGraphErrors*>(inputFile->Get(graphName_Acc.data()));
    if ( !graph_Acc ) {
      std::cerr << "Failed to load graph = '" << graphName_Acc << "' from file = '" << inputFileName.data() << "' !!" << std::endl;
      inputFile->ls();
      assert(0);
    }

    std::string graphName_Acc_cloned = Form("%s_cloned", graph_Acc->GetName());
    TGraphErrors* graph_Acc_cloned = (TGraphErrors*)graph_Acc->Clone(graphName_Acc_cloned.data());

    std::string fitFunction1_Acc_formula = "TMath::Erf((x - [0])/[1])*([2] + (x - [3])*([4] + (x - [3])*([5] + (x - [3])*([6] + (x - [3])*[7]))))";
    //std::string fitFunction1_Acc_formula = "TMath::Erf((x - [0])/[1])*[2]";
    TF1* fitFunction1_Acc = new TF1("fitFunction1_Acc", fitFunction1_Acc_formula.data(), fitFunction1_Acc_xMin, fitFunction1_Acc_xMax);
    fitFunction1_Acc->SetLineColor(2);
    fitFunction1_Acc->SetLineWidth(1);
    
    std::string fitFunction2_Acc_formula = "TMath::Erf(([0] - x)/[1])*([2] + ([3] - x)*([4] + ([3] - x)*([5] + ([3] - x)*([6] + ([3] - x)*[7])))) + [8]";
    TF1* fitFunction2_Acc = new TF1("fitFunction1_Acc", fitFunction2_Acc_formula.data(), fitFunction2_Acc_xMin, fitFunction2_Acc_xMax);
    fitFunction2_Acc->SetLineColor(6);
    fitFunction2_Acc->SetLineWidth(1);

    std::string fitFunction3_Acc_formula = "[0]";
    TF1* fitFunction3_Acc = new TF1("fitFunction3_Acc", fitFunction3_Acc_formula.data(), fitFunction3_Acc_xMin, fitFunction3_Acc_xMax);
    fitFunction3_Acc->SetLineColor(4);
    fitFunction3_Acc->SetLineWidth(1);

    for ( int idxPass = 0; idxPass <= 1; ++idxPass ) { // CV: add second pass to reject outliers (points on graph with large distance to fit)
      int numPoints = graph_Acc->GetN();
      for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
	double x, y;
	graph_Acc->GetPoint(idxPoint, x, y);
	double xErr = graph_Acc->GetErrorX(idxPoint);
	double yErr = graph_Acc->GetErrorY(idxPoint);
	if ( yErr < 1.e-2 || TMath::IsNaN(yErr) ) yErr = 1.e-2;
	if ( idxPass == 0 ) {
	  yErr = TMath::Sqrt(yErr*yErr + 2.5e-3);
	} else {
	  if ( x >= 1.e+3 ) {
	    double y_fitted = 0.;
	    if      ( x <= fitFunction1_Acc_xMax ) y_fitted = fitFunction1_Acc->Eval(x);
	    else if ( x <= fitFunction2_Acc_xMax ) y_fitted = fitFunction2_Acc->Eval(x);
	    else                                   y_fitted = fitFunction3_Acc->Eval(x);
	    yErr = TMath::Max(yErr, 0.1*TMath::Abs(y - y_fitted));
	    yErr = TMath::Max(yErr, 0.1*y_fitted);	    
	  }
	}
	graph_Acc_cloned->SetPointError(idxPoint, xErr, yErr);
      }

      if ( idxPass == 0 ) {
	double x0 = 50.;
	if ( (*channel) == "hadhad" ) x0 = 90.;
	fitFunction1_Acc->SetParameter(0,  x0);
	fitFunction1_Acc->SetParameter(1,  7.e+2);
	fitFunction1_Acc->SetParameter(2,  0.7);
	fitFunction1_Acc->FixParameter(3,  x0);
	fitFunction1_Acc->SetParameter(4,  0.e-3);
	fitFunction1_Acc->SetParameter(5,  0.e-6);
	fitFunction1_Acc->SetParameter(6,  0.e-9);
	fitFunction1_Acc->SetParameter(7,  0.e-12);
      }
      graph_Acc_cloned->Fit(fitFunction1_Acc, "", "", fitFunction1_Acc_xMin, fitFunction1_Acc_xMax);

      if ( idxPass == 0 ) {
	double x0 = 3.e+3;
	fitFunction2_Acc->SetParameter(0,  x0);
	fitFunction2_Acc->SetParameter(1,  1.5e+3);
	fitFunction2_Acc->SetParameter(2,  0.5);
	fitFunction2_Acc->FixParameter(3,  x0);
	fitFunction2_Acc->SetParameter(4,  0.e-3);
	fitFunction2_Acc->SetParameter(5,  0.e-6);
	fitFunction2_Acc->SetParameter(6,  0.e-9);
	fitFunction2_Acc->SetParameter(7,  0.e-12);
	fitFunction2_Acc->SetParameter(8,  0.2);
      }
      graph_Acc_cloned->Fit(fitFunction2_Acc, "", "", fitFunction2_Acc_xMin, fitFunction2_Acc_xMax);

      if ( idxPass == 0 ) {
	fitFunction3_Acc->SetParameter(0,  0.2);
      }
      graph_Acc_cloned->Fit(fitFunction3_Acc, "", "", fitFunction3_Acc_xMin, fitFunction3_Acc_xMax);
    }

    TGraphErrors* graph_Acc_fitted = makeFittedGraph(graph_Acc, fitFunction1_Acc, fitFunction2_Acc, fitFunction3_Acc);

    std::string outputFileName_Acc = Form("./plots/fitSVfitMEM_lo_Acc_%s.pdf", channel->data());
    showGraph(800, 650,
	      graph_Acc, 
	      fitFunction1_Acc, fitFunction2_Acc, fitFunction3_Acc, 
	      graph_Acc_fitted, 
	      false, 50., 5000., "m_{H} [GeV]", 1.30,
	      false, 0., 0.8, "Acceptance", 1.35,
	      outputFileName_Acc);
    
    graphs_Acc_fitted.push_back(graph_Acc_fitted);
  }

  std::string outputFileName = TString(inputFileName.data()).ReplaceAll(".root", "_fitted.root").Data();
  TFile* outputFile = new TFile(outputFileName.data(), "RECREATE");
  for ( std::vector<TGraph*>::iterator graph = graphs_xSection_woAcc_fitted.begin();
	graph != graphs_xSection_woAcc_fitted.end(); ++graph ) {
    (*graph)->Write();
  }
  for ( std::vector<TGraph*>::iterator graph = graphs_Acc_fitted.begin();
	graph != graphs_Acc_fitted.end(); ++graph ) {
    (*graph)->Write();
  }
  delete outputFile;

  delete inputFile;
}
