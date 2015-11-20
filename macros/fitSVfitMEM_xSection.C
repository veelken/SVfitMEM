
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

  std::string graphName_fitted = Form("%s_fitted", graph->GetName());
  TGraphErrors* graph_fitted = (TGraphErrors*)graph->Clone(graphName_fitted.data());
  graph_fitted->SetMarkerColor(8);
  graph_fitted->SetLineColor(8);
  graph_fitted->SetLineStyle(7);
  graph_fitted->SetLineWidth(1);

  double fitFunction1_xMax = findIntersectionX(graph, fitFunction1, fitFunction2);
  std::cout << "fitFunction1_xMax = " << fitFunction1_xMax << ": fitFunction1_y = " << fitFunction1->Eval(fitFunction1_xMax) << ", fitFunction2_y = " << fitFunction2->Eval(fitFunction1_xMax) << std::endl;
  double fitFunction2_xMax = findIntersectionX(graph, fitFunction2, fitFunction3);
  std::cout << "fitFunction2_xMax = " << fitFunction2_xMax << ": fitFunction2_y = " << fitFunction2->Eval(fitFunction2_xMax) << ", fitFunction3_y = " << fitFunction3->Eval(fitFunction2_xMax) << std::endl;

  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    double y_fitted = 0.;
    if      ( x <= fitFunction1_xMax ) y_fitted = fitFunction1->Eval(x);
    else if ( x <= fitFunction2_xMax ) y_fitted = fitFunction2->Eval(x);
    else                               y_fitted = fitFunction3->Eval(x);
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

  fitFunction1->Draw("same");
  fitFunction2->Draw("same");
  fitFunction3->Draw("same");

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

class fitFunction2_fcn 
{
 public:
  Double_t Evaluate(Double_t* x, Double_t* par) 
  {
    double x0 = par[0];

    double xMin = par[1];
    double y_xMin = par[2];
    double xMax = par[3];
    double y_xMax = par[4];
    
    double c2  = par[5];
    double c3  = par[6];
    double c4  = par[7];
    double c5  = par[8];
    double c6  = par[9];
    double c7  = par[10];
    double c8  = par[11];
    double c9  = par[12];
    double c10 = par[13];
    double c11 = par[14];
    double c12 = par[15];
    double c13 = par[16];
    
    double dxMin = xMin - x0;
    double f_xMin = dxMin*dxMin*(c2 + dxMin*(c3 + dxMin*(c4 + dxMin*(c5 + dxMin*(c6 + dxMin*(c7 + dxMin*(c8 + dxMin*(c9 + dxMin*(c10 + dxMin*(c11 + dxMin*(c12 + dxMin*c13)))))))))));
    double dxMax = xMax - x0;
    double f_xMax = dxMax*dxMax*(c2 + dxMax*(c3 + dxMax*(c4 + dxMax*(c5 + dxMax*(c6 + dxMax*(c7 + dxMax*(c8 + dxMax*(c9 + dxMax*(c10 + dxMax*(c11 + dxMax*(c12 + dxMax*c13)))))))))));
    double c1 = ((y_xMax - f_xMax) - (y_xMin - f_xMin))/(xMax - xMin);
    double c0 = y_xMin - c1*xMin - f_xMin;
    
    double dx = x[0] - x0;
    double f_x = c0 + dx*(c1 + dx*(c2 + dx*(c3 + dx*(c4 + dx*(c5 + dx*(c6 + dx*(c7 + dx*(c8 + dx*(c9 + dx*(c10 + dx*(c11 + dx*(c12 + dx*c13))))))))))));
    
    return f_x;
  }
};

void fitSVfitMEM_xSection()
{
  gROOT->SetBatch(true);

  std::string inputFilePath = "/afs/cern.ch/user/v/veelken/scratch0/SVfitMEM_with_vamp/CMSSW_7_4_6/src/TauAnalysis/SVfitMEM/data/";
  std::string inputFileName = "svFitMEM_xSection_and_AccCorr_13TeV.root";
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

  std::map<std::string, std::string> graphNames; // key = channel
  graphNames["hadhad"] = "graph_Xsection_woAcc_13TeV_hadhad_vamp";
  graphNames["lephad"] = "graph_Xsection_woAcc_13TeV_lephad_vamp";
  graphNames["hadlep"] = "graph_Xsection_woAcc_13TeV_hadlep_vamp";
  graphNames["leplep"] = "graph_Xsection_woAcc_13TeV_leplep_vamp";
  
  double fitFunction1_xMin = 5.e+1;
  double fitFunction1_xMax = 3.e+2;
  //double fitFunction2_xMin = 3.e+2;
  double fitFunction2_xMin = 2.8e+2;
  double fitFunction2_xMax = 4.6e+2;
  //double fitFunction3_xMin = 4.5e+2;
  double fitFunction3_xMin = 4.4e+2;
  double fitFunction3_xMax = 5.e+3;

  std::vector<TGraph*> graphs_fitted;

  for ( std::vector<std::string>::const_iterator channel = channels.begin();
	channel != channels.end(); ++channel ) {
    if ( graphNames.find(*channel) == graphNames.end() ) {
      std::cerr << "No graph defined for channel = " << (*channel) << " !!" << std::endl;
      assert(0);
    }
    std::string graphName = graphNames[*channel];
    assert(graphName != "");
    TGraphErrors* graph = dynamic_cast<TGraphErrors*>(inputFile->Get(graphName.data()));
    if ( !graph ) {
      std::cerr << "Failed to load graph = '" << graphName << "' from file = '" << inputFileName.data() << "' !!" << std::endl;
      inputFile->ls();
      assert(0);
    }

    std::string graphName_cloned = Form("%s_cloned", graph->GetName());
    TGraphErrors* graph_cloned = (TGraphErrors*)graph->Clone(graphName_cloned.data());

    std::string fitFunction1_formula = "[0]*TMath::Power(x, [1]) + [2]*TMath::Power(x, [3]) + [4]*TMath::Exp([5]*(TMath::Power(x, [6]) + TMath::Power(x, [7])))";
    TF1* fitFunction1 = new TF1("fitFunction1", fitFunction1_formula.data(), fitFunction1_xMin, fitFunction1_xMax);
    fitFunction1->SetLineColor(2);
    fitFunction1->SetLineWidth(1);

    std::string fitFunction2_formula = "[1] + (x - [0])*([2] + (x - [0])*([3] + (x - [0])*([4] + (x - [0])*([5] + (x - [0])*([6] + (x - [0])*([7]";
    fitFunction2_formula.append(" + (x - [0])*([8] + (x - [0])*([9] + (x - [0])*([10] + (x - [0])*([11] + (x - [0])*([12] + (x - [0])*([13] + (x - [0])*[14]))))))))))))");
    TF1* fitFunction2 = new TF1("fitFunction2", fitFunction2_formula.data(), fitFunction2_xMin, fitFunction2_xMax);
    //fitFunction2_fcn* fitFunction2_ptr = new fitFunction2_fcn();
    //TF1* fitFunction2 = new TF1("fitFunction2", fitFunction2_ptr, &fitFunction2_fcn::Evaluate, fitFunction2_xMin, fitFunction2_xMax, 17);
    fitFunction2->SetLineColor(6);
    fitFunction2->SetLineWidth(1);

    std::string fitFunction3_formula = "[0]/(x*x*x) + [1]/(x*x) + [2]/x + [3] + [4]*x + [5]*x*x + [6]*x*x*x + [7]*x*x*x*x + [8]*x*x*x*x*x + [9]*x*x*x*x*x*x";
    TF1* fitFunction3 = new TF1("fitFunction3", fitFunction3_formula.data(), fitFunction3_xMin, fitFunction3_xMax);
    fitFunction3->SetLineColor(4);
    fitFunction3->SetLineWidth(1);

    for ( int idxPass = 0; idxPass <= 1; ++idxPass ) { // CV: add second pass to reject outliers (points on graph with large distance to fit)
      int numPoints = graph->GetN();
      for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
	double x, y;
	graph->GetPoint(idxPoint, x, y);
	if ( x >= 5.e+2 ) {
	  double xErr = graph->GetErrorX(idxPoint);
	  double yErr = graph->GetErrorY(idxPoint);
	  if ( idxPass == 0 ) {
	    yErr = TMath::Max(yErr, 1.e-4);
	    yErr = TMath::Max(yErr, 0.1*y);
	  } else {
	    double y_fitted = 0.;
	    if      ( x <= fitFunction1_xMax ) y_fitted = fitFunction1->Eval(x);
	    else if ( x <= fitFunction2_xMax ) y_fitted = fitFunction2->Eval(x);
	    else                               y_fitted = fitFunction3->Eval(x);
	    yErr = TMath::Max(yErr, 0.1*TMath::Abs(y - y_fitted));
	    yErr = TMath::Max(yErr, 0.1*y_fitted);	    
	  }
	  graph_cloned->SetPointError(idxPoint, xErr, yErr);
	}
      }
  
      if ( idxPass == 0 ) {
	fitFunction1->SetParameter(0,  1.e+2);
	fitFunction1->SetParameter(1, -1.);
	fitFunction1->SetParameter(2,  1.e+1);
	fitFunction1->SetParameter(3, -1.e-2);
	fitFunction1->SetParameter(4,  1.e-4);
	fitFunction1->SetParameter(5, -1.e-3);
	fitFunction1->SetParameter(6,  1.);
	fitFunction1->SetParameter(7,  2.);
      }
      graph_cloned->Fit(fitFunction1, "", "", fitFunction1_xMin, fitFunction1_xMax);
      
      if ( idxPass == 0 ) {
	fitFunction3->SetParameter(0,  1.);
	fitFunction3->SetParameter(1,  1.e-3);
	fitFunction3->SetParameter(2,  1.e-6);
	fitFunction3->SetParameter(3,  1.e-9);
	fitFunction3->SetParameter(4,  1.e-12);
	fitFunction3->SetParameter(5,  1.e-15);
	fitFunction3->SetParameter(6,  1.e-18);
	fitFunction3->SetParameter(7,  1.e-21);
	fitFunction3->SetParameter(8,  1.e-24);
	fitFunction3->SetParameter(9,  1.e-27);
      }
      graph_cloned->Fit(fitFunction3, "", "", fitFunction3_xMin, fitFunction3_xMax);

      if ( idxPass == 0 ) {
	fitFunction2->FixParameter(0,  3.5e+2);
	fitFunction2->SetParameter(1,  1.e+6);
	fitFunction2->SetParameter(2,  1.e+3);
	fitFunction2->SetParameter(3,  1.);
	fitFunction2->SetParameter(4,  1.e-3);
	fitFunction2->SetParameter(5,  1.e-6);
	fitFunction2->SetParameter(6,  1.e-9);
	fitFunction2->SetParameter(7,  1.e-12);
	fitFunction2->SetParameter(8,  1.e-15);
	fitFunction2->SetParameter(9,  1.e-18);
	fitFunction2->SetParameter(10, 1.e-21);
	fitFunction2->SetParameter(11, 1.e-24);
	fitFunction2->SetParameter(12, 1.e-27);
	fitFunction2->SetParameter(13, 1.e-30);
	fitFunction2->SetParameter(14, 1.e-33);
	//fitFunction2->FixParameter(0,  3.5e+2);
	//fitFunction2->FixParameter(1,  fitFunction2_xMin);
	//fitFunction2->FixParameter(2,  fitFunction1->Eval(fitFunction2_xMin));
	//fitFunction2->FixParameter(3,  fitFunction2_xMax);
	//fitFunction2->FixParameter(4,  fitFunction3->Eval(fitFunction2_xMax));
	//fitFunction2->SetParameter(5,  1.);
	//fitFunction2->SetParameter(6,  1.e-3);
	//fitFunction2->SetParameter(7,  1.e-6);
	//fitFunction2->SetParameter(8,  1.e-9);
	//fitFunction2->SetParameter(9,  1.e-12);
	//fitFunction2->SetParameter(10, 1.e-15);
	//fitFunction2->SetParameter(11, 1.e-18);
	//fitFunction2->SetParameter(12, 1.e-21);
	//fitFunction2->SetParameter(13, 1.e-24);
	//fitFunction2->SetParameter(14, 1.e-27);
	//fitFunction2->SetParameter(15, 1.e-30);
	//fitFunction2->SetParameter(16, 1.e-33);
      }
      graph_cloned->Fit(fitFunction2, "", "", fitFunction2_xMin, fitFunction2_xMax);
    }

    //TGraphErrors* graph_fitted = makeFittedGraph(graph, fitFunction1, fitFunction1_xMax, fitFunction2, fitFunction2_xMax, fitFunction3);
    TGraphErrors* graph_fitted = makeFittedGraph(graph, fitFunction1, fitFunction2, fitFunction3);

    std::string outputFileName = Form("./plots/fitSVfitMEM_xSection_%s.pdf", channel->data());
    showGraph(800, 650,
	      graph, 
	      fitFunction1, fitFunction2, fitFunction3,
	      graph_fitted, 
	      false, 50., 5000., "m_{H} [GeV]", 1.30,
	      true, 5.e-7, 1.e+2, "#sigma(gg #rightarrow H) [pb]", 1.35,
	      outputFileName);
    
    graphs_fitted.push_back(graph_fitted);

    //delete fitFunction2_ptr;
  }

  std::string outputFileName = TString(inputFileName.data()).ReplaceAll(".root", "_fitted.root").Data();
  TFile* outputFile = new TFile(outputFileName.data(), "RECREATE");
  for ( std::vector<TGraph*>::iterator graph = graphs_fitted.begin();
	graph != graphs_fitted.end(); ++graph ) {
    (*graph)->Write();
  }
  delete outputFile;

  delete inputFile;
}
