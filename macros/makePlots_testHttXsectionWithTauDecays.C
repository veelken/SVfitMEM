
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

TGraphErrors* getGraph(TFile* inputFile, const std::string& graphName)
{
  TGraphErrors* graph = dynamic_cast<TGraphErrors*>(inputFile->Get(graphName.data()));
  if ( !graph ) {
    std::cerr << "<getGraph>: Failed to find graph = " << graphName << " in file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  return graph;
}

void fillMap(const std::string& inputFileName, const std::string& graphName, std::map<int, double>& xSection, std::map<int, double>& xSectionErr)
{
  TFile* inputFile = new TFile(inputFileName.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName << " !!" << std::endl;
    assert(0);
  }

  TGraphErrors* graph = getGraph(inputFile, graphName);

  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    double yErr = graph->GetErrorY(idxPoint);
    int x_int = TMath::Nint(x);
    xSection[x_int] = y;
    xSectionErr[x_int] = yErr;
  }
 
  delete inputFile;
}

//-------------------------------------------------------------------------------
void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& drawOption1, const std::string& legendEntry1, const std::string& legendOption1, 
		TGraph* graph2, const std::string& drawOption2, const std::string& legendEntry2, const std::string& legendOption2, 
		TGraph* graph3, const std::string& drawOption3, const std::string& legendEntry3, const std::string& legendOption3, 
		TGraph* graph4, const std::string& drawOption4, const std::string& legendEntry4, const std::string& legendOption4, 
		TGraph* graph5, const std::string& drawOption5, const std::string& legendEntry5, const std::string& legendOption5, 
		TGraph* graph6, const std::string& drawOption6, const std::string& legendEntry6, const std::string& legendOption6, 
		bool useLogScaleX, double xMin, double xMax, unsigned numBinsX, const std::string& xAxisTitle, double xAxisOffset,
		bool useLogScaleY, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		double legendX0, double legendY0, 
		const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.16);
  canvas->SetBottomMargin(0.12);
  canvas->SetLogx(useLogScaleX);
  canvas->SetLogy(useLogScaleY);

  //int colors[6] = { kMagenta - 7, kBlue - 7, kGreen - 6, kCyan - 6, kRed - 6, kBlack };
  int colors[6] = { kMagenta - 7, kBlue - 7, kCyan - 6, kGreen - 6, kRed - 6, kBlack };
  int markerStyles[6] = { 20, 21, 22, 23, 24, 25 };
  int markerSizes[6] = { 2, 2, 2, 2, 2, 2 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.21, legendY0 + 0.18, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numBinsX, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw("axis");

  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetMarkerSize(markerSizes[0]);
  graph1->Draw(drawOption1.data());
  legend->AddEntry(graph1, legendEntry1.data(), legendOption1.data());

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetMarkerSize(markerSizes[1]);
    graph2->Draw(drawOption2.data());
    legend->AddEntry(graph2, legendEntry2.data(), legendOption2.data());
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetMarkerSize(markerSizes[2]);
    graph3->Draw(drawOption3.data());
    legend->AddEntry(graph3, legendEntry3.data(), legendOption3.data());
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetMarkerSize(markerSizes[3]);
    graph4->Draw(drawOption4.data());
    legend->AddEntry(graph4, legendEntry4.data(), legendOption4.data());
  }

  if ( graph5 ) {
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetMarkerSize(markerSizes[4]);
    graph5->Draw(drawOption5.data());
    legend->AddEntry(graph5, legendEntry5.data(), legendOption5.data());
  }

  if ( graph6 ) {
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetMarkerSize(markerSizes[5]);
    graph6->Draw(drawOption6.data());
    legend->AddEntry(graph6, legendEntry6.data(), legendOption6.data());
  }
  
  legend->Draw();

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScaleY ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete legend;
  delete dummyHistogram;
  delete canvas;  
}
//-------------------------------------------------------------------------------

double square(double x)
{
  return x*x;
}

void makePlots_testHttXsectionWithTauDecays()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  typedef std::vector<std::string> vstring;
  vstring channels;
  channels.push_back("hadhad");
  channels.push_back("lephad");
  channels.push_back("hadlep");
  channels.push_back("leplep");

  std::map<std::string, std::string> inputFileNames; // key = channel
  inputFileNames["hadhad"] = "../testHttXsectionWithTauDecays_hadhad.root";
  inputFileNames["lephad"] = "../testHttXsectionWithTauDecays_lephad.root";
  inputFileNames["hadlep"] = "../testHttXsectionWithTauDecays_hadlep.root";
  inputFileNames["leplep"] = "../testHttXsectionWithTauDecays_leplep.root";

  typedef std::vector<double> vdouble;
  vdouble mH;
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
  
  typedef std::vector<unsigned> vunsigned;
  vunsigned numCalls;
  numCalls.push_back(10000);
  numCalls.push_back(20000);
  numCalls.push_back(50000);
  numCalls.push_back(100000);
  numCalls.push_back(200000);
  numCalls.push_back(500000);
  numCalls.push_back(1000000);
  numCalls.push_back(2000000);
  numCalls.push_back(5000000);
  
  enum { kVEGAS = 1, kVAMP = 2 };

  typedef std::vector<int> vint;
  vint intModes;
  intModes.push_back(kVEGAS); // VEGAS integration algorithm
  intModes.push_back(kVAMP);  // VAMP  

  typedef std::map<int, double> map1;
  typedef std::map<int, map1> map2;
  typedef std::map<int, map2> map3;
  typedef std::map<std::string, map3> map4;
  map4 xSections;               // key = channel, intMode, numCalls, mH
  map4 xSectionErrs;            // key = channel, intMode, numCalls, mH
  map4 xSections_literature;    // key = channel, intMode, numCalls, mH
  map4 xSectionErrs_literature; // key = channel, intMode, numCalls, mH

  std::string graphName = "graph_Xsection_times_BR_woAcc_numCalls%u_intMode%i";
  std::string graphName_literature = "graph_Xsection_times_BR_woAcc_literature";

  for ( vstring::const_iterator channel = channels.begin();
	channel != channels.end(); ++channel ) {
    for ( std::vector<unsigned>::const_iterator numCalls_i = numCalls.begin();
	  numCalls_i != numCalls.end(); ++numCalls_i ) {
      for ( std::vector<int>::const_iterator intMode = intModes.begin();
	    intMode != intModes.end(); ++intMode ) {
	for ( vdouble::const_iterator mH_i = mH.begin();
	      mH_i != mH.end(); ++mH_i ) {
	  std::string graphName_i = Form(graphName.data(), *numCalls_i, *intMode);
	  fillMap(inputFileNames[*channel], graphName_i, xSections[*channel][*intMode][*numCalls_i], xSectionErrs[*channel][*intMode][*numCalls_i]);
	  fillMap(inputFileNames[*channel], graphName_literature, xSections_literature[*channel][*intMode][*numCalls_i], xSectionErrs_literature[*channel][*intMode][*numCalls_i]);
	}
      }
    }
  }

  typedef std::map<int, TGraphErrors*> map1gr;
  typedef std::map<int, map1gr> map2gr;
  typedef std::map<std::string, map2gr> map3gr;
  map3gr graphs_xSection_div_literature;                   // key = channel, intMode, mH
  map3gr graphs_xSectionErr_div_literature;                // key = channel, intMode, mH
  map3gr graphs_xSection_minus_literature_div_xSectionErr; // key = channel, intMode, mH
  map3gr graphs_xSection_minus_ref_div_xSectionErr;        // key = channel, intMode, mH

  for ( vstring::const_iterator channel = channels.begin();
	channel != channels.end(); ++channel ) {
    for ( vdouble::const_iterator mH_i = mH.begin();
	      mH_i != mH.end(); ++mH_i ) {
      int mH_int = TMath::Nint(*mH_i);
      for ( std::vector<int>::const_iterator intMode = intModes.begin();
	    intMode != intModes.end(); ++intMode ) {
	TGraphErrors* graph_xSection_div_literature = new TGraphErrors(numCalls.size());
	TGraphErrors* graph_xSectionErr_div_literature = new TGraphErrors(numCalls.size());
	TGraphErrors* graph_xSection_minus_literature_div_xSectionErr = new TGraphErrors(numCalls.size());
	TGraphErrors* graph_xSection_minus_ref_div_xSectionErr = new TGraphErrors(numCalls.size());
	int idxPoint = 0;
	for ( std::vector<unsigned>::const_iterator numCalls_i = numCalls.begin();
	      numCalls_i != numCalls.end(); ++numCalls_i ) {
	  double xSection = xSections[*channel][*intMode][*numCalls_i][mH_int];
	  double xSectionErr = xSectionErrs[*channel][*intMode][*numCalls_i][mH_int];
	  double xSectionErr2 = square(xSectionErr);
	  double xSection_literature = xSections_literature[*channel][*intMode][*numCalls_i][mH_int];
	  double xSectionErr_literature = xSectionErrs_literature[*channel][*intMode][*numCalls_i][mH_int];
	  double xSectionErr2_literature = square(xSectionErr_literature);
	  //double xSection_ref = xSection_literature;
	  //double xSectionErr_ref = xSectionErr_literature;
	  unsigned numCalls_ref = numCalls[numCalls.size() - 1];
	  double xSection_ref = xSections[*channel][*intMode][numCalls_ref][mH_int];
	  double xSectionErr_ref = xSectionErrs[*channel][*intMode][numCalls_ref][mH_int];
	  double xSectionErr2_ref = square(xSectionErr_ref);
	  if ( xSection_literature > 0. ) {
	    graph_xSection_div_literature->SetPoint(idxPoint, *numCalls_i, xSection/xSection_literature);
	    graph_xSectionErr_div_literature->SetPoint(idxPoint, *numCalls_i, TMath::Sqrt(xSectionErr2)/xSection_literature);
	  }
	  if ( xSectionErr2 > 0. ) {
	    graph_xSection_minus_literature_div_xSectionErr->SetPoint(idxPoint, *numCalls_i, TMath::Abs(xSection - xSection_literature)/TMath::Sqrt(xSectionErr2 + xSectionErr2_literature));
	    graph_xSection_minus_ref_div_xSectionErr->SetPoint(idxPoint, *numCalls_i, TMath::Abs(xSection - xSection_ref)/TMath::Sqrt(xSectionErr2 + xSectionErr2_ref));
	  }
	  ++idxPoint;
	}
	graphs_xSection_div_literature[*channel][*intMode][mH_int] = graph_xSection_div_literature;
	graphs_xSectionErr_div_literature[*channel][*intMode][mH_int] = graph_xSectionErr_div_literature;
	graphs_xSection_minus_literature_div_xSectionErr[*channel][*intMode][mH_int] = graph_xSection_minus_literature_div_xSectionErr;
	graphs_xSection_minus_ref_div_xSectionErr[*channel][*intMode][mH_int] = graph_xSection_minus_ref_div_xSectionErr;
      }
    }
  }

  for ( vstring::const_iterator channel = channels.begin();
	channel != channels.end(); ++channel ) {
    for ( vdouble::const_iterator mH_i = mH.begin();
	  mH_i != mH.end(); ++mH_i ) {
      int mH_int = TMath::Nint(*mH_i);

      TGraphErrors* graph_xSection_div_literature_VEGAS = graphs_xSection_div_literature[*channel][kVEGAS][mH_int];
      TGraphErrors* graph_xSection_div_literature_VAMP = graphs_xSection_div_literature[*channel][kVAMP][mH_int];
      std::string outputFileName_xSection_div_literature = Form("plots/xSection_div_literature_%s_%i.png", channel->data(), mH_int);
      showGraphs(800, 600,
		 graph_xSection_div_literature_VEGAS, "pL", "VEGAS", "p", 
		 graph_xSection_div_literature_VAMP, "pL", "VAMP", "p", 
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 true, 5.e+3, 1.e+7, 10, "N_{calls}", 1.2, 
		 false, 0., 2., "#frac{Integral}{Literature}", 1.5,
		 0.68, 0.71, 
		 outputFileName_xSection_div_literature.data());

      TGraphErrors* graph_xSectionErr_div_literature_VEGAS = graphs_xSectionErr_div_literature[*channel][kVEGAS][mH_int];
      TGraphErrors* graph_xSectionErr_div_literature_VAMP = graphs_xSectionErr_div_literature[*channel][kVAMP][mH_int];
      std::string outputFileName_xSectionErr_div_literature = Form("plots/xSectionErr_div_literature_%s_%i.png", channel->data(), mH_int);
      showGraphs(800, 600,
		 graph_xSectionErr_div_literature_VEGAS, "pL", "VEGAS", "p", 
		 graph_xSectionErr_div_literature_VAMP, "pL", "VAMP", "p", 
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 true, 5.e+3, 1.e+7, 10, "N_{calls}", 1.2, 
		 true, 1.e-4, 1.e0, "#frac{#sigma(Integral)}{Literature}", 1.5,
		 0.68, 0.71, 
		 outputFileName_xSectionErr_div_literature.data());

      TGraphErrors* graph_xSection_minus_literature_div_xSectionErr_VEGAS = graphs_xSection_minus_literature_div_xSectionErr[*channel][kVEGAS][mH_int];
      TGraphErrors* graph_xSection_minus_literature_div_xSectionErr_VAMP = graphs_xSection_minus_literature_div_xSectionErr[*channel][kVAMP][mH_int];
      std::string outputFileName_xSection_minus_literature_div_xSectionErr = Form("plots/xSection_minus_literature_div_xSectionErr_%s_%i.png", channel->data(), mH_int);
      showGraphs(800, 600,
		 graph_xSection_minus_literature_div_xSectionErr_VEGAS, "pL", "VEGAS", "p", 
		 graph_xSection_minus_literature_div_xSectionErr_VAMP, "pL", "VAMP", "p", 
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 true, 5.e+3, 1.e+7, 10, "N_{calls}", 1.2, 
		 false, 0., 10., "#frac{Integral - Literature}{#sigma(Integral)}", 1.5,
		 0.68, 0.71, 
		 outputFileName_xSection_minus_literature_div_xSectionErr.data());

      TGraphErrors* graph_xSection_minus_ref_div_xSectionErr_VEGAS = graphs_xSection_minus_ref_div_xSectionErr[*channel][kVEGAS][mH_int];
      TGraphErrors* graph_xSection_minus_ref_div_xSectionErr_VAMP = graphs_xSection_minus_ref_div_xSectionErr[*channel][kVAMP][mH_int];
      std::string outputFileName_xSection_minus_ref_div_xSectionErr = Form("plots/xSection_minus_ref_div_xSectionErr_%s_%i.png", channel->data(), mH_int);
      showGraphs(800, 600,
		 graph_xSection_minus_ref_div_xSectionErr_VEGAS, "pL", "VEGAS", "p", 
		 graph_xSection_minus_ref_div_xSectionErr_VAMP, "pL", "VAMP", "p", 
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 0, "", "", "",
		 true, 5.e+3, 1.e+7, 10, "N_{calls}", 1.2, 
		 true, 1.e-2, 1.e+3, "#frac{Integral - Integral(N_{call} = 5mio)}{#sigma(Integral)}", 1.5,
		 0.68, 0.71, 
		 outputFileName_xSection_minus_ref_div_xSectionErr.data());
    }
  }
}


