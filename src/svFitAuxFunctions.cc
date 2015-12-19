#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"

#include <TMath.h>
#include <TF1.h>
#include <TFitResult.h>

namespace svFitMEM
{

double roundToNdigits(double x, int n)
{
  double tmp = TMath::Power(10., n);
  if ( x != 0. ) {
    tmp /= TMath::Power(10., TMath::Floor(TMath::Log10(TMath::Abs(x))));
  }
  double x_rounded = TMath::Nint(x*tmp)/tmp;
  //std::cout << "<roundToNdigits>: x = " << x << ", x_rounded = " << x_rounded << std::endl;
  return x_rounded;
}
  
TGraphErrors* makeGraph(const std::string& graphName, const std::vector<GraphPoint>& graphPoints)
{
  //std::cout << "<makeGraph>:" << std::endl;
  size_t numPoints = graphPoints.size();
  //std::cout << " numPoints = " << numPoints << std::endl;
  TGraphErrors* graph = new TGraphErrors(numPoints);
  graph->SetName(graphName.data());
  for ( size_t iPoint = 0; iPoint < numPoints; ++iPoint ) {
    const GraphPoint& graphPoint = graphPoints[iPoint];
    graph->SetPoint(iPoint, graphPoint.x_, graphPoint.y_);
    graph->SetPointError(iPoint, graphPoint.xErr_, graphPoint.yErr_);
  }
  return graph;
}
  
void extractResult(TGraphErrors* graph, double& mass, double& massErr, double& Lmax, int verbosity)
{
  // determine range of mTest values that are within ~2 sigma interval within maximum of likelihood function
  double x_Lmax = 0.;
  double y_Lmax = 0.;
  double idxPoint_Lmax = -1;
  for ( int iPoint = 0; iPoint < graph->GetN(); ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    if ( y > y_Lmax ) {
      x_Lmax = x;
      y_Lmax = y;
      idxPoint_Lmax = iPoint;
    }    
  }

  double xMin = 1.e+6;
  double xMax = 0.;
  for ( int iPoint = 0; iPoint < graph->GetN(); ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    if ( x < xMin ) xMin = x;
    if ( x > xMax ) xMax = x;
  }
  
  // fit log-likelihood function within ~2 sigma interval within maximum 
  // with parabola
  std::vector<GraphPoint> graphPoints_forFit;
  double xMin_fit = 1.e+6;
  double xMax_fit = 0.;
  for ( int iPoint = 0; iPoint < graph->GetN(); ++iPoint ) {
    double x, y;
    graph->GetPoint(iPoint, x, y);
    double xErr = graph->GetErrorX(iPoint);
    double yErr = graph->GetErrorY(iPoint);
    //std::cout << "point #" << iPoint << ": x = " << x << " +/- " << xErr << ", y = " << y << " +/- " << yErr << std::endl;
    if ( y > (1.e-1*y_Lmax) && TMath::Abs(iPoint - idxPoint_Lmax) <= 5 ) {
      GraphPoint graphPoint;
      graphPoint.x_ = x;
      graphPoint.xErr_ = xErr;
      if ( (x - xErr) < xMin_fit ) xMin_fit = x - xErr; 
      if ( (x + xErr) > xMax_fit ) xMax_fit = x + xErr; 
      graphPoint.y_ = -TMath::Log(y);
      graphPoint.yErr_ = yErr/y;
      graphPoints_forFit.push_back(graphPoint);
    }
  }

  TGraphErrors* likelihoodGraph_forFit = makeGraph("svFitLikelihoodGraph_forFit", graphPoints_forFit);
  int numPoints = likelihoodGraph_forFit->GetN();
  bool useFit = false;
  if ( numPoints >= 3 ) {  
    TF1* fitFunction = new TF1("fitFunction", "TMath::Power((x - [0])/[1], 2.) + [2]", xMin_fit, xMax_fit);
    fitFunction->SetParameter(0, x_Lmax);
    fitFunction->SetParameter(1, 0.20*x_Lmax);
    fitFunction->SetParameter(2, -TMath::Log(y_Lmax));
    
    std::string fitOptions = "NSQ";
    //if ( !verbosity ) fitOptions.append("Q");
    TFitResultPtr fitResult = likelihoodGraph_forFit->Fit(fitFunction, fitOptions.data());
    if ( fitResult.Get() ) {
      if ( verbosity >= 1 ) {
	std::cout << "fitting graph of p versus M(test) in range " << xMin_fit << ".." << xMax_fit << ", result:" << std::endl;
	std::cout << " parameter #0 = " << fitFunction->GetParameter(0) << " +/- " << fitFunction->GetParError(0) << std::endl;
	std::cout << " parameter #1 = " << fitFunction->GetParameter(1) << " +/- " << fitFunction->GetParError(1) << std::endl;
	std::cout << " parameter #2 = " << fitFunction->GetParameter(2) << " +/- " << fitFunction->GetParError(2) << std::endl;
	std::cout << "chi^2 = " << fitResult->Chi2() << std::endl;
      }
      if ( fitResult->Chi2() < (10.*numPoints) && 
	   fitFunction->GetParameter(0) > xMin && fitFunction->GetParameter(0) < xMax &&
	   TMath::Abs(fitFunction->GetParameter(0) - x_Lmax) < (0.10*x_Lmax) ) {
	mass = fitFunction->GetParameter(0);
	massErr = TMath::Sqrt(square(fitFunction->GetParameter(1)) + square(fitFunction->GetParError(0)));
	Lmax = TMath::Exp(-fitFunction->GetParameter(2));
	//std::cout << "fit: mass = " << mass << " +/- " << massErr << " (Lmax = " << Lmax << ")" << std::endl;
	useFit = true;
      }
    } else {
      std::cerr << "Warning in <extractResult>: Fit did not converge !!" << std::endl;
    }
    delete fitFunction;
  }
  if ( !useFit ) {
    mass = x_Lmax;
    massErr = TMath::Sqrt(0.5*(square(x_Lmax - xMin_fit) + square(xMax_fit - x_Lmax)))/TMath::Sqrt(2.*TMath::Log(10.));
    Lmax = y_Lmax;
    //std::cout << "graph: mass = " << mass << " +/- " << massErr << " (Lmax = " << Lmax << ")" << std::endl;
  }
  
  delete likelihoodGraph_forFit;
}

Vector normalize(const Vector& p)
{
  double p_x = p.x();
  double p_y = p.y();
  double p_z = p.z();
  double mag2 = square(p_x) + square(p_y) + square(p_z);
  if ( mag2 <= 0. ) return p;
  double mag = TMath::Sqrt(mag2);
  return Vector(p_x/mag, p_y/mag, p_z/mag);
}

double compScalarProduct(const Vector& p1, const Vector& p2)
{
  return (p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z());
}

Vector compCrossProduct(const Vector& p1, const Vector& p2)
{
  double p3_x = p1.y()*p2.z() - p1.z()*p2.y();
  double p3_y = p1.z()*p2.x() - p1.x()*p2.z();
  double p3_z = p1.x()*p2.y() - p1.y()*p2.x();
  return Vector(p3_x, p3_y, p3_z);
}

}
