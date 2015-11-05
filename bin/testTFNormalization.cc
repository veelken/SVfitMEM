
/**
   \class testTFNormalization testTFNormalization.cc "TauAnalysis/SVfitStandalone/bin/testTFNormalization.cc"
   \brief Verify that normalization of transfer functions is correct
*/

#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"
#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <TMath.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TH1.h>
#include <TAxis.h>

using namespace svFitMEM;

enum { kTF_met, kTF_lepTauDecay, kTF_hadTauDecay, kTF_hadTauEn };

enum { kElectron, kMuon };

double compIntegral(double (*g)(double* x, size_t dim, void* params), int numDimensions, double* xl, double* xu, void* params, double& integralErr)
{
  gsl_monte_function* vegasIntegrand = new gsl_monte_function;
  vegasIntegrand->f = g;
  vegasIntegrand->dim = numDimensions;
  vegasIntegrand->params = params;
  gsl_monte_vegas_state* vegasWorkspace = gsl_monte_vegas_alloc(numDimensions);
  gsl_rng_env_setup();
  gsl_rng* vegasRnd = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(vegasRnd, 12345); 
  gsl_monte_vegas_init(vegasWorkspace);
  double integral = 0.;
  integralErr = 0.;
  vegasWorkspace->stage = 0;
  const int numCallsGridOpt = 100000; 
  gsl_monte_vegas_integrate(vegasIntegrand, xl, xu, numDimensions, numCallsGridOpt, vegasRnd, vegasWorkspace, &integral, &integralErr);
  integral = 0.;
  integralErr = 0.;
  vegasWorkspace->stage = 1;
  const int numCallsIntEval = 400000;
  gsl_monte_vegas_integrate(vegasIntegrand, xl, xu, numDimensions, numCallsIntEval, vegasRnd, vegasWorkspace, &integral, &integralErr);
  delete vegasIntegrand;
  gsl_monte_vegas_free(vegasWorkspace);
  gsl_rng_free(vegasRnd);
  return integral;
}

double gTF_met(double* x, size_t dim, void* params_void)
{
  double* params = (double*)params_void;
  double residualX = x[0] - params[0];
  double residualY = x[1] - params[1];
  double invCovMETxx = params[2];
  double invCovMETxy = params[3];
  double invCovMETyx = params[4];
  double invCovMETyy = params[5];
  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) + residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  double const_MET = params[6];
  double prob = const_MET*TMath::Exp(-0.5*pull2);
  return prob;
}

double gTF_hadTauDecay(double* x, size_t dim, void* params_void)
{
  double* params = (double*)params_void;
  double tauP = params[0];
  double tauP2 = square(tauP);
  double tauEn2 = tauP2 + tauLeptonMass2;
  double tauEn = TMath::Sqrt(tauEn2);
  //double tauEta = params[1];
  //double tauPhi = params[2];
  double visMass = params[3];
  double visMass2 = square(visMass);
  double X = x[0];
  double nuEn = (1. - X)*tauEn;
  double nuP = nuEn; // nuMass = 0
  double nuP2 = square(nuP);
  //double phiNu = x[1];
  double cosThetaNu = ((1. - X)*tauEn2 - 0.5*(tauLeptonMass2 - visMass2))/(tauP*nuP);
  if ( !(cosThetaNu >= -1. && cosThetaNu <= +1.) ) return 0.;
  double visEn = TMath::Sqrt(tauP2 + nuP2 - 2.*tauP*nuP*cosThetaNu + visMass2);
  double prob = (1./(32.*square(TMath::Pi())*tauEn))*(1./visEn)*(X*tauEn2/tauP);
  //-------------------------------------------------------------------------
  // CV: multiply by constant matrix element, 
  //     chosen such that the branching fraction of the tau to decay into hadrons is reproduced      
  const double M2 = 16.*TMath::Pi()*cube(tauLeptonMass)*GammaTauToHad/(tauLeptonMass2 - visMass2);
  prob *= M2;
  //-------------------------------------------------------------------------
  return prob;
}

double GammaLep(int lepType) 
{
  double GammaLep = ( lepType == kMuon ) ? GammaTauToMu : GammaTauToElec;
  return GammaLep;
}

double gTF_lepTauDecay(double* x, size_t dim, void* params_void)
{
  double* params = (double*)params_void;
  double tauP = params[0];
  double tauP2 = square(tauP);
  double tauEn2 = tauP2 + tauLeptonMass2;
  double tauEn = TMath::Sqrt(tauEn2);
  //std::cout << "tauEn = " << tauEn << std::endl;
  //double tauEta = params[1];
  //double tauPhi = params[2];
  double visMass = params[3];
  double visMass2 = square(visMass);
  double X = x[0];
  double nunuEn = (1. - X)*tauEn;
  //std::cout << "X = " << X << ": nunuEn = " << nunuEn << std::endl;
  double nunuEn2 = square(nunuEn);
  double nunuMass2 = x[2];
  double nunuMass = TMath::Sqrt(nunuMass2);
  double nunuP2 = nunuEn2 - nunuMass2;
  if ( !(nunuP2 > 0.) ) return 0.;
  double nunuP = TMath::Sqrt(nunuP2);
  //double phiNu = x[1];
  double cosThetaNuNu = ((1. - X)*tauEn2 - 0.5*(tauLeptonMass2 + nunuMass2 - visMass2))/(tauP*nunuP);
  if ( !(cosThetaNuNu >= -1. && cosThetaNuNu <= +1.) ) return 0.;
  double visEn = TMath::Sqrt(tauP2 + nunuP2 - 2.*tauP*nunuP*cosThetaNuNu + visMass2);
  //std::cout << "cosThetaNuNu = " << cosThetaNuNu << ", visEn = " << visEn << std::endl;
  const double GFfactor = square(GF)/(2.*square(TMath::Pi()));
  double tauEn_rf = (tauLeptonMass2 + nunuMass2 - visMass2)/(2.*nunuMass);
  double visEn_rf = tauEn_rf - nunuMass;
  if ( !(tauEn_rf > 0. && visEn_rf > 0.) ) return 0.;
  double I = GFfactor*nunuMass2*(2.*tauEn_rf*visEn_rf - (2./3.)*TMath::Sqrt((square(tauEn_rf) - tauLeptonMass2)*(square(visEn_rf) - visMass2)));
  double prob = (1./(32.*square(TMath::Pi())*tauEn))*(1./visEn)*I*(X*tauEn2/tauP);
  return prob;
}

double gTF_hadTauEn(double* x, size_t dim, void* params_void)
{
  static HadTauTFCrystalBall* gHadTauTFCrystalBall = 0;
  if ( !gHadTauTFCrystalBall ) {
    gHadTauTFCrystalBall = new HadTauTFCrystalBall();
  }
  double* params = (double*)params_void;
  double genTauPt = params[0];
  double genTauEta = params[1];
  int decayMode = TMath::Nint(params[2]);
  gHadTauTFCrystalBall->setDecayMode(decayMode);
  double rec_div_genTauPt = x[0];
  double recTauPt = rec_div_genTauPt*genTauPt;
  double prob = (*gHadTauTFCrystalBall)(recTauPt, genTauPt, genTauEta);
  //static int idxCall = 0;
  //std::cout << "<gTF_hadTauEn (call #" << idxCall << "): prob = " << prob << std::endl;
  //++idxCall;
  return prob;
}

void fillWithOverFlow(TH1* histogram, double x, double weight = 1.)
{
  TAxis* xAxis = histogram->GetXaxis();
  int bin = xAxis->FindBin(x);
  int numBins = xAxis->GetNbins();
  if ( bin < 1       ) bin = 1;
  if ( bin > numBins ) bin = numBins;
  double xBin = xAxis->GetBinCenter(bin);
  histogram->Fill(xBin, weight);
}

int main(int argc, char* argv[]) 
{
  std::vector<int> checksToRun;
  //checksToRun.push_back(kTF_met);
  //checksToRun.push_back(kTF_lepTauDecay);
  //checksToRun.push_back(kTF_hadTauDecay);
  checksToRun.push_back(kTF_hadTauEn);

  TH1* norm_met = new TH1D("norm_met", "norm_met", 10000, 0., 10.);
  TH1* norm_lepTauDecay = new TH1D("norm_lepTauDecay", "norm_lepTauDecay", 10000, 0., 10.);
  TH1* norm_hadTauDecay = new TH1D("norm_hadTauDecay", "norm_hadTauDecay", 10000, 0., 10.);
  TH1* norm_hadTauEn = new TH1D("norm_hadTauEn", "norm_hadTauEn", 10000, 0., 10.);

  int idxCheck = 0;
  for ( std::vector<int>::const_iterator checkToRun = checksToRun.begin();
	checkToRun != checksToRun.end(); ++checkToRun ) {
    std::cout << "running check #" << idxCheck << ":" << std::endl;
    if ( (*checkToRun) == kTF_met ) {
      std::cout << "checking MET transfer function..." << std::endl;
      int numToys = 1000;
      double metResolution = 10.;
      TRandom3 rnd;
      for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
	double metPx_gen = rnd.Gaus(0., 10.);
	double metPy_gen = rnd.Gaus(0., 10.);
	//std::cout << "met (gen): Px = " << metPx_gen << ", Py = " << metPy_gen << std::endl;
	double metCov00  = square(rnd.Gaus(0., 10.));
	double metCov11  = square(rnd.Gaus(0., 10.));
	double metCov01  = TMath::Sqrt(metCov00*metCov11)*rnd.Uniform(-1., +1.);
	TMatrixD covMET(2,2);
	covMET(0,0) = metCov00;
	covMET(0,1) = metCov01;
	covMET(1,0) = metCov01;
	covMET(1,1) = metCov11;
	//std::cout << "covMET:" << std::endl;
	//covMET.Print();
	double covDet = covMET.Determinant();
	if ( covDet == 0 ) {
	  std::cerr << "Warning: MET covariance matrix cannot be inverted --> skipping !!" << std::endl;
	  continue;
	}
	TMatrixD invCovMET(covMET);
	invCovMET.Invert(); 
	double* xl = new double[2];
	double* xu = new double[2];
	xl[0] = -10.*metResolution;
	xu[0] = +10.*metResolution;
	xl[1] = -10.*metResolution;
	xu[1] = +10.*metResolution;
	double* params = new double[7];
	params[0] = metPx_gen;
	params[1] = metPy_gen;
	params[2] = invCovMET(0,0);
	params[3] = invCovMET(0,1);
	params[4] = invCovMET(1,0);
	params[5] = invCovMET(1,1);
	params[6] = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));
	double normalizationErr;
	double normalization = compIntegral(&gTF_met, 2, xl, xu, params, normalizationErr);
	delete [] xl;
	delete [] xu;
	delete [] params;
	fillWithOverFlow(norm_met, normalization);
	std::cout << " toy #" << idxToy << ": normalization = " << normalization << " +/- " << normalizationErr 
		  << " (expected = 1.0, ratio = " << normalization << " +/- " << normalizationErr << ")" << std::endl;
      }
    } else if ( (*checkToRun) == kTF_lepTauDecay ) {
      std::cout << "checking transfer function for leptonic tau decays..." << std::endl;
      int numToys = 1000;
      TRandom3 rnd;
      for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
	double tauP = rnd.Uniform(0., 100.);
	double tauEta = rnd.Uniform(-2.3, +2.3);
	double tauPhi = rnd.Uniform(-TMath::Pi(), +TMath::Pi());
	//double tauP = 100.;
	//double tauEta = 0.;
	//double tauPhi = 0.;
	std::cout << "tau: P = " << tauP << ", eta = " << tauEta << ", phi = " << tauPhi << std::endl;		
	int lepType = ( rnd.Uniform(0., 1.) > 0.5 ) ? kMuon : kElectron;
	double visMass = ( lepType == kMuon ) ? muonMass : electronMass;
	//double visMass = muonMass;
	std::cout << "visMass = " << visMass << std::endl;
	double* xl = new double[3];
	double* xu = new double[3];
	xl[0] = 0.;
	xu[0] = 1.;
	xl[1] = -TMath::Pi();
	xu[1] = +TMath::Pi();
	xl[2] = 0.;
	xu[2] = tauLeptonMass2;
	double* params = new double[4];
	params[0] = tauP;
	params[1] = tauEta;
	params[2] = tauPhi;
	params[3] = visMass;
	double normalizationErr;
	double normalization = compIntegral(&gTF_lepTauDecay, 3, xl, xu, params, normalizationErr);
	double tauEn = TMath::Sqrt(square(tauP) + tauLeptonMass2);
	double gamma = tauEn/tauLeptonMass;
	double GammaTauToLep_div_gamma = GammaLep(lepType)/gamma;
	delete [] xl;
	delete [] xu;
	delete [] params;
	fillWithOverFlow(norm_lepTauDecay, normalization/GammaTauToLep_div_gamma);
	std::cout << " toy #" << idxToy << ": normalization = " << normalization << " +/- " << normalizationErr 
		  << " (expected = " << GammaTauToLep_div_gamma << ", ratio = " << (normalization/GammaTauToLep_div_gamma) << " +/- " << (normalizationErr/GammaTauToLep_div_gamma) << ")" << std::endl;
      }
    } else if ( (*checkToRun) == kTF_hadTauDecay ) {
      std::cout << "checking transfer function for hadronic tau decays..." << std::endl;
      int numToys = 1000;
      TRandom3 rnd;
      for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
	double tauP = rnd.Uniform(0., 100.);
	double tauEta = rnd.Uniform(-2.3, +2.3);
	double tauPhi = rnd.Uniform(-TMath::Pi(), +TMath::Pi());
	//double tauP = 100.;
	//double tauEta = 0.;
	//double tauPhi = 0.;
	std::cout << "tau: P = " << tauP << ", eta = " << tauEta << ", phi = " << tauPhi << std::endl;
	double visMass = rnd.Uniform(chargedPionMass, 1.5);
	//double visMass = chargedPionMass;
	std::cout << "visMass = " << visMass << std::endl;
	double* xl = new double[2];
	double* xu = new double[2];
	xl[0] = 0.;
	xu[0] = 1.;
	xl[1] = -TMath::Pi();
	xu[1] = +TMath::Pi();
	double* params = new double[4];
	params[0] = tauP;
	params[1] = tauEta;
	params[2] = tauPhi;
	params[3] = visMass;
	double normalizationErr;
	double normalization = compIntegral(&gTF_hadTauDecay, 2, xl, xu, params, normalizationErr);
	double tauEn = TMath::Sqrt(square(tauP) + tauLeptonMass2);
	double gamma = tauEn/tauLeptonMass;
	double GammaTauToHad_div_gamma = GammaTauToHad/gamma;
	delete [] xl;
	delete [] xu;
	delete [] params;
	fillWithOverFlow(norm_hadTauDecay, normalization/GammaTauToHad_div_gamma);
	std::cout << " toy #" << idxToy << ": normalization = " << normalization << " +/- " << normalizationErr 
		  << " (expected = " << GammaTauToHad_div_gamma << ", ratio = " << (normalization/GammaTauToHad_div_gamma) << " +/- " << (normalizationErr/GammaTauToHad_div_gamma) << ")" << std::endl;
      }
    } else if ( (*checkToRun) == kTF_hadTauEn ) {
      std::cout << "checking transfer function for hadronic tau energy reconstruction..." << std::endl;
      int numToys = 1000;
      TRandom3 rnd;
      for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
	double tauPt = rnd.Uniform(20., 100.);
	double tauEta = rnd.Uniform(-2.3, +2.3);
	double tauPhi = rnd.Uniform(-TMath::Pi(), +TMath::Pi());
	std::cout << "tau: Pt = " << tauPt << ", eta = " << tauEta << ", phi = " << tauPhi << std::endl;
	double visMass = rnd.Uniform(chargedPionMass, 1.5);
	std::cout << "visMass = " << visMass << std::endl;
	int decayMode = -1;
	while ( !(decayMode == 0 || decayMode == 1 || decayMode == 2 || decayMode == 10) ) {
	  decayMode = TMath::Nint(rnd.Uniform(-0.5, +10.5));
	}
	std::cout << "decayMode = " << decayMode << std::endl;
	double* xl = new double[1];
	double* xu = new double[1];
	xl[0] = 0.;
	xu[0] = 2.;
	double* params = new double[3];
	params[0] = tauPt;
	params[1] = tauEta;
	params[2] = decayMode;
	double normalizationErr;
	double normalization = compIntegral(&gTF_hadTauEn, 1, xl, xu, params, normalizationErr);
	delete [] xl;
	delete [] xu;
	delete [] params;
	fillWithOverFlow(norm_hadTauEn, normalization);
	std::cout << " toy #" << idxToy << ": normalization = " << normalization << " +/- " << normalizationErr 
		  << " (expected = 1.0, ratio = " << normalization << " +/- " << normalizationErr << ")" << std::endl;
      }
    } else {
      std::cerr << "Check of type = " << (*checkToRun) << " is not defined !!" << std::endl;
      assert(0);
    }
    ++idxCheck;
  }

  TFile* outputFile = new TFile("testTFNormalization.root", "RECREATE");
  norm_met->Write();
  norm_lepTauDecay->Write();
  norm_hadTauDecay->Write();
  norm_hadTauEn->Write();
  delete outputFile;

  return 0;
}
