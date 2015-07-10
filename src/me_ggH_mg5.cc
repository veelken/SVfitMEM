//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.1.2, 2014-07-03
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "TauAnalysis/SVfitMEM/interface/me_ggH_mg5.h"
#include "TauAnalysis/SVfitMEM/interface/HelAmps_heft.h"
#include "TauAnalysis/SVfitMEM/interface/svFitAuxFunctions.h"

#include <TMath.h>

using namespace MG5_heft; 
using namespace svFitMEM;

const double mtop = 172.5;

namespace LHAPDF {
  double alphasPDF(int nset, double Q);
}

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > h WEIGHTED=2 HIG=1 HIW=1
// *   Decay: h > ta- ta+ WEIGHTED=2 HIW=1 HIG=1

//--------------------------------------------------------------------------
// Constructor

me_ggH_mg5::me_ggH_mg5(bool applyNWA, bool includeHtoTauTauDecay)
  : applyNWA_(applyNWA),
    includeHtoTauTauDecay_(includeHtoTauTauDecay)
{}

//--------------------------------------------------------------------------
// Initialize process.

void me_ggH_mg5::initProc(const string& param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_heft::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  //pars->printIndependentParameters(); 
  //pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_MTA); 
  mME.push_back(pars->mdl_MTA); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void me_ggH_mg5::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    //pars->printDependentParameters(); 
    //pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 16; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1,
      -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1},
      {-1, 1, 1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1,
      1, -1}, {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1,
      1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {256}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_gg_h_h_tamtap(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_gg_h_h_tamtap(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ ) {
    matrix_element[i] /= denominators[i];
  }

  double mH = this->getHiggsMass();
  double GammaH = this->getHiggsWidth();
  const double one_over_Pi = 1./TMath::Pi();
  double GammaH_times_mH = GammaH*mH;
  if ( applyNWA_ ) {
    for (int i = 0; i < nprocesses; i++ ) {
      matrix_element[i] *= (one_over_Pi*GammaH_times_mH);
    }
  }
  if ( !includeHtoTauTauDecay_ ) {
    for (int i = 0; i < nprocesses; i++ ) {
      matrix_element[i] /= 2.*(tauLeptonMass2/v2)*square(mH)*(1. - 4.*tauLeptonMass2/square(mH));    
      matrix_element[i] *= (one_over_Pi*GammaH_times_mH);
    }
  }
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double me_ggH_mg5::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
    return matrix_element[0]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void me_ggH_mg5::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  FFS2_3(w[3], w[2], pars->GC_86, pars->mdl_MH, pars->mdl_WH, w[4]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVS3_0(w[0], w[1], w[4], pars->GC_13, amp[0]); 

}
double me_ggH_mg5::matrix_gg_h_h_tamtap() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 1; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{2}}; 

  // Calculate color flows
  jamp[0] = +2. * (-amp[0]); 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



