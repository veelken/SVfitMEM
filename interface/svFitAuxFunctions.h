#ifndef TauAnalysis_SVfitMEM_svFitAuxFunctions_h
#define TauAnalysis_SVfitMEM_svFitAuxFunctions_h

#include <TGraphErrors.h>
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include <vector>
#include <string>

namespace svFitMEM
{
  //-----------------------------------------------------------------------------
  // define masses, widths and lifetimes of particles
  // relevant for computing values of likelihood functions in SVfit algorithm
  //
  // NOTE: the values are taken from
  //        K. Nakamura et al. (Particle Data Group),
  //        J. Phys. G 37, 075021 (2010)
  //
  const double electronMass = 0.51100e-3; // GeV
  const double electronMass2 = electronMass*electronMass;
  const double muonMass = 0.10566; // GeV
  const double muonMass2 = muonMass*muonMass; 
  
  const double chargedPionMass = 0.13957; // GeV
  const double chargedPionMass2 = chargedPionMass*chargedPionMass;
  const double neutralPionMass = 0.13498; // GeV
  const double neutralPionMass2 = neutralPionMass*neutralPionMass;

  const double rhoMesonMass = 0.77526; // GeV
  const double rhoMesonMass2 = rhoMesonMass*rhoMesonMass;
  const double a1MesonMass = 1.230; // GeV
  const double a1MesonMass2 = a1MesonMass*a1MesonMass;

  const double tauLeptonMass = 1.77685; // GeV
  const double tauLeptonMass2 = tauLeptonMass*tauLeptonMass;
  const double tauLeptonMass3 = tauLeptonMass2*tauLeptonMass;
  const double tauLeptonMass4 = tauLeptonMass3*tauLeptonMass;
  const double cTauLifetime = 8.711e-3; // centimeters
  //-----------------------------------------------------------------------------

  inline double square(double x)
  {
    return x*x;
  }

  inline double cube(double x)
  {
    return x*x*x;
  }

  inline double fourth(double x)
  {
    return x*x*x*x;
  }
  
  inline double fifth(double x)
  {
    return x*x*x*x*x;
  }

  inline double sixth(double x)
  {
    return x*x*x*x*x*x;
  }

  /**
     \typedef SVfitStandalone::Vector
     \brief   spacial momentum vector (equivalent to reco::Candidate::Vector)
  */
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;
  /**
     \typedef SVfitStandalone::LorentzVector
     \brief   lorentz vector (equivalent to reco::Candidate::LorentzVector)
  */
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  double roundToNdigits(double, int = 3);

  TGraphErrors* makeGraph(const std::string&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

  void extractResult(TGraphErrors*, double&, double&, double&, int = 0);

  Vector normalize(const Vector&);
  double compScalarProduct(const Vector&, const Vector&);
  Vector compCrossProduct(const Vector&, const Vector&);
}

#endif
