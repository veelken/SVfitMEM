//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.1.2, 2014-07-03
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_heft_H
#define Parameters_heft_H

#include <complex> 

#include "read_slha.h"
using namespace std; 

class Parameters_heft
{
  public:

    static Parameters_heft * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_WH1, mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt,
        mdl_ymb, aS, mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB,
        mdl_MP, mdl_conjg__CKM3x3, mdl_CKM3x3, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_sqrt__2, mdl_MH__exp__4, mdl_MT__exp__4, mdl_MH__exp__2,
        mdl_MT__exp__2, mdl_MH__exp__12, mdl_MH__exp__10, mdl_MH__exp__8,
        mdl_MH__exp__6, mdl_MT__exp__6, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee,
        mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw,
        mdl_v, mdl_ee__exp__2, mdl_MW__exp__12, mdl_MW__exp__10,
        mdl_MW__exp__8, mdl_MW__exp__6, mdl_MW__exp__4, mdl_AH, mdl_v__exp__2,
        mdl_lam, mdl_yb, mdl_yt, mdl_ytau, mdl_muH, mdl_gw__exp__2,
        mdl_cw__exp__2, mdl_sw__exp__2;
    std::complex<double> mdl_complexi; 
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2, mdl_GH, mdl_Gphi; 
    // Model couplings independent of aS
    std::complex<double> GC_86; 
    // Model couplings dependent on aS
    std::complex<double> GC_13; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_heft * instance; 
}; 

#endif  // Parameters_heft_H

