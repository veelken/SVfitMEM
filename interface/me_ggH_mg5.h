//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.1.2, 2014-07-03
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_heft_gg_tamtap_H
#define MG5_Sigma_heft_gg_tamtap_H

#include <complex> 
#include <vector> 

#include "Parameters_heft.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > h WEIGHTED=2 HIG=1 HIW=1
// *   Decay: h > ta- ta+ WEIGHTED=2 HIW=1 HIG=1
//--------------------------------------------------------------------------

class me_ggH_mg5
{
  public:

    // Constructor.
    me_ggH_mg5() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "g g > ta- ta+ (heft)";}

    //void setHiggsMass(double mH) 
    //{ 
    //  pars->mdl_MH = mH; 
    //
    //  pars->mdl_MP = 1. * pars->mdl_MH; 
    //  pars->mdl_MH__exp__2 = pow(pars->mdl_MH, 2.); 
    //  pars->mdl_MH__exp__4 = pow(pars->mdl_MH, 4.); 
    //  pars->mdl_MH__exp__6 = pow(pars->mdl_MH, 6.); 
    //  pars->mdl_MH__exp__8 = pow(pars->mdl_MH, 8.); 
    //  pars->mdl_MH__exp__10 = pow(pars->mdl_MH, 10.); 
    //  pars->mdl_MH__exp__12 = pow(pars->mdl_MH, 12.); 
    //  pars->mdl_AH = (47. * pars->mdl_ee__exp__2 * (1. - (2. * pars->mdl_MH__exp__4)/(987. *
    //    pars->mdl_MT__exp__4) - (14. * pars->mdl_MH__exp__2)/(705. * pars->mdl_MT__exp__2) + (213.
    //    * pars->mdl_MH__exp__12)/(2.634632e7 * pars->mdl_MW__exp__12) + (5. *
    //    pars->mdl_MH__exp__10)/(119756. * pars->mdl_MW__exp__10) + (41. *
    //    pars->mdl_MH__exp__8)/(180950. * pars->mdl_MW__exp__8) + (87. *
    //    pars->mdl_MH__exp__6)/(65800. * pars->mdl_MW__exp__6) + (57. * pars->mdl_MH__exp__4)/(6580.
    //    * pars->mdl_MW__exp__4) + (33. * pars->mdl_MH__exp__2)/(470. * pars->mdl_MW__exp__2)))/(72.
    //    * pow(M_PI, 2.) * pars->mdl_v);
    //  pars->mdl_lam = pars->mdl_MH__exp__2/(2. * pars->mdl_v__exp__2); 
    //
    //  pars->setDependentParameters();
    //  pars->setDependentCouplings();
    //}
    //void setHiggsWidth(double width) 
    //{
    //  pars->mdl_WH1 = width;
    //  pars->mdl_WH = width;
    //
    //  pars->setDependentParameters();
    //  pars->setDependentCouplings();
    //}

    double getHiggsMass() const
    {
      return pars->mdl_MH;
    } 
    double getHiggsWidth() const
    {
      return pars->mdl_WH;
    } 

    virtual int code() const {return 0;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 4; 
    static const int nprocesses = 1; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 5; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 1; 
    std::complex<double> amp[namplitudes]; 
    double matrix_gg_h_h_tamtap(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_heft * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_heft_gg_tamtap_H
