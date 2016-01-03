#ifndef APFELINIT_H
#define APFELINIT_H

#include <stdlib.h>
#include <string>
using std::size_t;
using std::string;

namespace APFEL
{
  /**
   * @brief A container for the APFEL setup parameters
   *
   * Data structure for one time initialization
   */
  class apfel_param
  {
  public:
    size_t evol_pto;  //!< Perturbative order of the evolution    
    string FNS;       //!< Flavour number scheme to use in the evolution
    bool   damp;      //!< FONLL damping switch
    
    string MODEV;     //!< Mode of DGLAP solution
    double xiR;       //!< Xi_R renormalisation scale factor
    double xiF;       //!< Xi_F factorisation scale factor
    
    size_t nff;       //!< Number of flavours in FFN
    size_t nf_as;     //!< Number of flavours in alphas running
    size_t nf_pdf;    //!< Number of flavours in PDF running
    
    double Q0;        //!< Initial Q
    double QM;        //!< Maximum Q
    double alphas;    //!< Value of alpha_s
    double QREF;      //!< Reference QCD scale (Typically M_Z)
    
    bool QED;         //!< Flag to enable QED evolution
    double alpha_qed; //!< Value of alpha_qed
    double QEDREF;    //!< Reference QED scale
    
    bool SxRes;       //!< Small-x resummation switch
    string SxOrd;     //!< Small-x resummation order
    
    string HQMASS;    //!< Heavy quark mass (POLE/DGLAP)
    
    double mc;        //!< Mass of charm
    double mb;        //!< Mass of bottom
    double mt;        //!< Mass of top
    double Qmc;       //!< Charm mass reference scale
    double Qmb;       //!< Bottom mass reference scale
    double Qmt;       //!< Top mass reference scale

    double CKM[3][3]; //!< CKM matrix
    
    double mz;        //!< Z mass
    double mw;        //!< W mass
    
    double gf;        //!< G_Fermi
    double sin2tw;    //!< Sin^2 theta_w
    
    bool TMC;         //!< Target mass corrections  
    double Mt;        //!< Target mass

    bool SIA;         //!< Enable SIA evolution
    double truc_eps;  //!< Set epsilon for the trucated solution
  };

  //! Setup APFEL with the previous data structure
  void SetParam(apfel_param const& par);
};

#endif
