************************************************************************
*
*     TabulationFrag.f:
*
*     Example program used for fragmentation function evolution.
*     The fragmentation functions are hardcoded and for the following
*     options are available:
*                          
*     DSS                  HKNS:
*                          
*     - "dss_pip_lo"       - "hknsff07_pip_lo"
*     - "dss_pip_nlo"      - "hknsff07_pip_nlo"
*     - "dss_pim_lo"       - "hknsff07_pim_lo"
*     - "dss_pim_nlo"      - "hknsff07_pim_nlo"
*     - "dss_Kp_lo"        - "hknsff07_Kp_lo"
*     - "dss_Kp_nlo"       - "hknsff07_Kp_nlo"
*     - "dss_Km_lo"        - "hknsff07_Km_lo"
*     - "dss_Km_nlo"       - "hknsff07_Km_nlo"
*     - "dss_p_lo"         - "hknsff07_p_lo"
*     - "dss_p_nlo"        - "hknsff07_p_nlo"
*     - "dss_pb_lo"        - "hknsff07_pb_lo"
*     - "dss_pb_nlo"       - "hknsff07_pb_nlo"
*     - "dss_h_lo"    
*     - "dss_h_nlo"   
*     - "dss_hb_lo"   
*     - "dss_hb_nlo"  
*     
************************************************************************
      program TabulationFrag
*
      implicit none
*
      integer ilha
      double precision Q0,Q
      double precision Q02,Q2
      double precision AlphaQCD
      double precision xPDF
      double precision eps
      double precision xlha(11)
      parameter(eps=1d-10)
      data xlha / 1d-2, 5d-2, 1d-1, 2d-1, 3d-1, 4d-1,
     1            5d-1, 6d-1, 7d-1, 8d-1, 9d-1 /
*
*     Settings
*
      call SetTimeLikeEvolution(.true.)
*
      call SetPDFSet("kretzer")
      call SetPerturbativeOrder(0)
      call SetMaxFlavourPDFs(5)
      call SetMaxFlavourAlpha(5)
      call SetPoleMasses(1.43d0,4.3d0,175d0)
      call SetAlphaEvolution("lambda")
      call SetLambdaQCDRef(0.220d0,4)
c      call SetLambdaQCDRef(0.323d0,4)
      call SetNumberOfGrids(2)
      call SetGridParameters(1,50,3,1d-2)
      call SetGridParameters(2,40,3,7d-1)
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Evolve PDFs on the grids
*
      write(6,*) "Enter initial and final scale in GeV^2"
      read(5,*) Q02,Q2
*
      Q0 = dsqrt(Q02) - eps
      Q  = dsqrt(Q2)
      call EvolveAPFEL(Q0,Q)
*
*     Tabulate PDFs for the LHA x values
*
      write(6,*) "alpha_QCD(mu2F) =",AlphaQCD(Q)
      write(6,*) "  "
*
      write(6,*) "Standard evolution:"
      write(6,"(a5,2a12,a14,a10,a12)") "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
      do ilha=2,11
         write(6,"(es7.1,5es12.4)") 
     1         xlha(ilha),
     2         xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3         xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4         2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5         xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6         xPDF(0,xlha(ilha))
      enddo
      write(*,*) "  "
*
      end
