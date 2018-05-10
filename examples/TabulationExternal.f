************************************************************************
*
*     Tabulation.f:
*
*     Example program used for the LH benchmark.
*
************************************************************************
      program Tabulation
*
      implicit none
*
      integer ilha
      double precision Q0,Q
      double precision Q02,Q2
      double precision AlphaQCD,AlphaQED
      double precision xPDF,xgamma
      double precision eps
      double precision xlha(11)

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
*
*     Some examples ...
*
c      call SetFastEvolution(.true.)
c      call LockGrids(.true.)
c      call EnableEvolutionOperator(.true.)
c      call SetFFNS(4)
c      call SetTheory("QavDS")
c      call SetTheory("QUniD")
c      call SetPerturbativeOrder(0)
c      call SetPDFEvolution("exactalpha")
c      call SetPDFSet("NNPDF23_nlo_as_0119_qed")
c      call SetPDFSet("MRST2004qed")
      call SetPDFSet("external")
c      call SetNumberOfGrids(1)
c      call SetGridParameters(1,30,3,1d-5)
c      call SetGridParameters(2,30,3,2d-1)
c      call SetGridParameters(3,30,3,8d-1)
c      call SetPDFSet("NNPDF30_nnlo_as_0118")
c      call SetAlphaQCDRef(0.118d0,91.2d0)
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("expandalpha")
c      call SetPoleMasses(1.275d0,4.18d0,173.03d0)
c      call SetMaxFlavourPDFs(5)
c      call SetMaxFlavourAlpha(5)
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
      write(6,*) "alpha_QED(mu2F) =",AlphaQED(Q)
      write(6,*) "  "
*
      write(6,*) "Standard evolution:"
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
      do ilha=3,11
         write(6,'(es7.1,6es12.4)') 
     1         xlha(ilha),
     2         xPDF(2,xlha(ilha)) - xPDF(-2,xlha(ilha)),
     3         xPDF(1,xlha(ilha)) - xPDF(-1,xlha(ilha)),
     4         2d0 * ( xPDF(-1,xlha(ilha)) + xPDF(-2,xlha(ilha)) ),
     5         xPDF(4,xlha(ilha)) + xPDF(-4,xlha(ilha)),
     6         xPDF(0,xlha(ilha)),
     7         xgamma(xlha(ilha))
      enddo
      write(*,*) "  "
*
      end
*
************************************************************************
*
*     External PDF set that the user can feed to APFEL to use as initial
*     scale set.
*     The external set must have the name "ExternalPDFSetAPFEL" and the
*     following structure:
*
*     call ExternalSetAPFEL(x,xf)
*
*     where "x" is a double corresponding to the Bjorken's variable and
*     xf(-6:7) is an array of double that corresponds to the set of
*     distributions in the physical basis (tbar, bbar, ..., g, ..., b, t)
*     and where the component xf(7) is the photon distribution.
*     This subroutine is called only if the user sets:
*
*      call SetPDFSet("external")     
*
************************************************************************
      subroutine ExternalSetAPFEL(x,Q,xf)
*
      implicit none
**
*     Input Variables
*
      double precision x, Q
**
*     Internal Variables
*
      integer ipdf
      double precision xuv,xdv,xg,xdbar,xubar,xs,xsbar
      double precision N_uv,auv,buv,N_dv,adv,bdv,N_g,ag
      double precision bg,N_db,adb,bdb,fs
**
*     Output Variables
*
      double precision xf(-6:7)
*
*     Parameters of the User defined PDFs
*
      N_uv = 5.107200d0
      auv  = 0.8d0
      buv  = 3d0
      N_dv = 3.064320d0
      adv  = 0.8d0
      bdv  = 4d0
      N_g  = 1.7d0
      ag   = -0.1d0
      bg   = 5d0
      N_db = 0.1939875d0
      adb  = -0.1d0
      bdb  = 6d0
      fs   = 0.2d0
*
*     User defined PDFs
*
      xuv   = N_uv * x**auv * ( 1d0 - x )**buv
      xdv   = N_dv * x**adv * ( 1d0 - x )**bdv
      xg    = N_g  * x**ag  * ( 1d0 - x )**bg 
      xdbar = N_db * x**adb * ( 1d0 - x )**bdb
      xubar = xdbar * ( 1d0 - x )
      xs    = fs * ( xdbar + xubar )
      xsbar = xs
*
*     Initialize PDFs to zero
*
      do ipdf=-6,7
         xf(ipdf) = 0d0
      enddo
*
      if(x.gt.1d0) return
*
      xf(3)  = xs
      xf(2)  = xuv + xubar
      xf(1)  = xdv + xdbar
      xf(0)  = xg
      xf(-1) = xdbar
      xf(-2) = xubar
      xf(-3) = xsbar
*
      return
      end
