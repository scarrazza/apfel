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
      double precision xPDFj,xgamma,xLepton,xf(-6:6)
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
c      call SetFFNS(3)
c      call SetTheory("QavDS")
c      call SetTheory("QED")
c      call SetTheory("QUniD")
c      call EnableLeptonEvolution(.true.)
c      call SetTauMass(1d10)
c      call SetPerturbativeOrder(0)
c      call SetPDFEvolution("exactalpha")
c      call SetPDFSet("NNPDF23_nlo_as_0119_qed.LHgrid")
c      call SetPDFSet("MRST2004qed.LHgrid")
c      call SetNumberOfGrids(1)
c      call SetGridParameters(1,30,3,1d-5)
c      call SetGridParameters(2,30,3,2d-1)
c      call SetGridParameters(3,30,3,8d-1)
c      call SetPDFSet("NNPDF30_nnlo_as_0118.LHgrid")
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
      write(6,'(a5,2a12,a14,a10,3a12,a13,a14)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon",
     2         "e^-+e^+","mu^-+mu^+","tau^-+tau^+"
      do ilha=3,11
         write(6,'(es7.1,9es12.4)') 
     1         xlha(ilha),
     2         xPDFj(2,xlha(ilha)) - xPDFj(-2,xlha(ilha)),
     3         xPDFj(1,xlha(ilha)) - xPDFj(-1,xlha(ilha)),
     4         2d0 * ( xPDFj(-1,xlha(ilha)) + xPDFj(-2,xlha(ilha)) ),
     5         xPDFj(4,xlha(ilha)) + xPDFj(-4,xlha(ilha)),
     6         xPDFj(0,xlha(ilha)),
     7         xgamma(xlha(ilha)),
     8         xLepton(1,xlha(ilha))+xLepton(-1,xlha(ilha)),
     9         xLepton(2,xlha(ilha))+xLepton(-2,xlha(ilha)),
     1         xLepton(3,xlha(ilha))+xLepton(-3,xlha(ilha))
      enddo
      write(*,*) "  "
c$$$*
c$$$      write(6,*) "Standard evolution using the xPDFall function:"
c$$$      write(6,'(a5,2a12,a14,a10,2a12)') "x",
c$$$     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
c$$$      do ilha=3,11
c$$$         call xPDFall(xlha(ilha),xf)
c$$$         write(6,'(es7.1,6es12.4)') 
c$$$     1         xlha(ilha),
c$$$     2         xf(2) - xf(-2),
c$$$     3         xf(1) - xf(-1),
c$$$     4         2d0 * ( xf(-1) + xf(-2) ),
c$$$     5         xf(4) + xf(-4),
c$$$     6         xf(0),
c$$$     7         xgamma(xlha(ilha))
c$$$      enddo
c$$$      write(*,*) "  "
*
      end
