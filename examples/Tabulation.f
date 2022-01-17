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
      double precision xPDFj,xgammaj,xLeptonj,xf(-6:6),xPDFxQ
      double precision eps
      double precision xlha(11)

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
*
*     Some examples ...
*
c      call SetMassMatchingScales(2d0,1d0,1d0)
c      call SetMSbarMasses(dsqrt(2d0),4.5d0,175d0)
c      call SetMSbarMasses(1.275d0,4.18d0,173.03d0)
c      call SetFastEvolution(.true.)
c      call LockGrids(.true.)
c      call EnableEvolutionOperator(.true.)
c      call SetFFNS(3)
c      call SetTheory("QCD")
c      call SetSmallxResummation(0, "NLL")
c      call EnableNLOQEDCorrections(.true.)
c      call EnableLeptonEvolution(.true.)
c      call SetTauMass(1d10)
c      call SetPerturbativeOrder(0)
c      call SetPDFEvolution("exactalpha")
c      call SetPDFSet("")
c      call SetPDFSet("NNPDF31_nnlo_as_0118")
c      call SetNumberOfGrids(1)
c      call SetGridParameters(1,30,3,1d-5)
c      call SetGridParameters(2,30,3,2d-1)
c      call SetGridParameters(3,30,3,8d-1)
c      call SetPDFSet("MSHT20nnlo_as118")
c      call SetAlphaQCDRef(0.118d0,91.2d0)
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("exactalpha")
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
     7         xgammaj(xlha(ilha)),
     8         xLeptonj(1,xlha(ilha))+xLeptonj(-1,xlha(ilha)),
     9         xLeptonj(2,xlha(ilha))+xLeptonj(-2,xlha(ilha)),
     1         xLeptonj(3,xlha(ilha))+xLeptonj(-3,xlha(ilha))
      enddo
      write(*,*) "  "
*
      write(6,*) "Standard evolution using the xPDFall function:"
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
      do ilha=3,11
         call xPDFall(xlha(ilha),xf)
         write(6,'(es7.1,6es12.4)') 
     1         xlha(ilha),
     2         xf(2) - xf(-2),
     3         xf(1) - xf(-1),
     4         2d0 * ( xf(-1) + xf(-2) ),
     5         xf(4) + xf(-4),
     6         xf(0)
      enddo
      write(*,*) "  "
*
*     Cached PDFs
*
      call CachePDFsAPFEL(Q0)
*
      write(6,*) "Cached evolution:"
      write(6,'(a5,2a12,a14,a10,3a12,a13,a14)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon",
     2         "e^-+e^+","mu^-+mu^+","tau^-+tau^+"
      do ilha=3,11
         write(6,'(es7.1,9es12.4)') 
     1         xlha(ilha),
     2         xPDFxQ(2,xlha(ilha),Q) - xPDFxQ(-2,xlha(ilha),Q),
     3         xPDFxQ(1,xlha(ilha),Q) - xPDFxQ(-1,xlha(ilha),Q),
     4         2d0*(xPDFxQ(-1,xlha(ilha),Q) + xPDFxQ(-2,xlha(ilha),Q)),
     5         xPDFxQ(4,xlha(ilha),Q) + xPDFxQ(-4,xlha(ilha),Q),
     6         xPDFxQ(0,xlha(ilha),Q),
     7         xPDFxQ(22,xlha(ilha),Q),
     8         xPDFxQ(11,xlha(ilha),Q)+xPDFxQ(-11,xlha(ilha),Q),
     9         xPDFxQ(13,xlha(ilha),Q)+xPDFxQ(-13,xlha(ilha),Q),
     1         xPDFxQ(15,xlha(ilha),Q)+xPDFxQ(-15,xlha(ilha),Q)
      enddo
      write(*,*) "  "
*
      write(6,*) "Cached evolution using the xPDFxQall function:"
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
      do ilha=3,11
         call xPDFxQall(xlha(ilha),Q,xf)
         write(6,'(es7.1,6es12.4)') 
     1         xlha(ilha),
     2         xf(2) - xf(-2),
     3         xf(1) - xf(-1),
     4         2d0 * ( xf(-1) + xf(-2) ),
     5         xf(4) + xf(-4),
     6         xf(0)
      enddo
      write(*,*) "  "
*
      end
