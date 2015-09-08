************************************************************************
*
*     DISObservables.f:
*
*     Test program for the DIS observables.
*
************************************************************************
      program DISObservables
*
      implicit none
**
*     Variables
*
      integer ilha
      double precision Q0,Q
      double precision Q02,Q2
      double precision AlphaQCD,AlphaQED
      double precision F2light,F2charm,F2bottom,F2total
      double precision FLlight,FLcharm,FLbottom,FLtotal
      double precision F3light,F3charm,F3bottom,F3total
      double precision eps
      double precision xlha(11)

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /
*
*     Settings
*
c      call SetMassScheme("ZM-VFNS")
      call SetProcessDIS("NC")
c      call SetEWCouplings(1d0,1d0,1d0,1d0)
c      call SetPropagatorCorrection(0.1d0)
c      call SetPolarizationDIS(0d0)
c      call SetProjectileDIS("electron")
c      call SetTargetDIS("proton")
c      call EnableTargetMassCorrections(.false.)
c      call EnableDampingFONLL(.true.)
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
      call InitializeAPFEL_DIS
*
*     Evolve PDFs on the grids
*
      write(6,*) "Enter initial and final scale in GeV^2"
      read(5,*) Q02,Q2
*
      Q0 = dsqrt(Q02) - eps
      Q  = dsqrt(Q2)
      call ComputeStructureFunctionsAPFEL(Q0,Q)
*
*     Tabulate PDFs for the LHA x values
*
      write(6,*) "alpha_QCD(mu2F) =",AlphaQCD(Q)
      write(6,*) "alpha_QED(mu2F) =",AlphaQED(Q)
      write(6,*) "  "
*
      write(6,'(a5,4a12)') "x","F2light","F2charm","F2bottom","F2total"
      do ilha=3,11
         write(6,'(es7.1,9es12.4)') 
     1         xlha(ilha),
     2         F2light(xlha(ilha)),
     3         F2charm(xlha(ilha)),
     4         F2bottom(xlha(ilha)),
     5         F2total(xlha(ilha))
      enddo
      write(*,*) "  "
*
      write(6,'(a5,4a12)') "x","FLlight","FLcharm","FLbottom","FLtotal"
      do ilha=3,11
         write(6,'(es7.1,9es12.4)') 
     1         xlha(ilha),
     2         FLlight(xlha(ilha)),
     3         FLcharm(xlha(ilha)),
     4         FLbottom(xlha(ilha)),
     5         FLtotal(xlha(ilha))
      enddo
      write(*,*) "  "
*
      write(6,'(a5,4a12)') "x","F3light","F3charm","F3bottom","F3total"
      do ilha=3,11
         write(6,'(es7.1,9es12.4)') 
     1         xlha(ilha),
     2         F3light(xlha(ilha)),
     3         F3charm(xlha(ilha)),
     4         F3bottom(xlha(ilha)),
     5         F3total(xlha(ilha))
      enddo
      write(*,*) "  "
*
      end
