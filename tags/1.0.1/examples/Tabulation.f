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
      character*1 answer

      parameter(eps=1d-10)
      data xlha / 1d-7, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2,
     1            1d-1, 3d-1, 5d-1, 7d-1, 9d-1 /            
*
*     Some examples ...
*
c      call SetFFNS(3)
c      call SetTheory("QECDS")
c      call SetPerturbativeOrder(2)
c      call SetPDFSet("NNPDF23_nlo_as_0119_qed.LHgrid")
c      call SetPDFSet("MRST2004qed.LHgrid")
*
*     Initializes integrals on the grids
*
      call InitializeAPFEL
*
*     Evolve PDFs on the grids
*
 101  write(6,*) "Enter initial and final scale in GeV^2"
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
      write(6,'(a5,2a12,a14,a10,2a12)') "x",
     1         "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon","photon"
*
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
      write(6,*) "Do you want evolve PDFs again? [y/n]"
      read(5,*) answer
      if(answer.eq."y") goto 101
*
      end
