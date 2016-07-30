************************************************************************
*
*     Sum rules program
*
************************************************************************
      program SumRules
*
      implicit none
*
      integer i
      double precision Q0,Q
      double precision Q02,Q2
      double precision NPDF,Ngamma
      double precision momsr,uvsr,dvsr,svsr
      double precision eps
      parameter(eps=1d-10)
      character*1 answer
*
c      call SetFFNS(3)
      call SetPerturbativeOrder(1)
      call SetNumberOfGrids(3)
      call SetGridParameters(1,130,3,1d-9)
      call SetGridParameters(2,60,5,1d-1)
      call SetGridParameters(3,20,5,8d-1)
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
      momsr = 0d0
      do i=-6,6
         momsr = momsr + NPDF(i,2)
      enddo
      momsr = momsr + Ngamma(2)
*
      uvsr = NPDF(2,1) - NPDF(-2,1)
      dvsr = NPDF(1,1) - NPDF(-1,1)
      svsr = NPDF(3,1) - NPDF(-3,1)
*
      write(6,*) "Sum rules at Q =",Q," GeV:"
      write(6,*) "  "
      write(6,*) "- Momentum sum rule        =",momsr
      write(6,*) "- Up valence sum rule      =",uvsr
      write(6,*) "- Down valence sum rule    =",dvsr
      write(6,*) "- Strange valence sum rule =",svsr
      write(6,*) "  "
*
      write(6,*) "Do you want to evolve PDFs again? [y/n]"
      read(5,*) answer
      if(answer.eq."y") goto 101
*
      end
