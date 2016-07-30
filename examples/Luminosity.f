************************************************************************
*
*     Luminosity program
*
************************************************************************
      program Luminosity
*
      implicit none
*
      integer i, NMX
      double precision Q0,Q,Qmin,Qmax,S
      double precision LUMI
      double precision eps
      parameter(eps=1d-10)
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
      S = 8d3**2
      Q0 = dsqrt(2d0)
      Qmin = 10
      Qmax = 6d3
      NMX  = 30
*
      write(6,*) "Luminosities at sqrt(S) =",dsqrt(S)," GeV:"
      write(6,*) "  "

      write(6,*) "- gg luminosity "
      
      do i=1,30
         Q = Qmin * (Qmax/Qmin)**(float(i-1)/float(NMX-1))
         call EvolveAPFEL(Q0,Q)
         write(6,*) "MX = ",Q,"LUMI = ",LUMI(0,0,S)
      enddo

      write(6,*) "  "
*
      end
