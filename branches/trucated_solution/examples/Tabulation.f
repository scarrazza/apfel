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
      integer i,n
      double precision x,xmin,xmax,step
      double precision Q0,Q
      double precision AlphaQCD,AlphaQED
      double precision xPDF,xgamma
      double precision eps,epstrunc

      double precision SigAPz(1000),SigAPp(1000),SigAPm(1000)
      double precision gAPz(1000),gAPp(1000),gAPm(1000)
      double precision ValAPz(1000),ValAPp(1000),ValAPm(1000)
      double precision T3APz(1000),T3APp(1000),T3APm(1000)
      double precision V3APz(1000),V3APp(1000),V3APm(1000)
      double precision SigAPexact(1000),SigAPexpand(1000)
      double precision gAPexact(1000),gAPexpand(1000)
      double precision ValAPexact(1000),ValAPexpand(1000)
      double precision T3APexact(1000),T3APexpand(1000)
      double precision V3APexact(1000),V3APexpand(1000)
      double precision SigAPtrun
      double precision gAPtrun
      double precision ValAPtrun
      double precision T3APtrun
      double precision V3APtrun
      double precision SigNN
      double precision gNN
      double precision ValNN
      double precision T3NN
      double precision V3NN

      parameter(eps=1d-10)
*
      Q0 = dsqrt(2d0) - eps
      Q  = dsqrt(10000d0)
      epstrunc = 1d-7
*
      n = 100
      xmin = 1d-5
      xmax = 0.5d0
      step = dexp( ( dlog(xmax) - dlog(xmin) ) / dble( n - 1 ) )
*
      call SetEpsTrunc(0d0)
      call SetFastEvolution(.true.)
      call SetPerturbativeOrder(1)
      call SetAlphaQCDRef(0.118d0,91.2d0)
      call SetPoleMasses(dsqrt(2d0),4.75d0,175d0)
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("expandalpha")
      call InitializeAPFEL
      call EvolveAPFEL(Q0,Q)
      x = xmin
      do i=1,n
         SigAPz(i) = xPDF(1,x) + xPDF(-1,x)
     1        + xPDF(2,x) + xPDF(-2,x)
     2        + xPDF(3,x) + xPDF(-3,x)
     3        + xPDF(4,x) + xPDF(-4,x)
     4        + xPDF(5,x) + xPDF(-5,x)
     5        + xPDF(6,x) + xPDF(-6,x)
         ValAPz(i) = xPDF(1,x) - xPDF(-1,x)
     1        + xPDF(2,x) - xPDF(-2,x)
     2        + xPDF(3,x) - xPDF(-3,x)
     3        + xPDF(4,x) - xPDF(-4,x)
     4        + xPDF(5,x) - xPDF(-5,x)
     5        + xPDF(6,x) - xPDF(-6,x)
         T3APz(i)  = xPDF(2,x) + xPDF(-2,x)
     1        - xPDF(1,x) - xPDF(-1,x)
         V3APz(i)  = xPDF(2,x) - xPDF(-2,x)
     1        - xPDF(1,x) + xPDF(-1,x)
         gAPz(i) = xPDF(0,x)
         x = x * step
      enddo
      call CleanUp
*
      call SetEpsTrunc(epstrunc)
      call SetFastEvolution(.true.)
      call SetPerturbativeOrder(1)
      call SetAlphaQCDRef(0.118d0,91.2d0)
      call SetPoleMasses(dsqrt(2d0),4.75d0,175d0)
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("expandalpha")
      call InitializeAPFEL
      call EvolveAPFEL(Q0,Q)
      x = xmin
      do i=1,n
         SigAPp(i) = xPDF(1,x) + xPDF(-1,x)
     1        + xPDF(2,x) + xPDF(-2,x)
     2        + xPDF(3,x) + xPDF(-3,x)
     3        + xPDF(4,x) + xPDF(-4,x)
     4        + xPDF(5,x) + xPDF(-5,x)
     5        + xPDF(6,x) + xPDF(-6,x)
         ValAPp(i) = xPDF(1,x) - xPDF(-1,x)
     1        + xPDF(2,x) - xPDF(-2,x)
     2        + xPDF(3,x) - xPDF(-3,x)
     3        + xPDF(4,x) - xPDF(-4,x)
     4        + xPDF(5,x) - xPDF(-5,x)
     5        + xPDF(6,x) - xPDF(-6,x)
         T3APp(i)  = xPDF(2,x) + xPDF(-2,x)
     1        - xPDF(1,x) - xPDF(-1,x)
         V3APp(i)  = xPDF(2,x) - xPDF(-2,x)
     1        - xPDF(1,x) + xPDF(-1,x)
         gAPp(i) = xPDF(0,x)
         x = x * step
      enddo
      call CleanUp
*
      call SetEpsTrunc(-epstrunc)
      call SetFastEvolution(.true.)
      call SetPerturbativeOrder(1)
      call SetAlphaQCDRef(0.118d0,91.2d0)
      call SetPoleMasses(dsqrt(2d0),4.75d0,175d0)
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("expandalpha")
      call InitializeAPFEL
      call EvolveAPFEL(Q0,Q)
      x = xmin
      do i=1,n
         SigAPm(i) = xPDF(1,x) + xPDF(-1,x)
     1        + xPDF(2,x) + xPDF(-2,x)
     2        + xPDF(3,x) + xPDF(-3,x)
     3        + xPDF(4,x) + xPDF(-4,x)
     4        + xPDF(5,x) + xPDF(-5,x)
     5        + xPDF(6,x) + xPDF(-6,x)
         ValAPm(i) = xPDF(1,x) - xPDF(-1,x)
     1        + xPDF(2,x) - xPDF(-2,x)
     2        + xPDF(3,x) - xPDF(-3,x)
     3        + xPDF(4,x) - xPDF(-4,x)
     4        + xPDF(5,x) - xPDF(-5,x)
     5        + xPDF(6,x) - xPDF(-6,x)
         T3APm(i)  = xPDF(2,x) + xPDF(-2,x)
     1        - xPDF(1,x) - xPDF(-1,x)
         V3APm(i)  = xPDF(2,x) - xPDF(-2,x)
     1        - xPDF(1,x) + xPDF(-1,x)
         gAPm(i) = xPDF(0,x)
         x = x * step
      enddo
      call CleanUp
*
      call SetFastEvolution(.true.)
      call SetPerturbativeOrder(1)
      call SetAlphaQCDRef(0.118d0,91.2d0)
      call SetPoleMasses(dsqrt(2d0),4.75d0,175d0)
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
c      call SetAlphaEvolution("expanded")
c      call SetPDFEvolution("expandalpha")
      call InitializeAPFEL
      call EvolveAPFEL(Q0,Q)
      x = xmin
      do i=1,n
         SigAPexact(i) = xPDF(1,x) + xPDF(-1,x)
     1        + xPDF(2,x) + xPDF(-2,x)
     2        + xPDF(3,x) + xPDF(-3,x)
     3        + xPDF(4,x) + xPDF(-4,x)
     4        + xPDF(5,x) + xPDF(-5,x)
     5        + xPDF(6,x) + xPDF(-6,x)
         ValAPexact(i) = xPDF(1,x) - xPDF(-1,x)
     1        + xPDF(2,x) - xPDF(-2,x)
     2        + xPDF(3,x) - xPDF(-3,x)
     3        + xPDF(4,x) - xPDF(-4,x)
     4        + xPDF(5,x) - xPDF(-5,x)
     5        + xPDF(6,x) - xPDF(-6,x)
         T3APexact(i)  = xPDF(2,x) + xPDF(-2,x)
     1        - xPDF(1,x) - xPDF(-1,x)
         V3APexact(i)  = xPDF(2,x) - xPDF(-2,x)
     1        - xPDF(1,x) + xPDF(-1,x)
         gAPexact(i) = xPDF(0,x)
         x = x * step
      enddo
      call CleanUp
*
      call SetFastEvolution(.true.)
      call SetPerturbativeOrder(1)
      call SetAlphaQCDRef(0.118d0,91.2d0)
      call SetPoleMasses(dsqrt(2d0),4.75d0,175d0)
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
      call SetAlphaEvolution("expanded")
      call SetPDFEvolution("expandalpha")
      call InitializeAPFEL
      call EvolveAPFEL(Q0,Q)
      x = xmin
      do i=1,n
         SigAPexpand(i) = xPDF(1,x) + xPDF(-1,x)
     1        + xPDF(2,x) + xPDF(-2,x)
     2        + xPDF(3,x) + xPDF(-3,x)
     3        + xPDF(4,x) + xPDF(-4,x)
     4        + xPDF(5,x) + xPDF(-5,x)
     5        + xPDF(6,x) + xPDF(-6,x)
         ValAPexpand(i) = xPDF(1,x) - xPDF(-1,x)
     1        + xPDF(2,x) - xPDF(-2,x)
     2        + xPDF(3,x) - xPDF(-3,x)
     3        + xPDF(4,x) - xPDF(-4,x)
     4        + xPDF(5,x) - xPDF(-5,x)
     5        + xPDF(6,x) - xPDF(-6,x)
         T3APexpand(i)  = xPDF(2,x) + xPDF(-2,x)
     1        - xPDF(1,x) - xPDF(-1,x)
         V3APexpand(i)  = xPDF(2,x) - xPDF(-2,x)
     1        - xPDF(1,x) + xPDF(-1,x)
         gAPexpand(i) = xPDF(0,x)
         x = x * step
      enddo
      call CleanUp
*
      call SetPDFset("NNPDF23_nlo_as_0118.LHgrid")
      call InitializeAPFEL
      call EvolveAPFEL(Q,Q)
      x = xmin
      do i=1,n
         SigNN = xPDF(1,x) + xPDF(-1,x)
     1         + xPDF(2,x) + xPDF(-2,x)
     2         + xPDF(3,x) + xPDF(-3,x)
     3         + xPDF(4,x) + xPDF(-4,x)
     4         + xPDF(5,x) + xPDF(-5,x)
     5         + xPDF(6,x) + xPDF(-6,x)
         ValNN = xPDF(1,x) - xPDF(-1,x)
     1         + xPDF(2,x) - xPDF(-2,x)
     2         + xPDF(3,x) - xPDF(-3,x)
     3         + xPDF(4,x) - xPDF(-4,x)
     4         + xPDF(5,x) - xPDF(-5,x)
     5         + xPDF(6,x) - xPDF(-6,x)
         T3NN  = xPDF(2,x) + xPDF(-2,x)
     1         - xPDF(1,x) - xPDF(-1,x)
         V3NN = xPDF(2,x) - xPDF(-2,x)
     1        - xPDF(1,x) + xPDF(-1,x)
         gNN = xPDF(0,x)
         SigAPtrun = SigAPz(i) + ( SigAPp(i) - SigAPm(i) ) /2d0/epstrunc
         ValAPtrun = ValAPz(i) + ( ValAPp(i) - ValAPm(i) ) /2d0/epstrunc
         T3APtrun = T3APz(i) + ( T3APp(i) - T3APm(i) ) /2d0/epstrunc
         V3APtrun = V3APz(i) + ( V3APp(i) - V3APm(i) ) /2d0/epstrunc
         gAPtrun = gAPz(i) + ( gAPp(i) - gAPm(i) ) /2d0/epstrunc

c         write(6,*) x,SigAPtrun/SigNN,ValAPtrun/ValNN,T3APtrun/T3NN,
c     1                V3APtrun/V3NN,gAPtrun/gNN
         write(6,*) x,gAPexact(i)/gNN,gAPexpand(i)/gNN,gAPtrun/gNN
c         write(6,*) x,SigAPexact(i)/SigNN,SigAPexpand(i)/SigNN,
c     1                SigAPtrun/SigNN
         x = x * step
      enddo
      call CleanUp





      end
