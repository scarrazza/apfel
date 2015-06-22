*     -*-fortran-*-
*
*     Parameters of the LHAPDF evolution grid     
*
      integer nxmax,nq2max
      integer nxLHA,nxmLHA,nq2LHA
      double precision q2minLHA,q2maxLHA
      double precision xminLHA,xmLHA,xmaxLHA
      double precision Lambda2
      character*4 InLHgrid
*
c      parameter(nxLHA    = 100)
c      parameter(nxmLHA   = 50)
c      parameter(nq2LHA   = 50)
c      parameter(q2minLHA = 1d0)
c      parameter(q2maxLHA = 1d10)
c      parameter(xminLHA  = 1d-9)
c      parameter(xmLHA    = 1d-1)
c      parameter(xmaxLHA  = 1d0)
      parameter(Lambda2  = 0.0625d0)
      parameter(nxmax    = 300)
      parameter(nq2max   = 200)
*
      common / LHgridParamAPFEL / xminLHA,xmLHA,xmaxLHA,q2minLHA,
     1                            q2maxLHA,nxLHA,nxmLHA,nq2LHA,
     2                            InLHgrid
