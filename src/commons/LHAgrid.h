*     -*-fortran-*-
*
*     Parameters of the LHAPDF evolution grid     
*
      integer nxLHA,nxmLHA,nq2LHA
      double precision q2minLHA,q2maxLHA
      double precision xminLHA,xmLHA,xmaxLHA
      double precision Lambda2
*
      parameter(nxLHA    = 100)
      parameter(nxmLHA   = 50)
      parameter(nq2LHA   = 50)
      parameter(q2minLHA = 1d0)
      parameter(q2maxLHA = 1d10)
      parameter(xminLHA  = 1d-9)
      parameter(xmLHA    = 1d-1)
      parameter(xmaxLHA  = 1d0)
      parameter(Lambda2  = 0.0625d0)
