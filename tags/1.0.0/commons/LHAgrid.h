*     -*-fortran-*-
*
*     Parameters of the LHAPDF evolution grid     
*
      integer nxLHA,nxmLHA,nq2LHA
      double precision q2minLHA,q2maxLHA
      double precision xminLHA,xmLHA,xmaxLHA
*
      parameter(nxLHA    = 100)
      parameter(nxmLHA   = 50)
      parameter(nq2LHA   = 50)
      parameter(q2minLHA = 1d0)
      parameter(q2maxLHA = 1d8)
      parameter(xminLHA  = 1d-7)
      parameter(xmLHA    = 1d-1)
      parameter(xmaxLHA  = 1d0)
