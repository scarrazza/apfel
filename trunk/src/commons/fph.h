*     -*-fortran-*-
*
*     PDFs in the physical basis at the final scale on the interpolation grid
*
      double precision fph(0:ngrid_max,-6:6,0:nint_max)
      double precision fgamma(0:ngrid_max,0:nint_max)

      double precision dfph(0:ngrid_max,-6:6,0:nint_max)
      double precision dfgamma(0:ngrid_max,0:nint_max)
*
      common / pdff / fph,fgamma
      common / dpdff / dfph,dfgamma
