*     -*-fortran-*-
*
*     PDFs in the physical basis at the final scale on the interpolation grid
*
      double precision fph(ngrid_max,-6:6,0:nint_max)
      double precision fgamma(ngrid_max,0:nint_max)
*
      common / pdff / fph,fgamma
