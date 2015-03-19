*     -*-fortran-*-
*
*     PDFs in the physical basis at the initial scale on the interpolation grid
*
      double precision f0ph(-6:6,0:nint_max)
      double precision f0lep(-3:3,0:nint_max)
*
      common / pdf0APFEL / f0ph,f0lep
