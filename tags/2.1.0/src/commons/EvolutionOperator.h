*     -*-fortran-*-
*
*     Total Evolution Operator to apply to any PDF set
*
      double precision PhQCD(0:ngrid_max,-7:6,-7:6,
     1                       0:nint_max,0:nint_max)
*
      common / EvolOpQCDAPFEL / PhQCD
