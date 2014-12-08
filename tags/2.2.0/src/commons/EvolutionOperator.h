*     -*-fortran-*-
*
*     Total Evolution Operator to apply to any PDF set
*
      double precision
     1     PhQCD(0:ngrid_max,-7:6,-7:6,0:nint_max,0:nint_max)
*
      common / EvolOpQCDAPFEL / PhQCD
