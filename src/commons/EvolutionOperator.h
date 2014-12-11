*     -*-fortran-*-
*
*     Total Evolution Operator to apply to any PDF set
*
      double precision
     1     Ph2PhQCD(0:ngrid_max,-7:6,-7:6,0:nint_max,0:nint_max),
     2     Ev2PhQCD(0:ngrid_max,-7:6,0:13,0:nint_max,0:nint_max),
     3     Ev2EvQCD(0:ngrid_max,0:13,0:13,0:nint_max,0:nint_max)
*
      common / EvolOpQCDAPFEL / Ph2PhQCD,Ev2PhQCD,Ev2EvQCD
