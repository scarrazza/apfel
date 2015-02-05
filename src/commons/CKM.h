*     -*-fortran-*-

      double precision V_ud,V_us
      double precision V_cd,V_cs
      double precision V_ub,V_cb
      double precision V_td,V_ts,V_tb
*
      double precision V_ud2,V_us2
      double precision V_cd2,V_cs2
      double precision V_ub2,V_cb2
      double precision V_td2,V_ts2,V_tb2

      character*4 InCKM
*
      common / CKMMatrixAPFEL / V_ud,V_us,V_ub,
     1                          V_cd,V_cs,V_cb,
     2                          V_td,V_ts,V_tb,
     3                          InCKM

      common / CKM2MatrixAPFEL / V_ud2,V_us2,V_ub2,
     1                           V_cd2,V_cs2,V_cb2,
     2                           V_td2,V_ts2,V_tb2

c      parameter(V_ud = 0.97428d0, V_us = 0.22530d0)
c      parameter(V_cd = 0.22520d0, V_cs = 0.97345d0)
c      parameter(V_ub = 0.00347d0, V_cb = 0.041000d0)
c      parameter(V_td = 0.00862d0, V_ts = 0.04030d0, V_tb = 0.999152d0)

c      parameter(V_ud2 = 0.9492215184d0 , V_us2 = 5.0760090000d-2)
c      parameter(V_cd2 = 5.0715040000d-2, V_cs2 = 0.9476049025d0 )
c      parameter(V_ub2 = 1.2040900000d-5, V_cb2 = 1.6810000000d-3)
c      parameter(V_td2 = 7.4304399999d-5, V_ts2 = 1.6240900000d-3, 
c     1          V_tb2 = 0.998304719104d0)

