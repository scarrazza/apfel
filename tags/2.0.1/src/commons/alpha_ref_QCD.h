*     -*-fortran-*-

      double precision alpha_ref_qcd,q2_ref_qcd
      double precision LambdaQCD(3:6)
      character*4 InAlpQCD

      data LambdaQCD / 0.250d0, 0.250d0, 0.250d0, 0.250d0 /
*
      common / coupQCD / alpha_ref_qcd,q2_ref_qcd,InAlpQCD
c      common / lambda / LambdaQCD
