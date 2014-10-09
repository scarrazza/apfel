*     -*-fortran-*-

      integer n_ref_qcd
      double precision lambda_ref_qcd
      double precision LambdaQCD(3:6)
      character*4 InLambdaQCD
*
      common / LambdaQCDAPFEL / lambda_ref_qcd,LambdaQCD,n_ref_qcd,
     1                          InLambdaQCD
