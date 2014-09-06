*     -*-fortran-*-
*
*     Parameter of the ODE algorithm
*
*     DO NOT CHANGE!
*
      double precision safety
      double precision pgrow
      double precision pshrnk
      double precision errcon

      parameter(safety=0.9d0)
      parameter(pgrow=-0.2d0)
      parameter(pshrnk=-0.25d0)
      parameter(errcon=1.89d-4)
*
      double precision A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,
     1B51,B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,
     2DC6
*
      parameter(A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     1B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     2B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     3B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     4C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     5DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     6DC6=C6-.25)
