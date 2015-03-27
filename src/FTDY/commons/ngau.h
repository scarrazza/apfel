*     -*-fortran-*-
*
*     Parametrs for the Gauss integration in Drell-Yan
*
      integer ngau
      parameter(ngau=32)
*
      double precision y1(ngau),y2(ngau),w1(ngau),w2(ngau)
*
      common / yw1 / y1,w1
      common / yw2 / y2,w2
