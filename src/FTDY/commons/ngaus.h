*     -*-fortran-*-
*
*     Parametrs for the Gauss integration in Drell-Yan
*
      integer ngaus
      parameter(ngaus=128)
*
      double precision y1s(ngaus),y2s(ngaus),w1s(ngaus),w2s(ngaus)
*
      common / yw1s / y1s,w1s
      common / yw2s / y2s,w2s
