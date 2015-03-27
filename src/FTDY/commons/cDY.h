*     -*-fortran-*-
*
*     Hard cross sections fo the Fixed-Target Drell-Yan process
*
      integer ixp1DY(mxdata),ixp2DY(mxdata)
      real cDY_NS(0:1,mxgridsizeDY,mxgridsizeDY,mxdata)
      real cDY_QG(0:1,mxgridsizeDY,mxgridsizeDY,mxdata)
      real cDY_GQ(0:1,mxgridsizeDY,mxgridsizeDY,mxdata)
*
      common / ccDYNLO / cDY_NS,cDY_QG,cDY_GQ,ixp1DY,ixp2DY
