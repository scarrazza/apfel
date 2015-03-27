*     -*-fortran-*-
*
*     Kinematic of the data points for the FK tables
*
      integer ndata
      double precision q2dat(mxdata),ydat(mxdata)
      double precision x1dat(mxdata),x2dat(mxdata)
      character*15 obs(mxdata)
*
      common / kinematics / q2dat,ydat,x1dat,x2dat,obs,ndata
