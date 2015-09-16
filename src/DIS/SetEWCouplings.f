************************************************************************
*
*     SetEWCouplings.f:
*
*     This subroutine sets the values of the vector and axial couplings
*     of the up- and down-type quarks.
*
************************************************************************
      subroutine SetEWCouplings(vd,vu,ad,au)
*
      implicit none
*
      include "../commons/EWCouplings.h"
*
*     Variables
*
      double precision vd,vu,ad,au
*
      VectorD = vd
      VectorU = vu
      AxialD  = ad
      AxialU  = au
*
*     If they are all zero, tell the code to use the standard couplings
*
      ExtCoup = .true.
      if(vd.eq.0d0.and.vu.eq.0d0.and.
     1   ad.eq.0d0.and.au.eq.0d0) ExtCoup = .false.
*
      InEWCouplings = "done"
*
      return
      end
