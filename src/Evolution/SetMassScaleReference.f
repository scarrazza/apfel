************************************************************************
*
*     SetMassScaleReference.f:
*
*     This subroutine sets the reference scales at which heavy quark
*     masses are given. Inactive if the pole masses are used.
*
************************************************************************
      subroutine SetMassScaleReference(Qc,Qb,Qt)
*
      implicit none
*
      include "../commons/m2th.h"
*
*     Variables
*
      double precision Qc,Qb,Qt
*
      Q2th(4)   = Qc * Qc
      Q2th(5)   = Qb * Qb
      Q2th(6)   = Qt * Qt
      InMassRef = "done"
*
      return
      end
