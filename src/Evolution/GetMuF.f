************************************************************************
*
*     GetMuF.f:
*
*     GetMuF returns the final factorization scale.
*     GetMuF0 returns the initial factorization scale.
*
************************************************************************
      function GetMuF()
*
      implicit none
*
      include "../commons/scales.h"
*
*     Variables
*
      double precision GetMuF
*
      GetMuF = dsqrt(Q2fin)
*
      return
      end
*
************************************************************************
      function GetMuF0()
*
      implicit none
*
      include "../commons/scales.h"
*
*     Variables
*
      double precision GetMuF0
*
      GetMuF0 = dsqrt(Q2ini)
*
      return
      end
