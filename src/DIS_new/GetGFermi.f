************************************************************************
*
*     GetGFermi.f:
*
*     This function returns the Fermi constant.
*
************************************************************************
      function GetGFermi()
*
      implicit none
*
      include "../commons/GFermi.h"
*
*     Variables
*
      double precision GetGFermi
*
      GetGFermi = GFermi
*
      return
      end
