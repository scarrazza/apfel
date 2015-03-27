************************************************************************
*
*     wcDYlo.f:
*
*     Routine that returns the LO contribution to the double 
*     differential DY distributions.
*
************************************************************************
      subroutine cDY_LO(xi,zi,x0,cDY_NS_LO)
*
      implicit none
**
*     Input Variable
*
      integer xi,zi
      double precision x0(2)
**
*     Internal Variable
*
      double precision elin
**
*     Output Variable
*
      double precision cDY_NS_LO(0:1,0:1)
*
*     Linear interpolation
*
      cDY_NS_LO(0,0) = elin(zi,x0(2))   * elin(xi,x0(1))   ! \a=xi,\b=zi
      cDY_NS_LO(1,0) = elin(zi+1,x0(2)) * elin(xi,x0(1))   ! \a=xi,\b=zi+1
      cDY_NS_LO(0,1) = elin(zi,x0(2))   * elin(xi+1,x0(1)) ! \a=xi+1,\b=zi
      cDY_NS_LO(1,1) = elin(zi+1,x0(2)) * elin(xi+1,x0(1)) ! \a=xi+1,\b=zi+1
*
      return
      end
