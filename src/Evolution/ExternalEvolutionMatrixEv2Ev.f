************************************************************************
*
*     ExternalEvolutionMatrixEv2Ev.f:
*
*     This function returns the "external" evolution matrix on the
*     interpolation grid. This is possible only if the internal grids are 
*     used or if one single external grid is provided.
*
************************************************************************
      function ExternalEvolutionMatrixEv2Ev(i,j,alpha,beta)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/EvolutionOperator.h"
**
*     Input Variables
*
      integer i,j
      integer alpha,beta
**
*     Output Variables
*
      double precision ExternalEvolutionMatrixEv2Ev
*
      ExternalEvolutionMatrixEv2Ev = Ev2EvQCD(0,i,j,alpha,beta)
*
      return
      end
