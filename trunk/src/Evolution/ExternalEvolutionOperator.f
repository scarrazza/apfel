************************************************************************
*
*     ExternalEvolutionOperator.f:
*
*     This function returns the "external" evolution operator on the
*     interpolation grid. This is possible only if the internal grids are 
*     used of if one single external grid is provided.
*
************************************************************************
      function ExternalEvolutionOperator(i,j,alpha,beta)
*
      implicit none
*
      include "../commons/EvolOp.h"
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
      double precision ExternalEvolutionOperator
*
      if(.not.EvolOp)then
         write(6,*) "The evolution operator computation is disabled."
         write(6,*) "This function cannot be used."
         write(6,*) "   "
         call exit(-10)
      endif
*
      ExternalEvolutionOperator = PhQCD(0,i,j,alpha,beta)
*
      return
      end
