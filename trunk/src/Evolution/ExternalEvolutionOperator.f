************************************************************************
*
*     ExternalEvolutionOperator.f:
*
*     This function returns the "external" evolution operator on the
*     interpolation grid. This is possible only if the internal grids are 
*     used of if one single external grid is provided.
*
************************************************************************
      function ExternalEvolutionOperator(i,j,x,beta)
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
      integer beta
      double precision x
**
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int_gen
**
*     Output Variables
*
      double precision ExternalEvolutionOperator
*
*     Check whether the evolution operator has been actually computed
*
      if(.not.EvolOp)then
         write(6,*) "The evolution operator computation is disabled."
         write(6,*) "This function cannot be used."
         write(6,*) "   "
         call exit(-10)
      endif
*
*     Check consistency of the input variables
*
      if(i.lt.-7.or.i.gt.6)then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid index, i =",i
         call exit(-10)
      endif
      if(j.lt.-7.or.j.gt.6)then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid index, j =",j
         call exit(-10)
      endif
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
      if(beta.lt.0.or.beta.gt.nin(0))then
         write(6,*) "In ExternalEvolutionOperator.f:"
         write(6,*) "Invalid index, beta =",beta
         call exit(-10)
      endif
*
*     interpolate
*
      n = inter_degree(0)
      ExternalEvolutionOperator = 0d0
      do alpha=0,nin(0)
         ExternalEvolutionOperator = ExternalEvolutionOperator 
     1                             + w_int_gen(n,alpha,x)
     2                             * PhQCD(0,i,j,alpha,beta)
      enddo
*
      return
      end
