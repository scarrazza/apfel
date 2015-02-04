************************************************************************
*
*     ExternalDISOperator.f:
*
*     This function returns the "external" DIS operators on the
*     interpolation grid. This is possible only if the internal grids are 
*     used or if one single external grid is provided.
*
************************************************************************
      function ExternalDISOperator(SF,ihq,i,x,beta)
*
      implicit none
*
      include "../commons/EvolOp.h"
      include "../commons/grid.h"
      include "../commons/DISOperators.h"
**
*     Input Variables
*
      integer ihq,i
      double precision x
      integer beta
      character*2 SF
**
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int_gen
**
*     Output Variables
*
      double precision ExternalDISOperator
*
*     Check whether the evolution operator has been actually computed
*
      if(.not.EvolOp)then
         write(6,*) "The evolution operator computation is disabled."
         write(6,*) "The 'ExternalDISOperator' function cannot be used."
         write(6,*) "   "
         call exit(-10)
      endif
*
*     Check consistency of the input variables
*
      if(SF.ne."F2".and.
     1   SF.ne."FL".and.
     2   SF.ne."F3")then
         write(6,*) "In ExternalDISOperator.f:"
         write(6,*) "Invalid Structure Function, SF = ",SF
         write(6,*) "  "
         write(6,*) "The options are:"
         write(6,*) "- 'F2'"
         write(6,*) "- 'FL'"
         write(6,*) "- 'F3'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(ihq.lt.3.or.ihq.gt.7)then
         write(6,*) "In ExternalDISOperator.f:"
         write(6,*) "Invalid HQ index, ihq =",ihq
         write(6,*) "  "
         write(6,*) "Valid range: ihq in [3:7]"
         call exit(-10)
      endif
      if(i.lt.0.or.i.gt.13)then
         write(6,*) "In ExternalDISOperator.f:"
         write(6,*) "Invalid index, i =",i
         call exit(-10)
      endif
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In ExternalDISOperator.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
      if(beta.lt.0.or.beta.gt.nin(0))then
         write(6,*) "In ExternalDISOperator.f:"
         write(6,*) "Invalid index, beta =",beta
         call exit(-10)
      endif
*
*     interpolate
*
      n = inter_degree(0)
      ExternalDISOperator = 0d0
      if(SF.eq."F2")then
         do alpha=0,nin(0)
            ExternalDISOperator = ExternalDISOperator 
     1                          + w_int_gen(n,alpha,x)
     2                          * EvOpF2(ihq,i,alpha,beta)
         enddo
      elseif(SF.eq."FL")then
         do alpha=0,nin(0)
            ExternalDISOperator = ExternalDISOperator 
     1                          + w_int_gen(n,alpha,x)
     2                          * EvOpFL(ihq,i,alpha,beta)
         enddo
      elseif(SF.eq."F3")then
         do alpha=0,nin(0)
            ExternalDISOperator = ExternalDISOperator 
     1                          + w_int_gen(n,alpha,x)
     2                          * EvOpF3(ihq,i,alpha,beta)
         enddo
      endif
*
      return
      end
