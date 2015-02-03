************************************************************************
*
*     JoinDISOperators.f:
*
*     This routine joins the subgrids used in the calculation in one
*     single grid.
*
************************************************************************
      subroutine JoinDISOperators
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/DISOperators.h"
**
*     Internal Variables
*
      integer ihq
      integer i
      integer alpha,beta,alphap,betap
      integer jgrid,dgrid,istart,density,offset
      double precision eps
      parameter(eps=1d-14)
*     
*     Set DIS Operators to zero
*     
      do i=0,13
         do ihq=3,7
            do alpha=0,nin(0)
               do beta=0,nin(0)
                  OpF2(0,ihq,i,alpha,beta) = 0d0
                  OpFL(0,ihq,i,alpha,beta) = 0d0
                  OpF3(0,ihq,i,alpha,beta) = 0d0
               enddo
            enddo
         enddo
      enddo
*
*     Fill DIS Operators
*     
      do alpha=0,nin(0)
*     Determine starting grid
         do dgrid=1,ngrid
            if(alpha.ge.TransitionPoint(dgrid).and.
     1         alpha.lt.TransitionPoint(dgrid+1))then
               goto 102
            endif
         enddo
         if(alpha.eq.nin(0)) dgrid = ngrid
 102     density = 1
         istart  = alpha
         alphap  = alpha - TransitionPoint(dgrid)
         betap   = alphap
         do jgrid=dgrid,ngrid
            do beta=istart,TransitionPoint(jgrid+1),density
               do ihq=3,7
                  do i=0,13
                     if(IsExt(dgrid))then
                        OpF2(0,ihq,i,alpha,beta) = 
     1                       OpF2(dgrid,ihq,i,alphap,betap)
                        OpFL(0,ihq,i,alpha,beta) = 
     1                       OpFL(dgrid,ihq,i,alphap,betap)
                        OpF3(0,ihq,i,alpha,beta) = 
     1                       OpF3(dgrid,ihq,i,alphap,betap)
                     else
                        OpF2(0,ihq,i,alpha,beta) = 
     1                       OpF2(dgrid,ihq,i,0,betap-alphap)
                        OpFL(0,ihq,i,alpha,beta) = 
     1                       OpFL(dgrid,ihq,i,0,betap-alphap)
                        OpF3(0,ihq,i,alpha,beta) = 
     1                       OpF3(dgrid,ihq,i,0,betap-alphap)
                     endif
*
                     if(abs(OpF2(0,ihq,i,alpha,beta)).lt.eps)
     1                    OpF2(0,ihq,i,alpha,beta) = 0d0 
                     if(abs(OpFL(0,ihq,i,alpha,beta)).lt.eps)
     1                    OpFL(0,ihq,i,alpha,beta) = 0d0 
                     if(abs(OpF3(0,ihq,i,alpha,beta)).lt.eps)
     1                    OpF3(0,ihq,i,alpha,beta) = 0d0 
                  enddo
               enddo     
               betap = betap + 1
            enddo
            offset  = ( density 
     1              - mod(TransitionPoint(jgrid+1)-istart,density) )
     2              * DensityFactor(jgrid+1)
            istart  = offset + TransitionPoint(jgrid+1)
            density = density * DensityFactor(jgrid+1)
         enddo
      enddo
*
      return
      end
