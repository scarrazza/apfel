************************************************************************
*
*     ConvoluteEvolutionWithDISOperators.f:
*
*     This routine convolutes the evolution operators with the DIS
*     operators.
*
************************************************************************
      subroutine ConvoluteEvolutionWithDISOperators
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/DISOperators.h"
      include "../commons/EvolutionOperator.h"
      include "../commons/TMC.h"
**
*     Internal Variables
*
      integer i,j,ihq
      integer alpha,beta,gamma
      double precision t1,t2
*
      call cpu_time(t1)
*
      do ihq=3,7
         do i=0,13
            do alpha=0,nin(0)
               do beta=0,nin(0)
                  EvOpF2(ihq,i,alpha,beta) = 0d0
                  EvOpFL(ihq,i,alpha,beta) = 0d0
                  EvOpF3(ihq,i,alpha,beta) = 0d0
               enddo
            enddo
            do alpha=0,nin(0)
               do beta=alpha,nin(0)
                  do j=0,13
                     do gamma=alpha,beta
                        EvOpF2(ihq,i,alpha,beta) =
     1                       EvOpF2(ihq,i,alpha,beta)
     2                       + OpF2(0,ihq,j,alpha,gamma)
     3                       * Ev2EvQCD(0,j,i,gamma,beta)
                        EvOpFL(ihq,i,alpha,beta) =
     1                       EvOpFL(ihq,i,alpha,beta)
     2                       + OpFL(0,ihq,j,alpha,gamma)
     3                       * Ev2EvQCD(0,j,i,gamma,beta)
                        EvOpF3(ihq,i,alpha,beta) =
     1                       EvOpF3(ihq,i,alpha,beta)
     2                       + OpF3(0,ihq,j,alpha,gamma)
     3                       * Ev2EvQCD(0,j,i,gamma,beta)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
*
*     Now convolute the TMC operators if needed
*
      if(TMC)then
         do ihq=3,7
            do i=0,13
               do alpha=0,nin(0)
                  do beta=0,nin(0)
                     EvOpI2(ihq,i,alpha,beta) = 0d0
                     EvOpI3(ihq,i,alpha,beta) = 0d0
                  enddo
               enddo
               do alpha=0,nin(0)
                  do beta=alpha,nin(0)
                     do j=0,13
                        do gamma=alpha,beta
                           EvOpI2(ihq,i,alpha,beta) =
     1                          EvOpI2(ihq,i,alpha,beta)
     2                          + OpI2(0,ihq,j,alpha,gamma)
     3                          * Ev2EvQCD(0,j,i,gamma,beta)
                           EvOpI3(ihq,i,alpha,beta) =
     1                          EvOpI3(ihq,i,alpha,beta)
     2                          + OpI3(0,ihq,j,alpha,gamma)
     3                          * Ev2EvQCD(0,j,i,gamma,beta)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
*
      call cpu_time(t2)
*
c      write(6,"(a,a,f9.5,a)") " Convolution of the DIS operators with ",
c     1                        "the evolution operators completed in",
c     2                        t2-t1," s"
c      write(6,*) " "
*
      return
      end
