************************************************************************
*
*     ConvolutePDFsWithDISOperators.f:
*
*     This routine convolutes the final state PDFs with the DIS
*     operators.
*
************************************************************************
      subroutine ConvolutePDFsWithDISOperators
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/DISOperators.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/fph.h"
      include "../commons/StructureFunctions.h"
      include "../commons/TMC.h"
**
*     Internal Variables
*
      integer ipdf,ihq
      integer alpha,beta
      integer alphap,dgrid
      double precision t1,t2
      double precision fpy(-6:6,0:nint_max)
      double precision fev(0:13,0:nint_max)
*
      call cpu_time(t1)
*
      do igrid=1,ngrid
         do alpha=0,nin(igrid)
            do ipdf=-6,6
               fpy(ipdf,alpha) = fph(igrid,ipdf,alpha)
            enddo
         enddo
*
         call PDFphys2evQCD(fpy,fev)
*
         do alpha=0,nin(igrid)
            fev(0,alpha) = fgamma(igrid,alpha)
         enddo
*
         if(IsExt(igrid))then
            do ihq=3,7
               do alpha=0,nin(igrid)
                  F2(ihq,igrid,alpha) = 0d0
                  FL(ihq,igrid,alpha) = 0d0
                  F3(ihq,igrid,alpha) = 0d0
                  do beta=alpha,nin(igrid)
                     do ipdf=0,13
                        F2(ihq,igrid,alpha) = F2(ihq,igrid,alpha)
     1                       + OpF2(igrid,ihq,ipdf,alpha,beta)
     2                       * fev(ipdf,beta)
                        FL(ihq,igrid,alpha) = FL(ihq,igrid,alpha)
     1                       + OpFL(igrid,ihq,ipdf,alpha,beta)
     2                       * fev(ipdf,beta)
                        F3(ihq,igrid,alpha) = F3(ihq,igrid,alpha)
     1                       + OpF3(igrid,ihq,ipdf,alpha,beta)
     2                       * fev(ipdf,beta)
                     enddo
                  enddo
               enddo
            enddo
         else
            do ihq=3,7
               do alpha=0,nin(igrid)
                  F2(ihq,igrid,alpha) = 0d0
                  FL(ihq,igrid,alpha) = 0d0
                  F3(ihq,igrid,alpha) = 0d0
                  do beta=0,nin(igrid)-alpha
                     do ipdf=0,13
                        F2(ihq,igrid,alpha) = F2(ihq,igrid,alpha)
     1                       + OpF2(igrid,ihq,ipdf,0,beta)
     2                       * fev(ipdf,alpha+beta)
                        FL(ihq,igrid,alpha) = FL(ihq,igrid,alpha)
     1                       + OpFL(igrid,ihq,ipdf,0,beta)
     2                       * fev(ipdf,alpha+beta)
                        F3(ihq,igrid,alpha) = F3(ihq,igrid,alpha)
     1                       + OpF3(igrid,ihq,ipdf,0,beta)
     2                       * fev(ipdf,alpha+beta)
                     enddo
                  enddo
               enddo
            enddo
         endif
*
*     Terms needed for the target mass corrections
*
         if(TMC)then
            do ihq=3,7
               do alpha=0,nin(igrid)
                  I2(ihq,igrid,alpha) = 0d0
                  I3(ihq,igrid,alpha) = 0d0
                  do beta=alpha,nin(igrid)
                     I2(ihq,igrid,alpha) = I2(ihq,igrid,alpha)
     1                    + J_TMC(igrid,alpha,beta)
     2                    * F2(ihq,igrid,beta)
     3                    / xg(igrid,beta)**2 
                     I3(ihq,igrid,alpha) = I3(ihq,igrid,alpha)
     1                    + J_TMC(igrid,alpha,beta)
     2                    * F3(ihq,igrid,beta)
     3                    / xg(igrid,beta)**2 
                  enddo
               enddo
            enddo
         endif
      enddo
*
*     Fill joint structure functions operators
*
      do alpha=0,nin(0)+inter_degree(0)
*     Determine starting grid
         do dgrid=1,ngrid
            if(alpha.ge.TransitionPoint(dgrid).and.
     1         alpha.lt.TransitionPoint(dgrid+1))then
               goto 102
            endif
         enddo
         if(alpha.ge.nin(0)) dgrid = ngrid
 102     alphap  = alpha - TransitionPoint(dgrid)
         do ihq=3,7
            F2(ihq,0,alpha) = F2(ihq,dgrid,alphap)
            FL(ihq,0,alpha) = FL(ihq,dgrid,alphap)
            F3(ihq,0,alpha) = F3(ihq,dgrid,alphap)
            I2(ihq,0,alpha) = I2(ihq,dgrid,alphap)
            I3(ihq,0,alpha) = I3(ihq,dgrid,alphap)
         enddo
      enddo
*
      call cpu_time(t2)
*
c      write(6,"(a,a,f9.5,a)") " Convolution of the DIS operators",
c     1                        " with PDFs completed in",t2-t1," s"
c      write(6,*) " "
*
      return
      end 
