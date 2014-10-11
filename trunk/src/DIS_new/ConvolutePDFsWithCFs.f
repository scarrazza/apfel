************************************************************************
*
*     ConvolutePDFsWithCFs.f:
*
*     This routine combines the evolved PDFs withe DIS coefficient 
*     functions computed on the grid.
*
************************************************************************
      subroutine ConvolutePDFsWithCFs(Q)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/fph.h"
      include "../commons/integralsDIS.h"
      include "../commons/m2th.h"
      include "../commons/StructureFunctions.h"
      include "../commons/TargetDIS.h"
**
*     Internal Variables
*
      double precision Q
**
*     Internal Variables
*
      integer jgrid
      integer nf
      integer i
      integer alpha,beta
      integer pt
      double precision Q2
      double precision singlet
      double precision F2t,FLt,F3t
      double precision as,a_QCD
      double precision bq(6),dq(6)
      double precision fup,fub,fdw,fdb
*
      Q2 = Q * Q
*
*     Find number of active flavours
*
      if(Q2.ge.m2th(6))then
         nf = 6
      elseif(Q2.ge.m2th(5))then
         nf = 5
      elseif(Q2.ge.m2th(4))then
         nf = 4
      else
         nf = 3
      endif
*
*     Compute needed couplings
*
      as = a_QCD(Q2)
      call ComputeChargesDIS(Q2,bq,dq)
*
*     F2
*
      do jgrid=1,ngrid
         do alpha=0,nin(jgrid)
*
*     Rearrange PDFs in case the target is not a proton.
*     (This assumes isospin symmetry)
*
            if(TargetDIS(1:7).eq."neutron")then
*     Exchage up and down
               do beta=0,nin(jgrid)
                  fdw = fph(jgrid,1,beta)
                  fdb = fph(jgrid,-1,beta)
                  fup = fph(jgrid,2,beta)
                  fub = fph(jgrid,-2,beta)
*
                  fph(jgrid,1,beta)  = fup
                  fph(jgrid,-1,beta) = fub
                  fph(jgrid,2,beta)  = fdw
                  fph(jgrid,-2,beta) = fdb
               enddo
            elseif(TargetDIS.eq."isoscalar")then
*     Average up and down
               do beta=0,nin(jgrid)
                  fdw = fph(jgrid,1,beta)
                  fdb = fph(jgrid,-1,beta)
                  fup = fph(jgrid,2,beta)
                  fub = fph(jgrid,-2,beta)
*
                  fph(jgrid,1,beta)  = ( fup + fdw ) / 2d0
                  fph(jgrid,-1,beta) = ( fub + fdb ) / 2d0
                  fph(jgrid,2,beta)  = fph(jgrid,1,beta)
                  fph(jgrid,-2,beta) = fph(jgrid,-1,beta)
               enddo
            endif
*
            do i=1,7
               F2(i,jgrid,alpha) = 0d0
            enddo
            do beta=0,nin(jgrid)-alpha
               singlet = 0d0
               do i=1,nf
                  singlet = singlet 
     1                    + fph(jgrid,i,alpha+beta)
     2                    + fph(jgrid,-i,alpha+beta)
               enddo
*     F2 flavour by flavour
               do i=1,nf
                  do pt=0,ipt
                     F2t = bq(i) 
     1                   * ( SC2(jgrid,nf,1,pt,0,beta) ! Gluon
     2                   * fph(jgrid,0,alpha+beta)
     3                   + SC2(jgrid,nf,2,pt,0,beta)   ! Singlet
     4                   * singlet
     5                   + SC2(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     6                   * ( fph(jgrid,i,alpha+beta) 
     7                   + fph(jgrid,-i,alpha+beta) ) )
*
                     F2(i,jgrid,alpha) = F2(i,jgrid,alpha)
     1                                 + as**pt * F2t
                  enddo
               enddo
            enddo
*     F2 total
            do i=1,nf
               F2(7,jgrid,alpha) = F2(7,jgrid,alpha) + F2(i,jgrid,alpha)
            enddo
*
*     FL
*
            do i=1,7
               FL(i,jgrid,alpha) = 0d0
            enddo
            do beta=0,nin(jgrid)-alpha
               singlet = 0d0
               do i=1,nf
                  singlet = singlet 
     1                    + fph(jgrid,i,alpha+beta)
     2                    + fph(jgrid,-i,alpha+beta)
               enddo
*     FL flavour by flavour
               do i=1,nf
                  do pt=0,ipt
                     FLt = bq(i) 
     1                   * ( SCL(jgrid,nf,1,pt,0,beta) ! Gluon
     2                   * fph(jgrid,0,alpha+beta)
     3                   + SCL(jgrid,nf,2,pt,0,beta)   ! Singlet
     4                   * singlet
     5                   + SCL(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     6                   * ( fph(jgrid,i,alpha+beta) 
     7                   + fph(jgrid,-i,alpha+beta) ) )
*
                     FL(i,jgrid,alpha) = FL(i,jgrid,alpha)
     1                                 + as**pt * FLt
                  enddo
               enddo
            enddo
*     FL total
            do i=1,nf
               FL(7,jgrid,alpha) = FL(7,jgrid,alpha) + FL(i,jgrid,alpha)
            enddo
*
*     F3
*
            do i=1,7
               F3(i,jgrid,alpha) = 0d0
            enddo
            do beta=0,nin(jgrid)-alpha
*     F3 flavour by flavour
               do i=1,nf
                  do pt=0,ipt
                     F3t = dq(i) 
     1                   * SC3(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     2                   * ( fph(jgrid,i,alpha+beta) 
     3                   - fph(jgrid,-i,alpha+beta) )
*
                     F3(i,jgrid,alpha) = F3(i,jgrid,alpha)
     1                                 + as**pt * F3t
                  enddo
               enddo
            enddo
*     F3 total
            do i=1,nf
               F3(7,jgrid,alpha) = F3(7,jgrid,alpha) + F3(i,jgrid,alpha)
            enddo
         enddo
      enddo
*
      return
      end
