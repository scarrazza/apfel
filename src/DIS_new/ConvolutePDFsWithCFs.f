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
      include "../commons/ProcessDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/CKM.h"
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
      integer ipr
      double precision Q2
      double precision singlet
      double precision F2t,FLt,F3t
      double precision as,a_QCD
      double precision bq(6),dq(6)
      double precision fup,fub,fdw,fdb
      double precision frac
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
*     Compute alphas
*
      as = a_QCD(Q2)
*
*     Electromagnetic and Neutral current structure functions
*
      if(ProcessDIS.eq."EM".or.ProcessDIS.eq."NC")then
*
*     Compute needed couplings
*
         call ComputeChargesDIS(Q2,bq,dq)
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
               elseif(TargetDIS(1:4).eq."iron")then
*     Reweight up and down by frac
                  frac = 23.403d0 / 49.618d0
                  do beta=0,nin(jgrid)
                     fdw = fph(jgrid,1,beta)
                     fdb = fph(jgrid,-1,beta)
                     fup = fph(jgrid,2,beta)
                     fub = fph(jgrid,-2,beta)
*
                     fph(jgrid,1,beta)  = frac * fup
     1                                  + ( 1d0 - frac ) * fdw 
                     fph(jgrid,-1,beta) = frac * fub
     1                                  + ( 1d0 - frac ) * fdb
                     fph(jgrid,2,beta)  = fph(jgrid,1,beta)
                     fph(jgrid,-2,beta) = fph(jgrid,-1,beta)
                  enddo
               endif
*
*     F2
*
               do i=1,7
                  F2(i,jgrid,alpha) = 0d0
               enddo
               do beta=0,nin(jgrid)-alpha
                  singlet = 0d0
                  do i=1,nf
                     singlet = singlet 
     1                       + fph(jgrid,i,alpha+beta)
     2                       + fph(jgrid,-i,alpha+beta)
                  enddo
*     F2 flavour by flavour
                  do i=1,nf
                     do pt=0,ipt
                        F2t = bq(i) 
     1                      * ( SC2(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * fph(jgrid,0,alpha+beta)
     3                      + SC2(jgrid,nf,2,pt,0,beta)   ! Singlet
     4                      * singlet
     5                      + SC2(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     6                      * ( fph(jgrid,i,alpha+beta) 
     7                      + fph(jgrid,-i,alpha+beta) ) )
*
                        F2(i,jgrid,alpha) = F2(i,jgrid,alpha)
     1                                    + as**pt * F2t
                     enddo
                  enddo
               enddo
*     F2 total
               do i=1,nf
                  F2(7,jgrid,alpha) = F2(7,jgrid,alpha)
     1                              + F2(i,jgrid,alpha)
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
     1                       + fph(jgrid,i,alpha+beta)
     2                       + fph(jgrid,-i,alpha+beta)
                  enddo
*     FL flavour by flavour
                  do i=1,nf
                     do pt=0,ipt
                        FLt = bq(i) 
     1                      * ( SCL(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * fph(jgrid,0,alpha+beta)
     3                      + SCL(jgrid,nf,2,pt,0,beta)   ! Singlet
     4                      * singlet
     5                      + SCL(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     6                      * ( fph(jgrid,i,alpha+beta) 
     7                      + fph(jgrid,-i,alpha+beta) ) )
*
                        FL(i,jgrid,alpha) = FL(i,jgrid,alpha)
     1                                    + as**pt * FLt
                     enddo
                  enddo
               enddo
*     FL total
               do i=1,nf
                  FL(7,jgrid,alpha) = FL(7,jgrid,alpha)
     1                              + FL(i,jgrid,alpha)
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
     1                      * SC3(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     2                      * ( fph(jgrid,i,alpha+beta) 
     3                      - fph(jgrid,-i,alpha+beta) )
*
                        F3(i,jgrid,alpha) = F3(i,jgrid,alpha)
     1                                    + as**pt * F3t
                     enddo
                  enddo
               enddo
*     F3 total
               do i=1,nf
                  F3(7,jgrid,alpha) = F3(7,jgrid,alpha)
     1                              + F3(i,jgrid,alpha)
               enddo
            enddo
         enddo
      elseif(ProcessDIS.eq."CC")then
         do jgrid=1,ngrid
            do alpha=0,nin(jgrid)
*
*     Rearrange PDFs in case the target is not a proton.
*     (This assumes isospin symmetry)
*
               if(TargetDIS(1:7).eq."proton")then
                  frac = 1d0
               elseif(TargetDIS(1:7).eq."neutron")then
                  frac = 0d0
               elseif(TargetDIS.eq."isoscalar")then
                  frac = 0.5d0
               elseif(TargetDIS(1:4).eq."iron")then
                  frac = 0.47166350921037d0 !23.403d0 / 49.618d0
               endif
*     Reweight up and down by frac
               do beta=0,nin(jgrid)
                  fdw = fph(jgrid,1,beta)
                  fdb = fph(jgrid,-1,beta)
                  fup = fph(jgrid,2,beta)
                  fub = fph(jgrid,-2,beta)
*
                  fph(jgrid,1,beta)  = frac * fdw + ( 1d0 - frac ) * fup
                  fph(jgrid,-1,beta) = frac * fdb + ( 1d0 - frac ) * fub
                  fph(jgrid,2,beta)  = frac * fup + ( 1d0 - frac ) * fdw
                  fph(jgrid,-2,beta) = frac * fub + ( 1d0 - frac ) * fdb
               enddo
*
               if(ProjectileDIS(1:8).eq."neutrino".or.
     1            ProjectileDIS(1:8).eq."positron")then
                  ipr = 1
               elseif(ProjectileDIS.eq."antineutrino".or.
     1            ProjectileDIS(1:8).eq."electron")then
                  ipr = - 1
               endif
*
*     F2
*
               do i=1,7
                  F2(i,jgrid,alpha) = 0d0
               enddo
               do beta=0,nin(jgrid)-alpha
*     F2light (all in the component 3 of F2)
                  do pt=0,ipt
                     F2t = 2d0
     1                   * ( SC2(jgrid,nf,1,pt,0,beta) ! Gluon
     2                   * ( V_ud2 + V_us2 )
     3                   * fph(jgrid,0,alpha+beta)
c     4                   + SC2(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                   * singlet
     6                   + SC2(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
     8                   + ( V_ud2 + V_us2 )
     9                   * fph(jgrid,-2*ipr,alpha+beta) 
     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
*
                     F2(3,jgrid,alpha) = F2(3,jgrid,alpha)
     1                                 + as**pt * F2t
*     F2charm
                     if(nf.ge.4)then
                        F2t = 2d0
     1                      * ( SC2(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * ( V_cd2 + V_cs2 )
     3                      * fph(jgrid,0,alpha+beta)
c     4                      + SC2(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                      * singlet
     6                      + SC2(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_cd2
     8                      * ( fph(jgrid,1*ipr,alpha+beta) 
     9                      + fph(jgrid,-4*ipr,alpha+beta) )
     1                      + V_cs2
     2                      * ( fph(jgrid,3*ipr,alpha+beta) 
     3                      + fph(jgrid,-4*ipr,alpha+beta) ) ) )
*             
                        F2(4,jgrid,alpha) = F2(4,jgrid,alpha)
     1                                    + as**pt * F2t
                     endif
*     F2bottom
                     if(nf.ge.5)then
                        F2t = 2d0
     1                      * ( SC2(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * ( V_ub2 + V_cb2 )
     3                      * fph(jgrid,0,alpha+beta)
c     4                      + SC2(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                      * singlet
     6                      + SC2(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_ub2
     8                      * ( fph(jgrid,-2*ipr,alpha+beta) 
     9                      + fph(jgrid,5*ipr,alpha+beta) )
     1                      + V_cb2
     2                      * ( fph(jgrid,-4*ipr,alpha+beta) 
     3                      + fph(jgrid,5*ipr,alpha+beta) ) ) )
*                     
                        F2(5,jgrid,alpha) = F2(5,jgrid,alpha)
     1                                    + as**pt * F2t
                     endif
*     F2top
                     if(nf.ge.6)then
                        F2t = 2d0 
     1                      * ( SC2(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * ( V_td2 + V_ts2 + V_tb2 )
     3                      * fph(jgrid,0,alpha+beta)
c     4                      + SC2(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                      * singlet
     6                      + SC2(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_td2
     8                      * ( fph(jgrid,1*ipr,alpha+beta) 
     9                      + fph(jgrid,-6*ipr,alpha+beta) )
     1                      + V_ts2
     2                      * ( fph(jgrid,3*ipr,alpha+beta) 
     3                      + fph(jgrid,-6*ipr,alpha+beta) )
     4                      + V_tb2
     5                      * ( fph(jgrid,5*ipr,alpha+beta) 
     6                      + fph(jgrid,-6*ipr,alpha+beta) ) ) )
*                  
                        F2(6,jgrid,alpha) = F2(6,jgrid,alpha)
     1                                    + as**pt * F2t
                     endif
                  enddo
               enddo
*     F2 total
               do i=3,nf
                  F2(7,jgrid,alpha) = F2(7,jgrid,alpha)
     1                              + F2(i,jgrid,alpha)
               enddo
*
*     FL
*
               do i=1,7
                  FL(i,jgrid,alpha) = 0d0
               enddo
               do beta=0,nin(jgrid)-alpha
*     FLlight (all in the component 3 of FL)
                  do pt=0,ipt
                     FLt = 2d0
     1                   * ( SCL(jgrid,nf,1,pt,0,beta) ! Gluon
     2                   * ( V_ud2 + V_us2 )
     3                   * fph(jgrid,0,alpha+beta)
c     4                   + SCL(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                   * singlet
     6                   + SCL(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
     8                   + ( V_ud2 + V_us2 )
     9                   * fph(jgrid,-2*ipr,alpha+beta) 
     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
*
                     FL(3,jgrid,alpha) = FL(3,jgrid,alpha)
     1                                 + as**pt * FLt
*     FLcharm
                     if(nf.ge.4)then
                        FLt = 2d0
     1                      * ( SCL(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * ( V_cd2 + V_cs2 )
     3                      * fph(jgrid,0,alpha+beta)
c     4                      + SCL(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                      * singlet
     6                      + SCL(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_cd2
     8                      * ( fph(jgrid,1*ipr,alpha+beta) 
     9                      + fph(jgrid,-4*ipr,alpha+beta) )
     1                      + V_cs2
     2                      * ( fph(jgrid,3*ipr,alpha+beta) 
     3                      + fph(jgrid,-4*ipr,alpha+beta) ) ) )
*                    
                        FL(4,jgrid,alpha) = FL(4,jgrid,alpha)
     1                                    + as**pt * FLt
                     endif
*     FLbottom
                     if(nf.ge.5)then
                        FLt = 2d0
     1                      * ( SCL(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * ( V_ub2 + V_cb2 )
     3                      * fph(jgrid,0,alpha+beta)
c     4                      + SCL(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                      * singlet
     6                      + SCL(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_ub2
     8                      * ( fph(jgrid,-2*ipr,alpha+beta) 
     9                      + fph(jgrid,5*ipr,alpha+beta) )
     1                      + V_cb2
     2                      * ( fph(jgrid,-4*ipr,alpha+beta) 
     3                      + fph(jgrid,5*ipr,alpha+beta) ) ) )
*                      
                        FL(5,jgrid,alpha) = FL(5,jgrid,alpha)
     1                                    + as**pt * FLt
                     endif
*     FLtop
                     if(nf.ge.6)then
                        FLt = 2d0 
     1                      * ( SCL(jgrid,nf,1,pt,0,beta) ! Gluon
     2                      * ( V_td2 + V_ts2 + V_tb2 )
     3                      * fph(jgrid,0,alpha+beta)
c     4                      + SCL(jgrid,nf,2,pt,0,beta)   ! Singlet
c     5                      * singlet
     6                      + SCL(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_td2
     8                      * ( fph(jgrid,1*ipr,alpha+beta) 
     9                      + fph(jgrid,-6*ipr,alpha+beta) )
     1                      + V_ts2
     2                      * ( fph(jgrid,3*ipr,alpha+beta) 
     3                      + fph(jgrid,-6*ipr,alpha+beta) )
     4                      + V_tb2
     5                      * ( fph(jgrid,5*ipr,alpha+beta) 
     6                      + fph(jgrid,-6*ipr,alpha+beta) ) ) )
*                    
                        FL(6,jgrid,alpha) = FL(6,jgrid,alpha)
     1                                    + as**pt * FLt
                     endif
                  enddo
               enddo
*     FL total
               do i=3,nf
                  FL(7,jgrid,alpha) = FL(7,jgrid,alpha)
     1                              + FL(i,jgrid,alpha)
               enddo
*
*     F3
*
               do i=1,7
                  F3(i,jgrid,alpha) = 0d0
               enddo
               do beta=0,nin(jgrid)-alpha
*     F3light (all in the component 3 of F3)
                  do pt=0,ipt
                     F3t = 2d0
     1                   * ( SC3(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
     8                   - ( V_ud2 + V_us2 )
     9                   * fph(jgrid,-2*ipr,alpha+beta) 
     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
*
                     F3(3,jgrid,alpha) = F3(3,jgrid,alpha)
     1                                 + as**pt * F3t
*     F3charm
                     if(nf.ge.4)then
                        F3t = 2d0
     1                      * ( SC3(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_cd2
     8                      * ( fph(jgrid,1*ipr,alpha+beta) 
     9                      - fph(jgrid,-4*ipr,alpha+beta) )
     1                      + V_cs2
     2                      * ( fph(jgrid,3*ipr,alpha+beta) 
     3                      - fph(jgrid,-4*ipr,alpha+beta) ) ) )
*                 
                        F3(4,jgrid,alpha) = F3(4,jgrid,alpha)
     1                                    + as**pt * F3t
                     endif
*     F3bottom
                     if(nf.ge.5)then
                        F3t = 2d0
     1                      * ( SC3(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_ub2
     8                      * ( - fph(jgrid,-2*ipr,alpha+beta) 
     9                      + fph(jgrid,5*ipr,alpha+beta) )
     1                      + V_cb2
     2                      * ( - fph(jgrid,-4*ipr,alpha+beta) 
     3                      + fph(jgrid,5*ipr,alpha+beta) ) ) )
*                  
                        F3(5,jgrid,alpha) = F3(5,jgrid,alpha)
     1                                    + as**pt * F3t
                     endif
*     F3top
                     if(nf.ge.6)then
                        F3t = 2d0 
     1                      * ( SC3(jgrid,nf,3,pt,0,beta)   ! Non-singlet
     7                      * ( V_td2
     8                      * ( fph(jgrid,1*ipr,alpha+beta) 
     9                      - fph(jgrid,-6*ipr,alpha+beta) )
     1                      + V_ts2
     2                      * ( fph(jgrid,3*ipr,alpha+beta) 
     3                      - fph(jgrid,-6*ipr,alpha+beta) )
     4                      + V_tb2
     5                      * ( fph(jgrid,5*ipr,alpha+beta) 
     6                      - fph(jgrid,-6*ipr,alpha+beta) ) ) )
*                    
                        F3(6,jgrid,alpha) = F3(6,jgrid,alpha)
     1                                    + as**pt * F3t
                     endif
                  enddo
               enddo
*     F3 total
               do i=3,nf
                  F3(7,jgrid,alpha) = F3(7,jgrid,alpha)
     1                              + F3(i,jgrid,alpha)
               enddo
            enddo
         enddo
      endif
*
      return
      end
