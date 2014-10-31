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
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/m2th.h"
      include "../commons/StructureFunctions.h"
      include "../commons/TargetDIS.h"
      include "../commons/ProcessDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/CKM.h"
      include "../commons/MassScheme.h"
      include "../commons/Nf_FF.h"
**
*     Internal Variables
*
      double precision Q
**
*     Internal Variables
*
      integer jgrid
      integer nf
      integer i,j,k
      integer alpha,beta
      integer pt
      integer ipr
      integer ixi(4:6)
      double precision t1,t2
      double precision Q2,W2,xi(4:6),c0(4:6),c1(4:6)
      double precision singlet
      double precision F2t,FLt,F3t
      double precision as,a_QCD
      double precision bq(6),dq(6)
      double precision fup,fub,fdw,fdb
      double precision frac
      double precision diff(nxi),sgn
      double precision SC2(0:ngrid_max,6,3,0:2,0:nint_max,0:nint_max)
      double precision SCL(0:ngrid_max,6,3,0:2,0:nint_max,0:nint_max)
      double precision SC3(0:ngrid_max,6,3,0:2,0:nint_max,0:nint_max)
      double precision damp(4:6)
*
      call cpu_time(t1)
*
      Q2    = Q * Q
      Q2DIS = Q2
*
*     Dumping factor for FONLL
*
      do i=4,6
         if(q2.gt.m2th(i))then
            damp(i) = ( 1d0 - m2th(i) / Q2 )**2d0
         else
            damp(i) = 0d0
         endif
      enddo
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
*     Find "ixi" such that "xigrid(ixi)" < "xi" < "xigrid(ixi+1)"
*
      do j=4,6
         ixi(j) = 0
         xi(j)  = Q2 / m2th(j)
         if(xi(j).le.ximin)then
            ixi(j) = 1
         elseif(xi(j).ge.ximax)then
            ixi(j) = nxi
         else
            diff(1) = xi(j) - xigrid(1)
            do i=2,nxi
               diff(i) = xi(j) - xigrid(i)
               sgn = diff(i-1) * diff(i)
               if(sgn.lt.0.d0)then
                  ixi(j) = i - 1
               endif
            enddo
         endif
*
*     Coefficients of the linear interpolation on the xi grid
*
         c0(j) = dlog(xigrid(ixi(j)+1)/xi(j))
     1         / dlog(xigrid(ixi(j)+1)/xigrid(ixi(j)))
         c1(j) = dlog(xi(j)/xigrid(ixi(j)))
     1         / dlog(xigrid(ixi(j)+1)/xigrid(ixi(j)))
      enddo
*
*     Compute alphas
*
      as = a_QCD(Q2)
*
*     Electromagnetic and Neutral current structure functions
*
      if(ProcessDIS.eq."EM".or.ProcessDIS.eq."NC")then
*
*     Construct the coefficient functions according to the scheme chosen
*
         if(MassScheme.eq."ZM-VFNS")then
            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  do pt=0,ipt
                     do k=1,3
                        do i=1,nf
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,nf,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,nf,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) = 
     1                          SC3zm(jgrid,nf,k,pt,0,alpha)
                        enddo
                        if(nf.lt.6)then
                           do i=nf+1,6
                              SC2(jgrid,i,k,pt,0,alpha) = 0d0
                              SCL(jgrid,i,k,pt,0,alpha) = 0d0
                              SC3(jgrid,i,k,pt,0,alpha) = 0d0
                           enddo
                        endif
                     enddo
                  enddo
               enddo
            enddo
         elseif(MassScheme(1:4).eq."FFNS")then
            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
*
*     Light coefficient functions
*
                  do i=1,Nf_FF
                     do pt=0,ipt
                        do k=1,3
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) = 
     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
                        enddo
*
*     Add the massive part at NNLO
*
                        if(pt.ge.2)then
                           if(Nf_FF.lt.6)then
                              do j=Nf_FF+1,6
                                 if(W2.ge.4d0*m2th(j))then
                                    SC2(jgrid,i,3,pt,0,alpha) = 
     1                              SC2(jgrid,i,3,pt,0,alpha)
     2                     + c0(j) * SC2m(1,jgrid,ixi(j),3,pt,0,alpha)
     3                     + c1(j) * SC2m(1,jgrid,ixi(j)+1,3,pt,0,alpha)
                                    SCL(jgrid,i,3,pt,0,alpha) = 
     1                              SCL(jgrid,i,3,pt,0,alpha)
     2                     + c0(j) * SCLm(1,jgrid,ixi(j),3,pt,0,alpha)
     3                     + c1(j) * SCLm(1,jgrid,ixi(j)+1,3,pt,0,alpha)
                                 endif
                              enddo
                           endif
                        endif
                     enddo
                  enddo
*
*     Heavy coefficient functions
*
                  if(Nf_FF.lt.6)then
                     do i=Nf_FF+1,6
                        do pt=0,ipt
                           do k=1,2
                              if(W2.ge.4d0*m2th(i))then
                                 SC2(jgrid,i,k,pt,0,alpha) = 
     1                       c0(i) * SC2m(1,jgrid,ixi(i),k,pt,0,alpha)
     2                     + c1(i) * SC2m(1,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SCL(jgrid,i,k,pt,0,alpha) =
     1                       c0(i) * SCLm(1,jgrid,ixi(i),k,pt,0,alpha)
     2                     + c1(i) * SCLm(1,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SC3(jgrid,i,k,pt,0,alpha) = 
     1                       c0(i) * SC3m(1,jgrid,ixi(i),k,pt,0,alpha)
     2                     + c1(i) * SC3m(1,jgrid,ixi(i)+1,k,pt,0,alpha)
                             else
                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
                              endif
                           enddo
                           SC2(jgrid,i,3,pt,0,alpha) = 0d0
                           SCL(jgrid,i,3,pt,0,alpha) = 0d0
                           SC3(jgrid,i,3,pt,0,alpha) = 0d0
                        enddo
                     enddo
                  endif
               enddo
            enddo
         elseif(MassScheme(1:4).eq."FFN0")then
            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
*
*     Light coefficient functions
*
                  do i=1,Nf_FF
                     do pt=0,ipt
                        do k=1,3
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) =
     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
                        enddo
*
*     Add the massive part at NNLO
*
                        if(pt.ge.2)then
                           if(Nf_FF.lt.6)then
                              do j=Nf_FF+1,6
                                 if(W2.ge.4d0*m2th(j))then
                                    SC2(jgrid,i,3,pt,0,alpha) = 
     1                              SC2(jgrid,i,3,pt,0,alpha)
     2                    + c0(j) * SC2m0(1,jgrid,ixi(j),3,pt,0,alpha)
     3                    + c1(j) * SC2m0(1,jgrid,ixi(j)+1,3,pt,0,alpha)
                                    SCL(jgrid,i,3,pt,0,alpha) = 
     1                              SCL(jgrid,i,3,pt,0,alpha)
     2                    + c0(j) * SCLm0(1,jgrid,ixi(j),3,pt,0,alpha)
     3                    + c1(j) * SCLm0(1,jgrid,ixi(j)+1,3,pt,0,alpha)
                                 endif
                              enddo
                           endif
                        endif
                     enddo
                  enddo
*
*     Heavy coefficient functions
*
                  if(Nf_FF.lt.6)then
                     do i=Nf_FF+1,6
                        do pt=0,ipt
                           do k=1,2
                              if(W2.ge.4d0*m2th(i))then
                                 SC2(jgrid,i,k,pt,0,alpha) = 
     1                      c0(i) * SC2m0(1,jgrid,ixi(i),k,pt,0,alpha)
     2                    + c1(i) * SC2m0(1,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SCL(jgrid,i,k,pt,0,alpha) = 
     1                      c0(i) * SCLm0(1,jgrid,ixi(i),k,pt,0,alpha)
     2                    + c1(i) * SCLm0(1,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SC3(jgrid,i,k,pt,0,alpha) = 
     1                      c0(i) * SC3m0(1,jgrid,ixi(i),k,pt,0,alpha)
     2                    + c1(i) * SC3m0(1,jgrid,ixi(i)+1,k,pt,0,alpha)
                              else
                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
                              endif
                           enddo
                           SC2(jgrid,i,3,pt,0,alpha) = 0d0
                           SCL(jgrid,i,3,pt,0,alpha) = 0d0
                           SC3(jgrid,i,3,pt,0,alpha) = 0d0
                        enddo
                     enddo
                  endif
               enddo
            enddo
         elseif(MassScheme(1:5).eq."FONLL")then








            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
*
*     Light coefficient functions
*
                  do i=1,Nf_FF
                     do pt=0,ipt
                        do k=1,3
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,nf,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,nf,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) = 
     1                          SC3zm(jgrid,nf,k,pt,0,alpha)
                        enddo
c$$$*
c$$$*     Add the massive part at NNLO
c$$$*
c$$$                        if(pt.ge.2)then
c$$$                           if(Nf_FF.lt.6)then
c$$$                              do j=Nf_FF+1,6
c$$$                                 if(W2.ge.4d0*m2th(j))then
c$$$                                    SC2(jgrid,i,3,pt,0,alpha) = 
c$$$     1                              SC2(jgrid,i,3,pt,0,alpha)
c$$$     2                     + c0(j) * SC2m(1,jgrid,ixi(j),3,pt,0,alpha)
c$$$     3                     + c1(j) * SC2m(1,jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$                                    SCL(jgrid,i,3,pt,0,alpha) = 
c$$$     1                              SCL(jgrid,i,3,pt,0,alpha)
c$$$     2                     + c0(j) * SCLm(1,jgrid,ixi(j),3,pt,0,alpha)
c$$$     3                     + c1(j) * SCLm(1,jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$                                 endif
c$$$                              enddo
c$$$                           endif
c$$$                        endif
                     enddo
                  enddo
*
*     Heavy coefficient functions
*
                  if(Nf_FF.lt.6)then
                     do i=Nf_FF+1,6
                        do pt=0,ipt
                           do k=1,3
                              SC2(jgrid,i,k,pt,0,alpha) = 0d0
                              SCL(jgrid,i,k,pt,0,alpha) = 0d0
                              SC3(jgrid,i,k,pt,0,alpha) = 0d0
*
*     Zero Mass Part
*
                              if(Q2.ge.m2th(i))then
                                 SC2(jgrid,i,k,pt,0,alpha) = 
     1                           SC2(jgrid,i,k,pt,0,alpha)
     2                         + damp(i) * SC2zm(jgrid,nf,k,pt,0,alpha)
*
                                 SCL(jgrid,i,k,pt,0,alpha) = 
     1                           SCL(jgrid,i,k,pt,0,alpha)
     2                         + damp(i) * SCLzm(jgrid,nf,k,pt,0,alpha)
*
                                 SC3(jgrid,i,k,pt,0,alpha) = 
     1                           SC3(jgrid,i,k,pt,0,alpha)
     2                         + damp(i) * SC3zm(jgrid,nf,k,pt,0,alpha)
                              endif
                           enddo
*
*     Massive Parts
*
                           do k=1,2
                              if(W2.ge.4d0*m2th(i))then
                                 SC2(jgrid,i,k,pt,0,alpha) =
     1                           SC2(jgrid,i,k,pt,0,alpha)
     2               + c0(i) * ( SC2m(1,jgrid,ixi(i),k,pt,0,alpha)
     3               - damp(i) * SC2m0(1,jgrid,ixi(i),k,pt,0,alpha) )
     4               + c1(i) * ( SC2m(1,jgrid,ixi(i)+1,k,pt,0,alpha)
     5               - damp(i) * SC2m0(1,jgrid,ixi(i)+1,k,pt,0,alpha) )

                                 SCL(jgrid,i,k,pt,0,alpha) =
     1                           SCL(jgrid,i,k,pt,0,alpha)
     2               + c0(i) * ( SCLm(1,jgrid,ixi(i),k,pt,0,alpha)
     3               - damp(i) * SCLm0(1,jgrid,ixi(i),k,pt,0,alpha) )
     4               + c1(i) * ( SCLm(1,jgrid,ixi(i)+1,k,pt,0,alpha)
     5               - damp(i) * SCLm0(1,jgrid,ixi(i)+1,k,pt,0,alpha) )

                                 SC3(jgrid,i,k,pt,0,alpha) =
     1                           SC3(jgrid,i,k,pt,0,alpha)
     2               + c0(i) * ( SC3m(1,jgrid,ixi(i),k,pt,0,alpha)
     3               - damp(i) * SC3m0(1,jgrid,ixi(i),k,pt,0,alpha) )
     4               + c1(i) * ( SC3m(1,jgrid,ixi(i)+1,k,pt,0,alpha)
     5               - damp(i) * SC3m0(1,jgrid,ixi(i)+1,k,pt,0,alpha) )
                              endif
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo







         endif
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
                  do i=1,6
                     do pt=0,ipt
                        F2t = bq(i) 
     1                      * ( SC2(jgrid,i,1,pt,0,beta) ! Gluon
     2                      * fph(jgrid,0,alpha+beta)
     3                      + SC2(jgrid,i,2,pt,0,beta)   ! Singlet
     4                      * singlet
     5                      + SC2(jgrid,i,3,pt,0,beta)   ! Non-singlet
     6                      * ( fph(jgrid,i,alpha+beta) 
     7                      + fph(jgrid,-i,alpha+beta) ) )
*
                        F2(i,jgrid,alpha) = F2(i,jgrid,alpha)
     1                                    + as**pt * F2t
                     enddo
                  enddo
               enddo
*     F2 total (combine according to the massive scheme)
               do i=1,6
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
                  do i=1,6
                     do pt=0,ipt
                        FLt = bq(i) 
     1                      * ( SCL(jgrid,i,1,pt,0,beta) ! Gluon
     2                      * fph(jgrid,0,alpha+beta)
     3                      + SCL(jgrid,i,2,pt,0,beta)   ! Singlet
     4                      * singlet
     5                      + SCL(jgrid,i,3,pt,0,beta)   ! Non-singlet
     6                      * ( fph(jgrid,i,alpha+beta) 
     7                      + fph(jgrid,-i,alpha+beta) ) )
*
                        FL(i,jgrid,alpha) = FL(i,jgrid,alpha)
     1                                    + as**pt * FLt
                     enddo
                  enddo
               enddo
*     FL total (combine according to the massive scheme)
               do i=1,6
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
     1                      * SC3(jgrid,i,3,pt,0,beta)   ! Non-singlet
     2                      * ( fph(jgrid,i,alpha+beta) 
     3                      - fph(jgrid,-i,alpha+beta) )
*
                        F3(i,jgrid,alpha) = F3(i,jgrid,alpha)
     1                                    + as**pt * F3t
                     enddo
                  enddo
               enddo
*     F3 total (it is always in the ZM scheme)
               do i=1,6
                  F3(7,jgrid,alpha) = F3(7,jgrid,alpha)
     1                              + F3(i,jgrid,alpha)
               enddo
            enddo
         enddo
      elseif(ProcessDIS.eq."CC")then
*
*     Construct the coefficient functions according to the scheme chosen
*
         if(MassScheme.eq."ZM-VFNS")then
            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  do pt=0,ipt
                     do k=1,3
                        do i=3,nf
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,nf,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,nf,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) = 
     1                          SC3zm(jgrid,nf,k,pt,0,alpha)
                        enddo
                        if(nf.lt.6)then
                           do i=nf+1,6
                              SC2(jgrid,i,k,pt,0,alpha) = 0d0
                              SCL(jgrid,i,k,pt,0,alpha) = 0d0
                              SC3(jgrid,i,k,pt,0,alpha) = 0d0
                           enddo
                        endif
                     enddo
                  enddo
               enddo
            enddo
         elseif(MassScheme(1:4).eq."FFNS")then
            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
*
*     Light coefficient functions
*
                  do i=3,Nf_FF
                     do pt=0,ipt
                        do k=1,3
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) = 
     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
                        enddo

                     enddo
                  enddo
*
*     Heavy coefficient functions
*
                  if(Nf_FF.lt.6)then
                     do i=Nf_FF+1,6
                        do pt=0,ipt
                           do k=1,3
                              if(W2.ge.m2th(i))then
                                 SC2(jgrid,i,k,pt,0,alpha) = 
     1                       c0(i) * SC2m(2,jgrid,ixi(i),k,pt,0,alpha)
     2                     + c1(i) * SC2m(2,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SCL(jgrid,i,k,pt,0,alpha) =
     1                       c0(i) * SCLm(2,jgrid,ixi(i),k,pt,0,alpha)
     2                     + c1(i) * SCLm(2,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SC3(jgrid,i,k,pt,0,alpha) = 
     1                       c0(i) * SC3m(2,jgrid,ixi(i),k,pt,0,alpha)
     2                     + c1(i) * SC3m(2,jgrid,ixi(i)+1,k,pt,0,alpha)
                              else
                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
                              endif
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         elseif(MassScheme(1:4).eq."FFN0")then
            do jgrid=1,ngrid
               do alpha=0,nin(jgrid)
                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
*
*     Light coefficient functions
*
                  do i=1,Nf_FF
                     do pt=0,ipt
                        do k=1,3
                           SC2(jgrid,i,k,pt,0,alpha) = 
     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
                           SCL(jgrid,i,k,pt,0,alpha) = 
     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
                           SC3(jgrid,i,k,pt,0,alpha) =
     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
                        enddo
                     enddo
                  enddo
*
*     Heavy coefficient functions
*
                  if(Nf_FF.lt.6)then
                     do i=Nf_FF+1,6
                        do pt=0,ipt
                           do k=1,3
                              if(W2.ge.m2th(i))then
                                 SC2(jgrid,i,k,pt,0,alpha) = 
     1                      c0(i) * SC2m0(2,jgrid,ixi(i),k,pt,0,alpha)
     2                    + c1(i) * SC2m0(2,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SCL(jgrid,i,k,pt,0,alpha) =
     1                      c0(i) * SCLm0(2,jgrid,ixi(i),k,pt,0,alpha)
     2                    + c1(i) * SCLm0(2,jgrid,ixi(i)+1,k,pt,0,alpha)
                                 SC3(jgrid,i,k,pt,0,alpha) = 
     1                      c0(i) * SC3m0(2,jgrid,ixi(i),k,pt,0,alpha)
     2                    + c1(i) * SC3m0(2,jgrid,ixi(i)+1,k,pt,0,alpha)
                              else
                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
                              endif
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         endif





*
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
                  singlet = 0d0
                  do i=1,nf
                     singlet = singlet 
     1                       + fph(jgrid,i,alpha+beta)
     2                       + fph(jgrid,-i,alpha+beta)
                  enddo
*     F2light (all in the component 3 of F2)
                  do pt=0,ipt
                     F2t = 2d0
     1                   * ( ( V_ud2 + V_us2 )        ! Gluon
     2                   * ( SC2(jgrid,3,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SC2(jgrid,3,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SC2(jgrid,3,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
     8                   + ( V_ud2 + V_us2 )
     9                   * fph(jgrid,-2*ipr,alpha+beta) 
     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
*
                     F2(3,jgrid,alpha) = F2(3,jgrid,alpha)
     1                                 + as**pt * F2t
*     F2charm
                     F2t = 2d0
     1                   * ( ( V_cd2 + V_cs2 )         ! Gluon
     2                   * ( SC2(jgrid,4,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SC2(jgrid,4,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SC2(jgrid,4,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_cd2
     8                   * ( fph(jgrid,1*ipr,alpha+beta)
     9                   + fph(jgrid,-4*ipr,alpha+beta) )
     1                   + V_cs2
     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
     3                   + fph(jgrid,-4*ipr,alpha+beta) ) ) )
*             
                     F2(4,jgrid,alpha) = F2(4,jgrid,alpha)
     1                                 + as**pt * F2t
*     F2bottom
                     F2t = 2d0
     1                   * ( ( V_ub2 + V_cb2 )         ! Gluon
     2                   * ( SC2(jgrid,5,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SC2(jgrid,5,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SC2(jgrid,5,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ub2
     8                   * ( fph(jgrid,-2*ipr,alpha+beta) 
     9                   + fph(jgrid,5*ipr,alpha+beta) )
     1                   + V_cb2
     2                   * ( fph(jgrid,-4*ipr,alpha+beta) 
     3                   + fph(jgrid,5*ipr,alpha+beta) ) ) )
*                    
                     F2(5,jgrid,alpha) = F2(5,jgrid,alpha)
     1                                 + as**pt * F2t
*     F2top
                     F2t = 2d0 
     1                   * ( ( V_td2 + V_ts2 + V_tb2 ) ! Gluon
     2                   * ( SC2(jgrid,6,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SC2(jgrid,6,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SC2(jgrid,6,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_td2
     8                   * ( fph(jgrid,1*ipr,alpha+beta) 
     9                   + fph(jgrid,-6*ipr,alpha+beta) )
     1                   + V_ts2
     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
     3                   + fph(jgrid,-6*ipr,alpha+beta) )
     4                   + V_tb2
     5                   * ( fph(jgrid,5*ipr,alpha+beta) 
     6                   + fph(jgrid,-6*ipr,alpha+beta) ) ) )
*                  
                     F2(6,jgrid,alpha) = F2(6,jgrid,alpha)
     1                                 + as**pt * F2t
                  enddo
               enddo
*     F2 total
               do i=3,6
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
*     FLlight (all in the component 3 of FL)
                  do pt=0,ipt
                     FLt = 2d0
     1                   * ( ( V_ud2 + V_us2 )         ! Gluon
     2                   * ( SCL(jgrid,3,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SCL(jgrid,3,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SCL(jgrid,3,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
     8                   + ( V_ud2 + V_us2 )
     9                   * fph(jgrid,-2*ipr,alpha+beta) 
     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
*
                     FL(3,jgrid,alpha) = FL(3,jgrid,alpha)
     1                                 + as**pt * FLt
*     FLcharm
                     FLt = 2d0
     1                   * ( ( V_cd2 + V_cs2 )         ! Gluon
     2                   * ( SCL(jgrid,4,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SCL(jgrid,4,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SCL(jgrid,4,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_cd2
     8                   * ( fph(jgrid,1*ipr,alpha+beta) 
     9                   + fph(jgrid,-4*ipr,alpha+beta) )
     1                   + V_cs2
     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
     3                   + fph(jgrid,-4*ipr,alpha+beta) ) ) )
*                 
                     FL(4,jgrid,alpha) = FL(4,jgrid,alpha)
     1                                 + as**pt * FLt
*     FLbottom
                     FLt = 2d0
     1                   * ( ( V_ub2 + V_cb2 )         ! Gluon
     2                   * ( SCL(jgrid,5,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SCL(jgrid,5,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SCL(jgrid,5,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_ub2
     8                   * ( fph(jgrid,-2*ipr,alpha+beta) 
     9                   + fph(jgrid,5*ipr,alpha+beta) )
     1                   + V_cb2
     2                   * ( fph(jgrid,-4*ipr,alpha+beta) 
     3                   + fph(jgrid,5*ipr,alpha+beta) ) ) )
*                   
                     FL(5,jgrid,alpha) = FL(5,jgrid,alpha)
     1                                 + as**pt * FLt
*     FLtop
                     FLt = 2d0 
     1                   * ( ( V_td2 + V_ts2 + V_tb2 ) ! Gluon
     2                   * ( SCL(jgrid,6,1,pt,0,beta)
     3                   * fph(jgrid,0,alpha+beta)
     4                   + SCL(jgrid,6,2,pt,0,beta)   ! Singlet
     5                   * singlet )
     6                   + SCL(jgrid,6,3,pt,0,beta)   ! Non-singlet
     7                   * ( V_td2
     8                   * ( fph(jgrid,1*ipr,alpha+beta) 
     9                   + fph(jgrid,-6*ipr,alpha+beta) )
     1                   + V_ts2
     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
     3                   + fph(jgrid,-6*ipr,alpha+beta) )
     4                   + V_tb2
     5                   * ( fph(jgrid,5*ipr,alpha+beta) 
     6                   + fph(jgrid,-6*ipr,alpha+beta) ) ) )
*                 
                     FL(6,jgrid,alpha) = FL(6,jgrid,alpha)
     1                                 + as**pt * FLt
                  enddo
               enddo
*     FL total
               do i=3,6
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
     1                   * ( SC3(jgrid,3,3,pt,0,beta)      ! Non-singlet
     2                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
     3                   - ( V_ud2 + V_us2 )
     4                   * fph(jgrid,-2*ipr,alpha+beta) 
     5                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
*
                     F3(3,jgrid,alpha) = F3(3,jgrid,alpha)
     1                                 + as**pt * F3t
*     F3charm
                     F3t = 2d0 * ( ( V_cd2 + V_cs2 )         ! Gluon
     1                   * SC3(jgrid,4,1,pt,0,beta)
     2                   * fph(jgrid,0,alpha+beta)
     3                   + SC3(jgrid,4,3,pt,0,beta)          ! Non-singlet
     4                   * ( V_cd2
     5                   * ( fph(jgrid,1*ipr,alpha+beta) 
     6                   - fph(jgrid,-4*ipr,alpha+beta) )
     7                   + V_cs2
     8                   * ( fph(jgrid,3*ipr,alpha+beta) 
     9                   - fph(jgrid,-4*ipr,alpha+beta) ) ) )
*                 
                     F3(4,jgrid,alpha) = F3(4,jgrid,alpha)
     1                                 + as**pt * F3t
*     F3bottom
                     F3t = 2d0
     1                   * ( SC3(jgrid,5,3,pt,0,beta)   ! Non-singlet
     2                   * ( V_ub2
     3                   * ( - fph(jgrid,-2*ipr,alpha+beta) 
     4                   + fph(jgrid,5*ipr,alpha+beta) )
     5                   + V_cb2
     6                   * ( - fph(jgrid,-4*ipr,alpha+beta) 
     7                   + fph(jgrid,5*ipr,alpha+beta) ) ) )
*                 
                     F3(5,jgrid,alpha) = F3(5,jgrid,alpha)
     1                                 + as**pt * F3t
*     F3top
                     F3t = 2d0 
     1                   * ( SC3(jgrid,6,3,pt,0,beta)   ! Non-singlet
     2                   * ( V_td2
     3                   * ( fph(jgrid,1*ipr,alpha+beta) 
     4                   - fph(jgrid,-6*ipr,alpha+beta) )
     5                   + V_ts2
     6                   * ( fph(jgrid,3*ipr,alpha+beta) 
     7                   - fph(jgrid,-6*ipr,alpha+beta) )
     8                   + V_tb2
     9                   * ( fph(jgrid,5*ipr,alpha+beta) 
     1                   - fph(jgrid,-6*ipr,alpha+beta) ) ) )
*                 
                     F3(6,jgrid,alpha) = F3(6,jgrid,alpha)
     1                                 + as**pt * F3t
                  enddo
               enddo
*     F3 total
               do i=3,6
                  F3(7,jgrid,alpha) = F3(7,jgrid,alpha)
     1                              + F3(i,jgrid,alpha)
               enddo
            enddo
         enddo
      endif
*
      call cpu_time(t2)
*
      write(6,"(a,f7.3,a)") " Convolution completed in",t2-t1," s"
      write(6,*) " "
*
      return
      end
