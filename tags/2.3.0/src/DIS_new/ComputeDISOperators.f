************************************************************************
*
*     ComputeDISOperators.f:
*
*     This routine construct the DIS operators to be convoluted with the
*     final state PDFs to obtain the DIS structure functions.
*
************************************************************************
      subroutine ComputeDISOperators(Q)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/DISOperators.h"
      include "../commons/m2th.h"
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
      integer jgrid,ipdf,ihq,pt
      integer alpha,beta
      integer nf
      integer gbound
      integer ipr
      integer i,ixi(4:6)
      double precision Q2,W2
      double precision as(0:2),a_QCD
      double precision bq(6),dq(6)
      double precision frac,fr3
      double precision C2g(3:6),C2ps(3:6),C2nsp(3:6),C2nsm(3:6)
      double precision CLg(3:6),CLps(3:6),CLnsp(3:6),CLnsm(3:6)
      double precision C3g(3:6),C3nsp(3:6),C3nsm(3:6)
      double precision Kl,Kc,Kb,Kt
      double precision sgn,diff(nxi),xi(4:6),c0(4:6),c1(4:6)
      double precision t1,t2
*
      call cpu_time(t1)
*
*     Initialize operators
*
      do jgrid=1,ngrid
         gbound = 0
         if(IsExt(jgrid)) gbound = nin(jgrid)
         do alpha=0,gbound
            do beta=alpha,nin(jgrid)
               do ihq=3,7
                  do ipdf=0,13
                     OpF2(jgrid,ihq,ipdf,alpha,beta) = 0d0
                     OpFL(jgrid,ihq,ipdf,alpha,beta) = 0d0
                     OpF3(jgrid,ihq,ipdf,alpha,beta) = 0d0
                  enddo
               enddo
            enddo
         enddo
      enddo
*
*     Final scale
*
      Q2 = Q * Q
*
*     Compute alphas
*
      as(0) = 1d0
      as(1) = a_QCD(Q2)
      as(2) = as(1) * as(1)
*
*     Find number of active flavours at the scale Q2
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
c      if(MassScheme(1:3).eq."FFN") nf = Nf_FF
*
*     Define the ratio # of protons / # of nucleons of the target
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
*
*     Factor to use in front of T3 and V3
*
      fr3 = ( 2d0 * frac - 1d0 )
*
*     Find "ixi" such that "xigrid(ixi)" < "xi" < "xigrid(ixi+1)"
*
      do ihq=4,6
         ixi(ihq) = 0
         xi(ihq)  = Q2 / m2th(ihq)
         if(xi(ihq).le.ximin)then
            ixi(ihq) = 1
         elseif(xi(ihq).ge.ximax)then
            ixi(ihq) = nxi
         else
            diff(1) = xi(ihq) - xigrid(1)
            do i=2,nxi
               diff(i) = xi(ihq) - xigrid(i)
               sgn = diff(i-1) * diff(i)
               if(sgn.lt.0.d0)then
                  ixi(ihq) = i - 1
               endif
            enddo
         endif
*
*     Coefficients of the linear interpolation on the xi grid
*
         c0(ihq) = dlog(xigrid(ixi(ihq)+1)/xi(ihq))
     1           / dlog(xigrid(ixi(ihq)+1)/xigrid(ixi(ihq)))
         c1(ihq) = dlog(xi(ihq)/xigrid(ixi(ihq)))
     1           / dlog(xigrid(ixi(ihq)+1)/xigrid(ixi(ihq)))
      enddo
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
*     
*     If an external grid is found compute the whole operator
*     
            gbound = 0
            if(IsExt(jgrid)) gbound = nin(jgrid)
*     
            do alpha=0,gbound
               do beta=alpha,nin(jgrid)
*
*     Final state energy
*
                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
*     
*     Put together coefficient functions
*     
                  do ihq=3,6
                     C2g(ihq)   = 0d0
                     C2ps(ihq)  = 0d0
                     C2nsp(ihq) = 0d0
                     CLg(ihq)   = 0d0
                     CLps(ihq)  = 0d0
                     CLnsp(ihq) = 0d0
                     C3nsm(ihq) = 0d0
                  enddo
*     
*     Construct coefficient functions according to the scheme chosen
*     
                  if(MassScheme.eq."ZM-VFNS")then
*
*     Light coefficient functions
*
                     do pt=0,ipt
                        C2g(3)   = C2g(3)
     1                       + as(pt) * SC2zm(jgrid,nf,1,pt,alpha,beta)
                        C2ps(3)  = C2ps(3)
     1                       + as(pt) * SC2zm(jgrid,nf,2,pt,alpha,beta)
                        C2nsp(3) = C2nsp(3)
     1                       + as(pt) * SC2zm(jgrid,nf,3,pt,alpha,beta)
                        CLg(3)   = CLg(3)
     1                       + as(pt) * SCLzm(jgrid,nf,1,pt,alpha,beta)
                        CLps(3)  = CLps(3)
     1                       + as(pt) * SCLzm(jgrid,nf,2,pt,alpha,beta)
                        CLnsp(3) = CLnsp(3)
     1                       + as(pt) * SCLzm(jgrid,nf,3,pt,alpha,beta)
                        C3nsm(3) = C3nsm(3)
     1                       + as(pt) * SC3zm(jgrid,nf,4,pt,alpha,beta)
                     enddo
*
*     Heavy quark coefficient functions
*
                     if(nf.gt.3)then
                        do ihq=4,nf
                           C2g(ihq)   = C2g(3)
                           C2ps(ihq)  = C2ps(3)
                           C2nsp(ihq) = C2nsp(3)
                           CLg(ihq)   = CLg(3)
                           CLps(ihq)  = CLps(3)
                           CLnsp(ihq) = CLnsp(3)
                           C3nsm(ihq) = C3nsm(3)
                        enddo
                     endif
                  elseif(MassScheme.eq."FFN0")then
*
*     Light coefficient functions
*
                     do pt=0,ipt
                        C2g(3)   = C2g(3) + as(pt)
     1                       * SC2zm(jgrid,Nf_FF,1,pt,alpha,beta)
                        C2ps(3)  = C2ps(3) + as(pt)
     1                       * SC2zm(jgrid,Nf_FF,2,pt,alpha,beta)
                        C2nsp(3) = C2nsp(3) + as(pt)
     1                       * SC2zm(jgrid,Nf_FF,3,pt,alpha,beta)
                        CLg(3)   = CLg(3) + as(pt)
     1                       * SCLzm(jgrid,Nf_FF,1,pt,alpha,beta)
                        CLps(3)  = CLps(3) + as(pt)
     1                       * SCLzm(jgrid,Nf_FF,2,pt,alpha,beta)
                        CLnsp(3) = CLnsp(3) + as(pt)
     1                       * SCLzm(jgrid,Nf_FF,3,pt,alpha,beta)
                        C3nsm(3) = C3nsm(3) + as(pt)
     1                       * SC3zm(jgrid,Nf_FF,4,pt,alpha,beta)
                     enddo
                     if(Nf_FF.gt.3)then
                        do ihq=4,Nf_FF
                           C2g(ihq)   = C2g(3)
                           C2ps(ihq)  = C2ps(3)
                           C2nsp(ihq) = C2nsp(3)
                           CLg(ihq)   = CLg(3)
                           CLps(ihq)  = CLps(3)
                           CLnsp(ihq) = CLnsp(3)
                           C3nsm(ihq) = C3nsm(3)
                        enddo
                     endif
*
*     Add the massive part at NNLO to the light CFs
*
                     if(ipt.ge.2)then
                        if(Nf_FF.lt.6)then
                           do ihq=Nf_FF+1,6
                              if(W2.ge.4d0*m2th(ihq))then
                                 C2nsp(3) = C2nsp(3) + as(2) *
     1                                ( c0(ihq)
     2                                * SC2m0NC(jgrid,ixi(ihq),
     3                                3,2,alpha,beta)
     4                                + c1(ihq)
     5                                * SC2m0NC(jgrid,ixi(ihq)+1,
     6                                3,2,alpha,beta) )
                                 CLnsp(3) = CLnsp(3) + as(2) *
     1                                ( c0(ihq)
     2                                * SCLm0NC(jgrid,ixi(ihq),
     3                                3,2,alpha,beta)
     4                                + c1(ihq)
     5                                * SCLm0NC(jgrid,ixi(ihq)+1,
     6                                3,2,alpha,beta) )
                              endif
                           enddo
                        endif
                     endif
*
*     Heavy quark coefficient functions
*
                  if(Nf_FF.lt.6)then
                     do ihq=Nf_FF+1,6
                        if(W2.ge.4d0*m2th(ihq))then
                           do pt=0,ipt
                              C2g(ihq) = C2g(ihq) + as(pt)
     1                             * ( c0(ihq)
     2                             * SC2m0NC(jgrid,ixi(ihq),
     3                             1,pt,0,alpha)
     4                             + c1(ihq)
     5                             * SC2m0NC(jgrid,ixi(ihq)+1,
     6                             1,pt,0,alpha) )
                              C2ps(ihq) = C2ps(ihq) + as(pt)
     1                             * ( c0(ihq)
     2                             * SC2m0NC(jgrid,ixi(ihq),
     3                             2,pt,0,alpha)
     4                             + c1(ihq)
     5                             * SC2m0NC(jgrid,ixi(ihq)+1,
     6                             2,pt,0,alpha) )
                              CLg(ihq) = CLg(ihq) + as(pt)
     1                             * ( c0(ihq)
     2                             * SCLm0NC(jgrid,ixi(ihq),
     3                             1,pt,0,alpha)
     4                             + c1(ihq)
     5                             * SCLm0NC(jgrid,ixi(ihq)+1,
     6                             1,pt,0,alpha) )
                              CLps(ihq) = CLps(ihq) + as(pt)
     1                             * ( c0(ihq)
     2                             * SCLm0NC(jgrid,ixi(ihq),
     3                             2,pt,0,alpha)
     4                             + c1(ihq)
     5                             * SCLm0NC(jgrid,ixi(ihq)+1,
     6                             2,pt,0,alpha) )
                           enddo
                           C3nsm(ihq) = C3nsm(3)
                        endif
                     enddo
                  endif






















                  endif
*     
*     F2
*     
*     Light Component
*     
*     Singlet
                  OpF2(jgrid,3,1,alpha,beta)  = 
     1                 ( bq(1) + bq(2) + bq(3) )
     2                 * ( C2ps(3) + C2nsp(3) / 6d0 )
*     Gluon
                  OpF2(jgrid,3,2,alpha,beta)  = 
     1                 ( bq(1) + bq(2) + bq(3) ) * C2g(3)
*     T3
                  OpF2(jgrid,3,9,alpha,beta)  = 
     1                 fr3 * ( bq(2) - bq(1) ) * C2nsp(3) / 2d0
*     T8
                  OpF2(jgrid,3,10,alpha,beta) = 
     1                 ( bq(1) + bq(2) - 2d0 * bq(3) ) * C2nsp(3) / 6d0
*     T15
                  OpF2(jgrid,3,11,alpha,beta) = 
     1                 ( bq(1) + bq(2) + bq(3) ) * C2nsp(3) / 12d0
*     T24
                  OpF2(jgrid,3,12,alpha,beta) = 
     1                 ( bq(1) + bq(2) + bq(3) ) * C2nsp(3) / 20d0
*     T35
                  OpF2(jgrid,3,13,alpha,beta) = 
     1                 ( bq(1) + bq(2) + bq(3) ) * C2nsp(3) / 30d0
*     
*     Charm Component
*     
                  if(nf.ge.4)then
*     Singlet
                     OpF2(jgrid,4,1,alpha,beta)  = bq(4)
     1                    * ( C2ps(4) + C2nsp(4) / 6d0 )
*     Gluon
                     OpF2(jgrid,4,2,alpha,beta)  = bq(4) * C2g(4)
*     T15
                     OpF2(jgrid,4,11,alpha,beta) = - bq(4) * C2nsp(4)
     1                    / 4d0
*     T24
                     OpF2(jgrid,4,12,alpha,beta) = bq(4) * C2nsp(4)
     1                    / 20d0
*     T35
                     OpF2(jgrid,4,13,alpha,beta) = bq(4) * C2nsp(4)
     1                    / 30d0
                  endif
*     
*     Bottom Component
*     
                  if(nf.ge.5)then
*     Singlet
                     OpF2(jgrid,5,1,alpha,beta)  = bq(5)
     1                    * ( C2ps(5) + C2nsp(5) / 6d0 )
*     Gluon
                     OpF2(jgrid,5,2,alpha,beta)  = bq(5) * C2g(5)
*     T24
                     OpF2(jgrid,5,12,alpha,beta) = - bq(5) * C2nsp(5)
     1                    / 5d0
*     T35
                     OpF2(jgrid,5,13,alpha,beta) = bq(5) * C2nsp(5)
     1                    / 30d0
                  endif
*     
*     Top Component
*     
                  if(nf.ge.6)then
*     Singlet
                     OpF2(jgrid,6,1,alpha,beta)  = bq(6)
     1                    * ( C2ps(6) + C2nsp(6) / 6d0 )
*     Gluon
                     OpF2(jgrid,6,2,alpha,beta)  = bq(6) * C2g(6)
*     T35
                     OpF2(jgrid,6,13,alpha,beta) = - bq(6) * C2nsp(6)
     1                    / 6d0
                  endif
*     
*     Total
*     
                  do ipdf=0,13
                     do ihq=3,nf
                        OpF2(jgrid,7,ipdf,alpha,beta) = 
     1                       OpF2(jgrid,7,ipdf,alpha,beta)
     2                       + OpF2(jgrid,ihq,ipdf,alpha,beta)
                     enddo
                  enddo
*     
*     FL
*     
*     Light Component
*     
*     Singlet
                  OpFL(jgrid,3,1,alpha,beta)  = 
     1                 ( bq(1) + bq(2) + bq(3) )
     2                 * ( CLps(3) + CLnsp(3) / 6d0 )
*     Gluon
                  OpFL(jgrid,3,2,alpha,beta)  = 
     1                 ( bq(1) + bq(2) + bq(3) ) * CLg(3)
*     T3
                  OpFL(jgrid,3,9,alpha,beta)  = 
     1                 fr3 * ( bq(2) - bq(1) ) * CLnsp(3) / 2d0
*     T8
                  OpFL(jgrid,3,10,alpha,beta) = 
     1                 ( bq(1) + bq(2) - 2d0 * bq(3) ) * CLnsp(3) / 6d0
*     T15
                  OpFL(jgrid,3,11,alpha,beta) = 
     1                 ( bq(1) + bq(2) + bq(3) ) * CLnsp(3) / 12d0
*     T24
                  OpFL(jgrid,3,12,alpha,beta) = 
     1                 ( bq(1) + bq(2) + bq(3) ) * CLnsp(3) / 20d0
*     T35
                  OpFL(jgrid,3,13,alpha,beta) = 
     1                 ( bq(1) + bq(2) + bq(3) ) * CLnsp(3) / 30d0
*     
*     Charm Component
*     
                  if(nf.ge.4)then
*     Singlet
                     OpFL(jgrid,4,1,alpha,beta)  = bq(4)
     1                    * ( CLps(4) + CLnsp(4) / 6d0 )
*     Gluon
                     OpFL(jgrid,4,2,alpha,beta)  = bq(4) * CLg(4)
*     T15
                     OpFL(jgrid,4,11,alpha,beta) = - bq(4) * CLnsp(4)
     1                    / 4d0
*     T24
                     OpFL(jgrid,4,12,alpha,beta) = bq(4) * CLnsp(4)
     1                    / 20d0
*     T35
                     OpFL(jgrid,4,13,alpha,beta) = bq(4) * CLnsp(4)
     1                    / 30d0
                  endif
*     
*     Bottom Component
*     
                  if(nf.ge.5)then
*     Singlet
                     OpFL(jgrid,5,1,alpha,beta)  = bq(5)
     1                    * ( CLps(5) + CLnsp(5) / 6d0 )
*     Gluon
                     OpFL(jgrid,5,2,alpha,beta)  = bq(5) * CLg(5)
*     T24
                     OpFL(jgrid,5,12,alpha,beta) = - bq(5) * CLnsp(5)
     1                    / 5d0
*     T35
                     OpFL(jgrid,5,13,alpha,beta) = bq(5) * CLnsp(5)
     1                    / 30d0
                  endif
*     
*     Top Component
*     
                  if(nf.ge.6)then
*     Singlet
                     OpFL(jgrid,6,1,alpha,beta)  = bq(6)
     1                    * ( CLps(6) + CLnsp(6) / 6d0 )
*     Gluon
                     OpFL(jgrid,6,2,alpha,beta)  = bq(6) * CLg(6)
*     T35
                     OpFL(jgrid,6,13,alpha,beta) = - bq(6) * CLnsp(6)
     1                    / 6d0
                  endif
*     
*     Total
*     
                  do ipdf=0,13
                     do ihq=3,nf
                        OpFL(jgrid,7,ipdf,alpha,beta) = 
     1                       OpFL(jgrid,7,ipdf,alpha,beta)
     2                       + OpFL(jgrid,ihq,ipdf,alpha,beta)
                     enddo
                  enddo
*     
*     F3 (Only for neutral current processes)
*     
                  if(ProcessDIS.eq."NC")then
*     
*     Light Component
*     
*     Valence
                     OpF3(jgrid,3,3,alpha,beta) = 
     1                    ( dq(1) + dq(2) + dq(3) ) * C3nsm(3) / 6d0
*     V3
                     OpF3(jgrid,3,4,alpha,beta) = 
     1                    fr3 * ( dq(2) - dq(1) ) * C3nsm(3) / 2d0
*     V8
                     OpF3(jgrid,3,5,alpha,beta) = 
     1                    ( dq(1) + dq(2) - 2d0 * dq(3) ) * C3nsm(3)
     2                    / 6d0
*     V15
                     OpF3(jgrid,3,6,alpha,beta) = 
     1                    ( dq(1) + dq(2) + dq(3) ) * C3nsm(3) / 12d0
*     V24
                     OpF3(jgrid,3,7,alpha,beta) = 
     1                    ( dq(1) + dq(2) + dq(3) ) * C3nsm(3) / 20d0
*     V35
                     OpF3(jgrid,3,8,alpha,beta) = 
     1                    ( dq(1) + dq(2) + dq(3) ) * C3nsm(3) / 30d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Valence
                        OpF3(jgrid,4,3,alpha,beta) = dq(4) * C3nsm(4)
     1                       / 6d0
*     V15
                        OpF3(jgrid,4,6,alpha,beta) = - dq(4) * C3nsm(4)
     1                       / 4d0
*     V24
                        OpF3(jgrid,4,7,alpha,beta) = dq(4) * C3nsm(4)
     1                       / 20d0
*     V35
                        OpF3(jgrid,4,8,alpha,beta) = dq(4) * C3nsm(4)
     1                       / 30d0
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Valence
                        OpF3(jgrid,5,3,alpha,beta) = dq(5) * C3nsm(5)
     1                       / 6d0
*     V24
                        OpF3(jgrid,5,7,alpha,beta) = - dq(5) * C3nsm(5)
     1                       / 5d0
*     V35
                        OpF3(jgrid,5,8,alpha,beta) = dq(5) * C3nsm(5)
     1                       / 30d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Valence
                        OpF3(jgrid,6,3,alpha,beta) = dq(6) * C3nsm(6)
     1                       / 6d0
*     V35
                        OpF3(jgrid,6,8,alpha,beta) = - dq(6) * C3nsm(6)
     1                       / 6d0
                     endif
*     
*     Total
*     
                     do ipdf=0,13
                        do ihq=3,nf
                           OpF3(jgrid,7,ipdf,alpha,beta) = 
     1                          OpF3(jgrid,7,ipdf,alpha,beta)
     2                          + OpF3(jgrid,ihq,ipdf,alpha,beta)
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
*     
*     Charged current structure functions
*     
      elseif(ProcessDIS.eq."CC")then
*     
*     Different projectiles
*     
         if(ProjectileDIS(1:8).eq."neutrino".or.
     1        ProjectileDIS(1:8).eq."positron")then
            ipr = 1
         elseif(ProjectileDIS.eq."antineutrino".or.
     1           ProjectileDIS(1:8).eq."electron")then
            ipr = - 1
         endif
*     
*     Define useful constants
*     
         Kl = V_ud2 + V_us2
         Kc = V_cd2 + V_cs2
         Kb = V_ub2 + V_cb2
         Kt = V_td2 + V_ts2 + V_tb2
*     
         do jgrid=1,ngrid
*     
*     If an external grid is found compute the whole operator
*     
            gbound = 0
            if(IsExt(jgrid)) gbound = nin(jgrid)
*     
            do alpha=0,gbound
               do beta=alpha,nin(jgrid)
*     
*     Put together coefficient functions
*     
                  do ihq=3,6
                     C2g(ihq)   = 0d0
                     C2ps(ihq)  = 0d0
                     C2nsp(ihq) = 0d0
                     C2nsm(ihq) = 0d0
                     CLg(ihq)   = 0d0
                     CLps(ihq)  = 0d0
                     CLnsp(ihq) = 0d0
                     CLnsm(ihq) = 0d0
                     C3g(ihq)   = 0d0
                     C3nsp(ihq) = 0d0
                     C3nsm(ihq) = 0d0
                  enddo
*     
*     Construct coefficient functions according to the scheme chosen
*     
                  if(MassScheme.eq."ZM-VFNS")then
                     do pt=0,ipt
                        C2g(3)   = C2g(3)
     1                       + as(pt) * SC2zm(jgrid,nf,1,pt,alpha,beta)
                        C2ps(3)  = C2ps(3)
     1                       + as(pt) * SC2zm(jgrid,nf,2,pt,alpha,beta)
                        C2nsp(3) = C2nsp(3)
     1                       + as(pt) * SC2zm(jgrid,nf,3,pt,alpha,beta)
                        C2nsm(3) = C2nsm(3)
     1                       + as(pt) * SC2zm(jgrid,nf,4,pt,alpha,beta)
                        CLg(3)   = CLg(3)
     1                       + as(pt) * SCLzm(jgrid,nf,1,pt,alpha,beta)
                        CLps(3)  = CLps(3)
     1                       + as(pt) * SCLzm(jgrid,nf,2,pt,alpha,beta)
                        CLnsp(3) = CLnsp(3)
     1                       + as(pt) * SCLzm(jgrid,nf,3,pt,alpha,beta)
                        CLnsm(3) = CLnsm(3)
     1                       + as(pt) * SCLzm(jgrid,nf,4,pt,alpha,beta)
                        C3nsp(3) = C3nsp(3)
     1                       + as(pt) * SC3zm(jgrid,nf,3,pt,alpha,beta)
                        C3nsm(3) = C3nsm(3)
     1                       + as(pt) * SC3zm(jgrid,nf,4,pt,alpha,beta)
                     enddo
                     if(nf.gt.3)then
                        do ihq=4,nf
                           C2g(ihq)   = C2g(3)
                           C2ps(ihq)  = C2ps(3)
                           C2nsp(ihq) = C2nsp(3)
                           C2nsm(ihq) = C2nsm(3)
                           CLg(ihq)   = CLg(3)
                           CLps(ihq)  = CLps(3)
                           CLnsp(ihq) = CLnsp(3)
                           CLnsm(ihq) = CLnsm(3)
                           C3nsp(ihq) = C3nsp(3)
                           C3nsm(ihq) = C3nsm(3)
                        enddo
                     endif
                  endif
*
*     Change sign to the F3 coefficient functions according to the projectile
*
                  C3g   = ipr * C3g
                  C3nsp = ipr * C3nsp
                  C3nsm = ipr * C3nsm
*     
*     F2
*     
*     Light Component
*     
*     Singlet
                  OpF2(jgrid,3,1,alpha,beta)  = 
     1                 2d0 * Kl * ( C2ps(3) + C2nsp(3) / 6d0 )
*     Gluon
                  OpF2(jgrid,3,2,alpha,beta)  = 2d0 * Kl * C2g(3)
*     V3
                  OpF2(jgrid,3,4,alpha,beta)  = - ipr * fr3
     1                 * ( 2d0 * V_ud2 + V_us2 ) * C2nsm(3) / 2d0
*     V8
                  OpF2(jgrid,3,5,alpha,beta)  = - ipr
     1                 * V_us2 * C2nsm(3) / 2d0
*     T3
                  OpF2(jgrid,3,9,alpha,beta)  = fr3
     1                 * V_us2 * C2nsp(3) / 2d0
*     T8
                  OpF2(jgrid,3,10,alpha,beta) = 
     1                 ( 2d0 * V_ud2 - V_us2 ) * C2nsp(3) / 6d0
*     T15
                  OpF2(jgrid,3,11,alpha,beta) = Kl * C2nsp(3) / 6d0
*     T24
                  OpF2(jgrid,3,12,alpha,beta) = Kl * C2nsp(3) / 10d0
*     T35
                  OpF2(jgrid,3,13,alpha,beta) = Kl * C2nsp(3) / 15d0
*     
*     Charm Component
*     
                  if(nf.ge.4)then
*     Singlet
                     OpF2(jgrid,4,1,alpha,beta)  = 
     1                    2d0 * Kc * ( C2ps(4) + C2nsp(4) / 6d0 )
*     Gluon
                     OpF2(jgrid,4,2,alpha,beta)  = 2d0 * Kc * C2g(4)
*     V3
                     OpF2(jgrid,4,4,alpha,beta)  = 
     1                    - ipr * fr3 * V_cd2 * C2nsm(4) / 2d0
*     V8
                     OpF2(jgrid,4,5,alpha,beta)  =
     1                    ipr * ( V_cd2 - 2d0 * V_cs2 ) * C2nsm(4) / 6d0
*     V15
                     OpF2(jgrid,4,6,alpha,beta)  =
     1                    ipr * Kc * C2nsm(4) / 3d0
*     T3
                     OpF2(jgrid,4,9,alpha,beta)  = 
     1                    - fr3 * V_cd2 * C2nsp(4) / 2d0
*     T8
                     OpF2(jgrid,4,10,alpha,beta) = 
     1                    ( V_cd2 - 2d0 * V_cs2 ) * C2nsp(4) / 6d0
*     T15
                     OpF2(jgrid,4,11,alpha,beta) = - Kc * C2nsp(4) / 6d0
*     T24
                     OpF2(jgrid,4,12,alpha,beta) = Kc * C2nsp(4) / 10d0
*     T35
                     OpF2(jgrid,4,13,alpha,beta) = Kc * C2nsp(4) / 15d0
                  endif
*     
*     Bottom Component
*     
                  if(nf.ge.5)then
*     Singlet
                     OpF2(jgrid,5,1,alpha,beta)  = 
     1                    2d0 * Kb * ( C2ps(5) + C2nsp(5) / 6d0 )
*     Gluon
                     OpF2(jgrid,5,2,alpha,beta)  = 2d0 * Kb * C2g(5)
*     V3
                     OpF2(jgrid,5,4,alpha,beta)  =
     1                    - ipr * fr3 * V_ub2 * C2nsm(5) / 2d0
*     V8
                     OpF2(jgrid,5,5,alpha,beta)  = 
     1                    - ipr * V_ub2 * C2nsm(5) / 6d0
*     V15
                     OpF2(jgrid,5,6,alpha,beta)  = - ipr
     1                    * ( V_ub2 - 3d0 * V_cb2 ) * C2nsm(5) / 12d0
*     V24
                     OpF2(jgrid,5,7,alpha,beta)  = 
     1                    - ipr * Kb * C2nsm(5) / 4d0
*     T3
                     OpF2(jgrid,5,9,alpha,beta)  =
     1                    fr3 * V_ub2 * C2nsp(5) / 2d0 
*     T8
                     OpF2(jgrid,5,10,alpha,beta) = 
     1                    V_ub2 * C2nsp(5) / 6d0
*     T15
                     OpF2(jgrid,5,11,alpha,beta) = 
     1                    ( V_ub2 - 3d0 * V_cb2 ) * C2nsp(5) / 12d0
*     T24
                     OpF2(jgrid,5,12,alpha,beta) = 
     1                    - 3d0 * Kb * C2nsp(5) / 20d0
*     T35
                     OpF2(jgrid,5,13,alpha,beta) = Kb * C2nsp(5) / 15d0
                  endif
*     
*     Top Component
*     
                  if(nf.ge.6)then
*     Singlet
                     OpF2(jgrid,6,1,alpha,beta)  = 
     1                    2d0 * Kt * ( C2ps(6) + C2nsp(6) / 6d0 )
*     Gluon
                     OpF2(jgrid,6,2,alpha,beta)  = 2d0 * Kt * C2g(6)
*     V3
                     OpF2(jgrid,6,4,alpha,beta)  =
     1                    - ipr * fr3 * V_td2 * C2nsm(6) / 2d0
*     V8
                     OpF2(jgrid,6,5,alpha,beta)  = 
     1                    ipr * ( V_td2 - 2d0 * V_ts2 ) * C2nsm(6) / 6d0
*     V15
                     OpF2(jgrid,6,6,alpha,beta)  =
     1                    ipr * ( V_td2 + V_ts2 ) * C2nsm(6) / 12d0
*     V24
                     OpF2(jgrid,6,7,alpha,beta)  = 
     1                    ipr * ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                    * C2nsm(6) / 20d0
*     V35
                     OpF2(jgrid,6,8,alpha,beta)  =
     1                    ipr * Kt * C2nsm(6) / 5d0
*     T3
                     OpF2(jgrid,6,9,alpha,beta)  =
     1                    - fr3 * V_td2 * C2nsp(6) / 2d0 
*     T8
                     OpF2(jgrid,6,10,alpha,beta) =
     1                    ( V_td2 - 2d0 * V_ts2 ) * C2nsp(6) / 6d0
*     T15
                     OpF2(jgrid,6,11,alpha,beta) = 
     1                    ( V_td2 + V_ts2 ) * C2nsp(6) / 12d0
*     T24
                     OpF2(jgrid,6,12,alpha,beta) = 
     1                    ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                    * C2nsp(6) / 20d0
*     T35
                     OpF2(jgrid,6,13,alpha,beta) =
     1                    - 2d0 * Kt * C2nsp(6) / 15d0
                  endif
*     
*     Total
*     
                  do ipdf=0,13
                     do ihq=3,nf
                        OpF2(jgrid,7,ipdf,alpha,beta) = 
     1                       OpF2(jgrid,7,ipdf,alpha,beta)
     2                       + OpF2(jgrid,ihq,ipdf,alpha,beta)
                     enddo
                  enddo
*     
*     FL
*     
*     Light Component
*     
*     Singlet
                  OpFL(jgrid,3,1,alpha,beta)  = 
     1                 2d0 * Kl * ( CLps(3) + CLnsp(3) / 6d0 )
*     Gluon
                  OpFL(jgrid,3,2,alpha,beta)  = 2d0 * Kl * CLg(3)
*     V3
                  OpFL(jgrid,3,4,alpha,beta)  = - ipr * fr3
     1                 * ( 2d0 * V_ud2 + V_us2 ) * CLnsm(3) / 2d0
*     V8
                  OpFL(jgrid,3,5,alpha,beta)  = - ipr
     1                 * V_us2 * CLnsm(3) / 2d0
*     T3
                  OpFL(jgrid,3,9,alpha,beta)  = fr3
     1                 * V_us2 * CLnsp(3) / 2d0
*     T8
                  OpFL(jgrid,3,10,alpha,beta) = 
     1                 ( 2d0 * V_ud2 - V_us2 ) * CLnsp(3) / 6d0
*     T15
                  OpFL(jgrid,3,11,alpha,beta) = Kl * CLnsp(3) / 6d0
*     T24
                  OpFL(jgrid,3,12,alpha,beta) = Kl * CLnsp(3) / 10d0
*     T35
                  OpFL(jgrid,3,13,alpha,beta) = Kl * CLnsp(3) / 15d0
*     
*     Charm Component
*     
                  if(nf.ge.4)then
*     Singlet
                     OpFL(jgrid,4,1,alpha,beta)  = 
     1                    2d0 * Kc * ( CLps(4) + CLnsp(4) / 6d0 )
*     Gluon
                     OpFL(jgrid,4,2,alpha,beta)  = 2d0 * Kc * CLg(4)
*     V3
                     OpFL(jgrid,4,4,alpha,beta)  = 
     1                    - ipr * fr3 * V_cd2 * CLnsm(4) / 2d0
*     V8
                     OpFL(jgrid,4,5,alpha,beta)  =
     1                    ipr * ( V_cd2 - 2d0 * V_cs2 ) * CLnsm(4) / 6d0
*     V15
                     OpFL(jgrid,4,6,alpha,beta)  =
     1                    ipr * Kc * CLnsm(4) / 3d0
*     T3
                     OpFL(jgrid,4,9,alpha,beta)  = 
     1                    - fr3 * V_cd2 * CLnsp(4) / 2d0
*     T8
                     OpFL(jgrid,4,10,alpha,beta) = 
     1                    ( V_cd2 - 2d0 * V_cs2 ) * CLnsp(4) / 6d0
*     T15
                     OpFL(jgrid,4,11,alpha,beta) = - Kc * CLnsp(4) / 6d0
*     T24
                     OpFL(jgrid,4,12,alpha,beta) = Kc * CLnsp(4) / 10d0
*     T35
                     OpFL(jgrid,4,13,alpha,beta) = Kc * CLnsp(4) / 15d0
                  endif
*     
*     Bottom Component
*     
                  if(nf.ge.5)then
*     Singlet
                     OpFL(jgrid,5,1,alpha,beta)  = 
     1                    2d0 * Kb * ( CLps(5) + CLnsp(5) / 6d0 )
*     Gluon
                     OpFL(jgrid,5,2,alpha,beta)  = 2d0 * Kb * CLg(5)
*     V3
                     OpFL(jgrid,5,4,alpha,beta)  =
     1                    - ipr * fr3 * V_ub2 * CLnsm(5) / 2d0
*     V8
                     OpFL(jgrid,5,5,alpha,beta)  = 
     1                    - ipr * V_ub2 * CLnsm(5) / 6d0
*     V15
                     OpFL(jgrid,5,6,alpha,beta)  = - ipr
     1                    * ( V_ub2 - 3d0 * V_cb2 ) * CLnsm(5) / 12d0
*     V24
                     OpFL(jgrid,5,7,alpha,beta)  = 
     1                    - ipr * Kb * CLnsm(5) / 4d0
*     T3
                     OpFL(jgrid,5,9,alpha,beta)  =
     1                    fr3 * V_ub2 * CLnsp(5) / 2d0 
*     T8
                     OpFL(jgrid,5,10,alpha,beta) = 
     1                    V_ub2 * CLnsp(5) / 6d0
*     T15
                     OpFL(jgrid,5,11,alpha,beta) = 
     1                    ( V_ub2 - 3d0 * V_cb2 ) * CLnsp(5) / 12d0
*     T24
                     OpFL(jgrid,5,12,alpha,beta) = 
     1                    - 3d0 * Kb * CLnsp(5) / 20d0
*     T35
                     OpFL(jgrid,5,13,alpha,beta) = Kb * CLnsp(5) / 15d0
                  endif
*     
*     Top Component
*     
                  if(nf.ge.6)then
*     Singlet
                     OpFL(jgrid,6,1,alpha,beta)  = 
     1                    2d0 * Kt * ( CLps(6) + CLnsp(6) / 6d0 )
*     Gluon
                     OpFL(jgrid,6,2,alpha,beta)  = 2d0 * Kt * CLg(6)
*     V3
                     OpFL(jgrid,6,4,alpha,beta)  =
     1                    - ipr * fr3 * V_td2 * CLnsm(6) / 2d0
*     V8
                     OpFL(jgrid,6,5,alpha,beta)  = 
     1                    ipr * ( V_td2 - 2d0 * V_ts2 ) * CLnsm(6) / 6d0
*     V15
                     OpFL(jgrid,6,6,alpha,beta)  =
     1                    ipr * ( V_td2 + V_ts2 ) * CLnsm(6) / 12d0
*     V24
                     OpFL(jgrid,6,7,alpha,beta)  = 
     1                    ipr * ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                    * CLnsm(6) / 20d0
*     V35
                     OpFL(jgrid,6,8,alpha,beta)  =
     1                    ipr * Kt * CLnsm(6) / 5d0
*     T3
                     OpFL(jgrid,6,9,alpha,beta)  =
     1                    - fr3 * V_td2 * CLnsp(6) / 2d0 
*     T8
                     OpFL(jgrid,6,10,alpha,beta) =
     1                    ( V_td2 - 2d0 * V_ts2 ) * CLnsp(6) / 6d0
*     T15
                     OpFL(jgrid,6,11,alpha,beta) = 
     1                    ( V_td2 + V_ts2 ) * CLnsp(6) / 12d0
*     T24
                     OpFL(jgrid,6,12,alpha,beta) = 
     1                    ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                    * CLnsp(6) / 20d0
*     T35
                     OpFL(jgrid,6,13,alpha,beta) =
     1                    - 2d0 * Kt * CLnsp(6) / 15d0
                  endif
*     
*     Total
*     
                  do ipdf=0,13
                     do ihq=3,nf
                        OpFL(jgrid,7,ipdf,alpha,beta) = 
     1                       OpFL(jgrid,7,ipdf,alpha,beta)
     2                       + OpFL(jgrid,ihq,ipdf,alpha,beta)
                     enddo
                  enddo
*     
*     F3
*     
*     Light Component
*     
*     Gluon
                  OpF3(jgrid,3,2,alpha,beta)  = 2d0 * Kl * C3g(3)
*     Valence
                  OpF3(jgrid,3,3,alpha,beta)  = Kl * C3nsm(3) / 3d0
*     V3
                  OpF3(jgrid,3,4,alpha,beta)  = fr3
     1                 * V_us2 * C3nsp(3) / 2d0
*     V8
                  OpF3(jgrid,3,5,alpha,beta)  = 
     1                 ( 2d0 * V_ud2 - V_us2 ) * C3nsm(3) / 6d0
*     V15
                  OpF3(jgrid,3,6,alpha,beta)  = Kl * C3nsm(3) / 6d0
*     V24
                  OpF3(jgrid,3,7,alpha,beta)  = Kl * C3nsm(3) / 10d0
*     V35
                  OpF3(jgrid,3,8,alpha,beta)  = Kl * C3nsm(3) / 15d0
*     T3
                  OpF3(jgrid,3,9,alpha,beta)  = - ipr * fr3
     1                 * ( 2d0 * V_ud2 + V_us2 ) * C3nsp(3) / 2d0
*     T8
                  OpF3(jgrid,3,10,alpha,beta) = - ipr
     1                 * V_us2 * C3nsp(3) / 2d0
*     
*     Charm Component
*     
                  if(nf.ge.4)then
*     Gluon
                     OpF3(jgrid,4,2,alpha,beta)  = 2d0 * Kc * C3g(4)
*     Valence
                     OpF3(jgrid,4,3,alpha,beta)  = Kc * C3nsm(4) / 3d0
*     V3
                     OpF3(jgrid,4,4,alpha,beta)  = 
     1                    - fr3 * V_cd2 * C3nsm(4) / 2d0
*     V8
                     OpF3(jgrid,4,5,alpha,beta)  = 
     1                    ( V_cd2 - 2d0 * V_cs2 ) * C3nsm(4) / 6d0
*     V15
                     OpF3(jgrid,4,6,alpha,beta)  = - Kc * C3nsm(4) / 6d0
*     V24
                     OpF3(jgrid,4,7,alpha,beta)  = Kc * C3nsm(4) / 10d0
*     V35
                     OpF3(jgrid,4,8,alpha,beta)  = Kc * C3nsm(4) / 15d0
*     T3
                     OpF3(jgrid,4,9,alpha,beta)  = 
     1                    - ipr * fr3 * V_cd2 * C3nsp(4) / 2d0
*     T8
                     OpF3(jgrid,4,10,alpha,beta) =
     1                    ipr * ( V_cd2 - 2d0 * V_cs2 ) * C3nsp(4) / 6d0
*     T15
                     OpF3(jgrid,4,11,alpha,beta) =
     1                    ipr * Kc * C3nsp(4) / 3d0
                  endif
*     
*     Bottom Component
*     
                  if(nf.ge.5)then
*     Gluon
                     OpF3(jgrid,5,2,alpha,beta)  = 2d0 * Kb * C3g(5)
*     Valence
                     OpF3(jgrid,5,3,alpha,beta)  = Kb * C3nsm(5) / 3d0
*     V3
                     OpF3(jgrid,5,4,alpha,beta)  =
     1                    fr3 * V_ub2 * C3nsm(5) / 2d0 
*     V8
                     OpF3(jgrid,5,5,alpha,beta)  = 
     1                    V_ub2 * C3nsm(5) / 6d0
*     V15
                     OpF3(jgrid,5,6,alpha,beta)  = 
     1                    ( V_ub2 - 3d0 * V_cb2 ) * C3nsm(5) / 12d0
*     V24
                     OpF3(jgrid,5,7,alpha,beta)  = 
     1                    - 3d0 * Kb * C3nsm(5) / 20d0
*     V35
                     OpF3(jgrid,5,8,alpha,beta)  = Kb * C3nsm(5) / 15d0
*     T3
                     OpF3(jgrid,5,9,alpha,beta)  =
     1                    - ipr * fr3 * V_ub2 * C3nsp(5) / 2d0
*     T8
                     OpF3(jgrid,5,10,alpha,beta) = 
     1                    - ipr * V_ub2 * C3nsp(5) / 6d0
*     T15
                     OpF3(jgrid,5,11,alpha,beta) = - ipr
     1                    * ( V_ub2 - 3d0 * V_cb2 ) * C3nsp(5) / 12d0
*     T24
                     OpF3(jgrid,5,12,alpha,beta) = 
     1                    - ipr * Kb * C3nsp(5) / 4d0
                  endif
*     
*     Top Component
*     
                  if(nf.ge.6)then
*     Gluon
                     OpF3(jgrid,6,2,alpha,beta)  = 2d0 * Kt * C3g(6)
*     Valence
                     OpF3(jgrid,6,3,alpha,beta)  = Kt * C3nsm(6) / 3d0
*     V3
                     OpF3(jgrid,6,4,alpha,beta)  =
     1                    - fr3 * V_td2 * C3nsp(6) / 2d0 
*     V8
                     OpF3(jgrid,6,5,alpha,beta)  =
     1                    ( V_td2 - 2d0 * V_ts2 ) * C3nsp(6) / 6d0
*     V15
                     OpF3(jgrid,6,6,alpha,beta)  = 
     1                    ( V_td2 + V_ts2 ) * C3nsp(6) / 12d0
*     V24
                     OpF3(jgrid,6,7,alpha,beta)  = 
     1                    ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                    * C3nsp(6) / 20d0
*     V35
                     OpF3(jgrid,6,8,alpha,beta)  =
     1                    - 2d0 * Kt * C3nsp(6) / 15d0
*     T3
                     OpF3(jgrid,6,9,alpha,beta)  =
     1                    - ipr * fr3 * V_td2 * C3nsp(6) / 2d0
*     T8
                     OpF3(jgrid,6,10,alpha,beta) = 
     1                    ipr * ( V_td2 - 2d0 * V_ts2 ) * C3nsp(6) / 6d0
*     T15
                     OpF3(jgrid,6,11,alpha,beta) =
     1                    ipr * ( V_td2 + V_ts2 ) * C3nsp(6) / 12d0
*     T24
                     OpF3(jgrid,6,12,alpha,beta) = 
     1                    ipr * ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                    * C3nsp(6) / 20d0
*     T35
                     OpF3(jgrid,6,13,alpha,beta) =
     1                    ipr * Kt * C3nsp(6) / 5d0
                  endif
*     
*     Total
*     
                  do ipdf=0,13
                     do ihq=3,nf
                        OpF3(jgrid,7,ipdf,alpha,beta) = 
     1                       OpF3(jgrid,7,ipdf,alpha,beta)
     2                       + OpF3(jgrid,ihq,ipdf,alpha,beta)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
*
      call cpu_time(t2)
*
      write(6,"(a,a,f9.5,a)") " Computation of the DIS operators",
     1                        " completed in",t2-t1," s"
      write(6,*) " "
*
      return
      end 