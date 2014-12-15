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
c      include "../commons/Nf_FF.h"
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
      double precision Q2
      double precision ipr
      double precision as(0:2),a_QCD
      double precision bq(6),dq(6)
      double precision frac,fr3
      double precision C2g,C2ps,C2nsp,C2nsm
      double precision CLg,CLps,CLnsp,CLnsm
      double precision C3nsp,C3nsm
      double precision t1,t2



c      integer ixi(4:6)
c      double precision W2,xi(4:6),c0(4:6),c1(4:6)
c      double precision singlet
c      double precision F2t,FLt,F3t
c      double precision fup,fub,fdw,fdb
c      double precision diff(nxi),sgn
c      double precision SC2(0:ngrid_max,6,3,0:2,0:nint_max,0:nint_max)
c      double precision SCL(0:ngrid_max,6,3,0:2,0:nint_max,0:nint_max)
c      double precision SC3(0:ngrid_max,6,3,0:2,0:nint_max,0:nint_max)
c      double precision damp(4:6)


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
*     Construct the coefficient functions according to the scheme chosen
*
      if(MassScheme.eq."ZM-VFNS")then
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
*     F2
*
*     Put together coefficient functions
*
                     C2g   = 0d0
                     C2ps  = 0d0
                     C2nsp = 0d0
                     do pt=0,ipt
                        C2g   = C2g
     1                       + as(pt) * SC2zm(jgrid,nf,1,pt,alpha,beta)
                        C2ps  = C2ps
     1                       + as(pt) * SC2zm(jgrid,nf,2,pt,alpha,beta)
                        C2nsp = C2nsp
     1                       + as(pt) * SC2zm(jgrid,nf,3,pt,alpha,beta)
                     enddo
*     
*     Light Component
*     
*     Singlet
                     OpF2(jgrid,3,1,alpha,beta)  = 
     1                    ( bq(1) + bq(2) + bq(3) )
     2                    * ( C2ps + C2nsp / 6d0 )
*     Gluon
                     OpF2(jgrid,3,2,alpha,beta)  = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2g
*     T3
                     OpF2(jgrid,3,9,alpha,beta)  = 
     1                    fr3 * ( bq(2) - bq(1) ) * C2nsp / 2d0
*     T8
                     OpF2(jgrid,3,10,alpha,beta) = 
     1                    ( bq(1) + bq(2) - 2d0 * bq(3) ) * C2nsp / 6d0
*     T15
                     OpF2(jgrid,3,11,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2nsp / 12d0
*     T24
                     OpF2(jgrid,3,12,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2nsp / 20d0
*     T35
                     OpF2(jgrid,3,13,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2nsp / 30d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Singlet
                        OpF2(jgrid,4,1,alpha,beta)  = bq(4)
     1                       * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,4,2,alpha,beta)  = bq(4) * C2g
*     T15
                        OpF2(jgrid,4,11,alpha,beta) = - bq(4) * C2nsp 
     1                                              / 4d0
*     T24
                        OpF2(jgrid,4,12,alpha,beta) = bq(4) * C2nsp
     1                                              / 20d0
*     T35
                        OpF2(jgrid,4,13,alpha,beta) = bq(4) * C2nsp
     1                                              / 30d0
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Singlet
                        OpF2(jgrid,5,1,alpha,beta)  = bq(5)
     1                       * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,5,2,alpha,beta)  = bq(5) * C2g
*     T24
                        OpF2(jgrid,5,12,alpha,beta) = - bq(5) * C2nsp
     1                                              / 5d0
*     T35
                        OpF2(jgrid,5,13,alpha,beta) = bq(5) * C2nsp
     1                                              / 30d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Singlet
                        OpF2(jgrid,6,1,alpha,beta)  = bq(6)
     1                       * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,6,2,alpha,beta)  = bq(6) * C2g
*     T35
                        OpF2(jgrid,6,13,alpha,beta) = - bq(6) * C2nsp
     1                                              / 6d0
                     endif
*     
*     Total
*     
                     do ipdf=0,13
                        do ihq=3,nf
                           OpF2(jgrid,7,ipdf,alpha,beta) = 
     1                          OpF2(jgrid,7,ipdf,alpha,beta)
     2                          + OpF2(jgrid,ihq,ipdf,alpha,beta)
                        enddo
                     enddo
*     
*     FL
*     
*     Put together coefficient functions
*     
                     CLg   = 0d0
                     CLps  = 0d0
                     CLnsp = 0d0
                     do pt=0,ipt
                        CLg   = CLg
     1                       + as(pt) * SCLzm(jgrid,nf,1,pt,alpha,beta)
                        CLps  = CLps
     1                       + as(pt) * SCLzm(jgrid,nf,2,pt,alpha,beta)
                        CLnsp = CLnsp
     1                       + as(pt) * SCLzm(jgrid,nf,3,pt,alpha,beta)
                     enddo
*     
*     Light Component
*     
*     Singlet
                     OpFL(jgrid,3,1,alpha,beta)  = 
     1                    ( bq(1) + bq(2) + bq(3) )
     2                    * ( CLps + CLnsp / 6d0 )
*     Gluon
                     OpFL(jgrid,3,2,alpha,beta)  = 
     1                    ( bq(1) + bq(2) + bq(3) ) * CLg
*     T3
                     OpFL(jgrid,3,9,alpha,beta)  = 
     1                    fr3 * ( bq(2) - bq(1) ) * CLnsp / 2d0
*     T8
                     OpFL(jgrid,3,10,alpha,beta) = 
     1                    ( bq(1) + bq(2) - 2d0 * bq(3) ) * CLnsp / 6d0
*     T15
                     OpFL(jgrid,3,11,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * CLnsp / 12d0
*     T24
                     OpFL(jgrid,3,12,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * CLnsp / 20d0
*     T35
                     OpFL(jgrid,3,13,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * CLnsp / 30d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Singlet
                        OpFL(jgrid,4,1,alpha,beta)  = bq(4)
     1                       * ( CLps + CLnsp / 6d0 )
*     Gluon
                        OpFL(jgrid,4,2,alpha,beta)  = bq(4) * CLg
*     T15
                        OpFL(jgrid,4,11,alpha,beta) = - bq(4) * CLnsp
     1                                              / 4d0
*     T24
                        OpFL(jgrid,4,12,alpha,beta) = bq(4) * CLnsp
     1                                              / 20d0
*     T35
                        OpFL(jgrid,4,13,alpha,beta) = bq(4) * CLnsp
     1                                              / 30d0
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Singlet
                        OpFL(jgrid,5,1,alpha,beta)  = bq(5)
     1                       * ( CLps + CLnsp / 6d0 )
*     Gluon
                        OpFL(jgrid,5,2,alpha,beta)  = bq(5) * CLg
*     T24
                        OpFL(jgrid,5,12,alpha,beta) = - bq(5) * CLnsp
     1                                              / 5d0
*     T35
                        OpFL(jgrid,5,13,alpha,beta) = bq(5) * CLnsp
     1                                              / 30d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Singlet
                        OpFL(jgrid,6,1,alpha,beta)  = bq(6)
     1                       * ( CLps + CLnsp / 6d0 )
*     Gluon
                        OpFL(jgrid,6,2,alpha,beta)  = bq(6) * CLg
*     T35
                        OpFL(jgrid,6,13,alpha,beta) = - bq(6) * CLnsp
     1                                              / 6d0
                     endif
*     
*     Total
*     
                     do ipdf=0,13
                        do ihq=3,nf
                           OpFL(jgrid,7,ipdf,alpha,beta) = 
     1                          OpFL(jgrid,7,ipdf,alpha,beta)
     2                          + OpFL(jgrid,ihq,ipdf,alpha,beta)
                        enddo
                     enddo
*     
*     F3 (Only for neutral current processes)
*     
                     if(ProcessDIS.eq."NC")then
*     
*     Put together coefficient functions
*     
                        C3nsm = 0d0
                        do pt=0,ipt
                           C3nsm = C3nsm
     1                        + as(pt) * SC3zm(jgrid,nf,4,pt,alpha,beta)
                        enddo
*     
*     Light Component
*     
*     Valence
                        OpF3(jgrid,3,3,alpha,beta) = 
     1                       ( dq(1) + dq(2) + dq(3) ) * C3nsm / 6d0
*     V3
                        OpF3(jgrid,3,4,alpha,beta) = 
     1                       fr3 * ( dq(2) - dq(1) ) * C3nsm / 2d0
*     V8
                        OpF3(jgrid,3,5,alpha,beta) = 
     1                       ( dq(1) + dq(2) - 2d0 * dq(3) ) * C3nsm
     2                       / 6d0
*     V15
                        OpF3(jgrid,3,6,alpha,beta) = 
     1                       ( dq(1) + dq(2) + dq(3) ) * C3nsm / 12d0
*     V24
                        OpF3(jgrid,3,7,alpha,beta) = 
     1                       ( dq(1) + dq(2) + dq(3) ) * C3nsm / 20d0
*     V35
                        OpF3(jgrid,3,8,alpha,beta) = 
     1                       ( dq(1) + dq(2) + dq(3) ) * C3nsm / 30d0
*     
*     Charm Component
*     
                        if(nf.ge.4)then
*     Valence
                           OpF3(jgrid,4,3,alpha,beta) = dq(4) * C3nsm
     1                                                / 6d0
*     V15
                           OpF3(jgrid,4,6,alpha,beta) = - dq(4) * C3nsm
     1                                                / 4d0
*     V24
                           OpF3(jgrid,4,7,alpha,beta) = dq(4) * C3nsm
     1                                                / 20d0
*     V35
                           OpF3(jgrid,4,8,alpha,beta) = dq(4) * C3nsm
     1                                                / 30d0
                        endif
*     
*     Bottom Component
*     
                        if(nf.ge.5)then
*     Valence
                           OpF3(jgrid,5,5,alpha,beta) = dq(5) * C3nsm
     1                                                / 6d0
*     V24
                           OpF3(jgrid,5,7,alpha,beta) = - dq(5) * C3nsm
     1                                                / 5d0
*     V35
                           OpF3(jgrid,5,8,alpha,beta) = dq(5) * C3nsm
     1                                                / 30d0
                        endif
*     
*     Top Component
*     
                        if(nf.ge.6)then
*     Valence
                           OpF3(jgrid,6,5,alpha,beta) = dq(6) * C3nsm
     1                                                / 6d0
*     V35
                           OpF3(jgrid,6,8,alpha,beta) = - dq(6) * C3nsm
     1                                                / 6d0
                        endif
*     
*     Total
*     
                        do ipdf=0,13
                           do ihq=3,nf
                              OpF3(jgrid,7,ipdf,alpha,beta) = 
     1                             OpF3(jgrid,7,ipdf,alpha,beta)
     2                             + OpF3(jgrid,ihq,ipdf,alpha,beta)
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
     1         ProjectileDIS(1:8).eq."positron")then
               ipr = 1d0
            elseif(ProjectileDIS.eq."antineutrino".or.
     1             ProjectileDIS(1:8).eq."electron")then
               ipr = - 1d0
            endif
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
*     F2
*
*     Put together coefficient functions
*
                     C2g   = 0d0
                     C2ps  = 0d0
                     C2nsp = 0d0
                     do pt=0,ipt
                        C2g   = C2g
     1                       + as(pt) * SC2zm(jgrid,nf,1,pt,alpha,beta)
                        C2ps  = C2ps
     1                       + as(pt) * SC2zm(jgrid,nf,2,pt,alpha,beta)
                        C2nsp = C2nsp
     1                       + as(pt) * SC2zm(jgrid,nf,3,pt,alpha,beta)
                     enddo
*     
*     Light Component
*     
*     Singlet
                     OpF2(jgrid,3,1,alpha,beta)  = 
     1                    ( bq(1) + bq(2) + bq(3) )
     2                    * ( C2ps + C2nsp / 6d0 )
*     Gluon
                     OpF2(jgrid,3,2,alpha,beta)  = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2g
*     T3
                     OpF2(jgrid,3,9,alpha,beta)  = 
     1                    fr3 * ( bq(2) - bq(1) ) * C2nsp / 2d0
*     T8
                     OpF2(jgrid,3,10,alpha,beta) = 
     1                    ( bq(1) + bq(2) - 2d0 * bq(3) ) * C2nsp / 6d0
*     T15
                     OpF2(jgrid,3,11,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2nsp / 12d0
*     T24
                     OpF2(jgrid,3,12,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2nsp / 20d0
*     T35
                     OpF2(jgrid,3,13,alpha,beta) = 
     1                    ( bq(1) + bq(2) + bq(3) ) * C2nsp / 30d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Singlet
                        OpF2(jgrid,4,1,alpha,beta)  = bq(4)
     1                       * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,4,2,alpha,beta)  = bq(4) * C2g
*     T15
                        OpF2(jgrid,4,11,alpha,beta) = - bq(4) * C2nsp 
     1                                              / 4d0
*     T24
                        OpF2(jgrid,4,12,alpha,beta) = bq(4) * C2nsp
     1                                              / 20d0
*     T35
                        OpF2(jgrid,4,13,alpha,beta) = bq(4) * C2nsp
     1                                              / 30d0
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Singlet
                        OpF2(jgrid,5,1,alpha,beta)  = bq(5)
     1                       * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,5,2,alpha,beta)  = bq(5) * C2g
*     T24
                        OpF2(jgrid,5,12,alpha,beta) = - bq(5) * C2nsp
     1                                              / 5d0
*     T35
                        OpF2(jgrid,5,13,alpha,beta) = bq(5) * C2nsp
     1                                              / 30d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Singlet
                        OpF2(jgrid,6,1,alpha,beta)  = bq(6)
     1                       * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,6,2,alpha,beta)  = bq(6) * C2g
*     T35
                        OpF2(jgrid,6,13,alpha,beta) = - bq(6) * C2nsp
     1                                              / 6d0
                     endif
*     
*     Total
*     
                     do ipdf=0,13
                        do ihq=3,nf
                           OpF2(jgrid,7,ipdf,alpha,beta) = 
     1                          OpF2(jgrid,7,ipdf,alpha,beta)
     2                          + OpF2(jgrid,ihq,ipdf,alpha,beta)
                        enddo
                     enddo





                  enddo
               enddo
            enddo
         endif



      endif










c$$$*
c$$$*     Dumping factor for FONLL
c$$$*
c$$$      do i=4,6
c$$$         if(q2.gt.m2th(i))then
c$$$            damp(i) = ( 1d0 - m2th(i) / Q2 )**2d0
c$$$         else
c$$$            damp(i) = 0d0
c$$$         endif
c$$$      enddo
c$$$*
c$$$*     Find "ixi" such that "xigrid(ixi)" < "xi" < "xigrid(ixi+1)"
c$$$*
c$$$      do j=4,6
c$$$         ixi(j) = 0
c$$$         xi(j)  = Q2 / m2th(j)
c$$$         if(xi(j).le.ximin)then
c$$$            ixi(j) = 1
c$$$         elseif(xi(j).ge.ximax)then
c$$$            ixi(j) = nxi
c$$$         else
c$$$            diff(1) = xi(j) - xigrid(1)
c$$$            do i=2,nxi
c$$$               diff(i) = xi(j) - xigrid(i)
c$$$               sgn = diff(i-1) * diff(i)
c$$$               if(sgn.lt.0.d0)then
c$$$                  ixi(j) = i - 1
c$$$               endif
c$$$            enddo
c$$$         endif
c$$$*
c$$$*     Coefficients of the linear interpolation on the xi grid
c$$$*
c$$$         c0(j) = dlog(xigrid(ixi(j)+1)/xi(j))
c$$$     1         / dlog(xigrid(ixi(j)+1)/xigrid(ixi(j)))
c$$$         c1(j) = dlog(xi(j)/xigrid(ixi(j)))
c$$$     1         / dlog(xigrid(ixi(j)+1)/xigrid(ixi(j)))
c$$$      enddo
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$                  do pt=0,ipt
c$$$                     do k=1,3
c$$$                        do i=1,nf
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,nf,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,nf,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC3zm(jgrid,nf,k,pt,0,alpha)
c$$$                        enddo
c$$$                        if(nf.lt.6)then
c$$$                           do i=nf+1,6
c$$$                              SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                           enddo
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         elseif(MassScheme(1:4).eq."FFNS")then
c$$$            do jgrid=1,ngrid
c$$$               do alpha=0,nin(jgrid)
c$$$                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
c$$$*
c$$$*     Light coefficient functions
c$$$*
c$$$                  do i=1,Nf_FF
c$$$                     do pt=0,ipt
c$$$                        do k=1,3
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                        enddo
c$$$*
c$$$*     Add the massive part at NNLO
c$$$*
c$$$                        if(pt.ge.2)then
c$$$                           if(Nf_FF.lt.6)then
c$$$                              do j=Nf_FF+1,6
c$$$                                 if(W2.ge.4d0*m2th(j))then
c$$$                                    SC2(jgrid,i,3,pt,0,alpha) = 
c$$$     1                              SC2(jgrid,i,3,pt,0,alpha)
c$$$     2                     + c0(j) * SC2mNC(jgrid,ixi(j),3,pt,0,alpha)
c$$$     3                     + c1(j) * SC2mNC(jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$                                    SCL(jgrid,i,3,pt,0,alpha) = 
c$$$     1                              SCL(jgrid,i,3,pt,0,alpha)
c$$$     2                     + c0(j) * SCLmNC(jgrid,ixi(j),3,pt,0,alpha)
c$$$     3                     + c1(j) * SCLmNC(jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$                                 endif
c$$$                              enddo
c$$$                           endif
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$*
c$$$*     Heavy coefficient functions
c$$$*
c$$$                  if(Nf_FF.lt.6)then
c$$$                     do i=Nf_FF+1,6
c$$$                        do pt=0,ipt
c$$$                           do k=1,2
c$$$                              if(W2.ge.4d0*m2th(i))then
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                       c0(i) * SC2mNC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                     + c1(i) * SC2mNC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) =
c$$$     1                       c0(i) * SCLmNC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                     + c1(i) * SCLmNC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                       c0(i) * SC3mNC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                     + c1(i) * SC3mNC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                             else
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              endif
c$$$                           enddo
c$$$                           SC2(jgrid,i,3,pt,0,alpha) = 0d0
c$$$                           SCL(jgrid,i,3,pt,0,alpha) = 0d0
c$$$                           SC3(jgrid,i,3,pt,0,alpha) = 0d0
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$         elseif(MassScheme(1:4).eq."FFN0")then
c$$$            do jgrid=1,ngrid
c$$$               do alpha=0,nin(jgrid)
c$$$                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
c$$$*
c$$$*     Light coefficient functions
c$$$*
c$$$                  do i=1,Nf_FF
c$$$                     do pt=0,ipt
c$$$                        do k=1,3
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) =
c$$$     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                        enddo
c$$$*
c$$$*     Add the massive part at NNLO
c$$$*
c$$$                        if(pt.ge.2)then
c$$$                           if(Nf_FF.lt.6)then
c$$$                              do j=Nf_FF+1,6
c$$$                                 if(W2.ge.4d0*m2th(j))then
c$$$                                    SC2(jgrid,i,3,pt,0,alpha) = 
c$$$     1                              SC2(jgrid,i,3,pt,0,alpha)
c$$$     2                    + c0(j) * SC2m0NC(jgrid,ixi(j),3,pt,0,alpha)
c$$$     3                    + c1(j) * SC2m0NC(jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$                                    SCL(jgrid,i,3,pt,0,alpha) = 
c$$$     1                              SCL(jgrid,i,3,pt,0,alpha)
c$$$     2                    + c0(j) * SCLm0NC(jgrid,ixi(j),3,pt,0,alpha)
c$$$     3                    + c1(j) * SCLm0NC(jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$                                 endif
c$$$                              enddo
c$$$                           endif
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$*
c$$$*     Heavy coefficient functions
c$$$*
c$$$                  if(Nf_FF.lt.6)then
c$$$                     do i=Nf_FF+1,6
c$$$                        do pt=0,ipt
c$$$                           do k=1,2
c$$$                              if(W2.ge.4d0*m2th(i))then
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                      c0(i) * SC2m0NC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                    + c1(i) * SC2m0NC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                      c0(i) * SCLm0NC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                    + c1(i) * SCLm0NC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                      c0(i) * SC3m0NC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                    + c1(i) * SC3m0NC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                              else
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              endif
c$$$                           enddo
c$$$                           SC2(jgrid,i,3,pt,0,alpha) = 0d0
c$$$                           SCL(jgrid,i,3,pt,0,alpha) = 0d0
c$$$                           SC3(jgrid,i,3,pt,0,alpha) = 0d0
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$         elseif(MassScheme(1:5).eq."FONLL")then
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$            do jgrid=1,ngrid
c$$$               do alpha=0,nin(jgrid)
c$$$                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
c$$$*
c$$$*     Light coefficient functions
c$$$*
c$$$                  do i=1,Nf_FF
c$$$                     do pt=0,ipt
c$$$                        do k=1,3
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,nf,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,nf,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC3zm(jgrid,nf,k,pt,0,alpha)
c$$$                        enddo
c$$$c$$$*
c$$$c$$$*     Add the massive part at NNLO
c$$$c$$$*
c$$$c$$$                        if(pt.ge.2)then
c$$$c$$$                           if(Nf_FF.lt.6)then
c$$$c$$$                              do j=Nf_FF+1,6
c$$$c$$$                                 if(W2.ge.4d0*m2th(j))then
c$$$c$$$                                    SC2(jgrid,i,3,pt,0,alpha) = 
c$$$c$$$     1                              SC2(jgrid,i,3,pt,0,alpha)
c$$$c$$$     2                     + c0(j) * SC2mNC(jgrid,ixi(j),3,pt,0,alpha)
c$$$c$$$     3                     + c1(j) * SC2mNC(jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$c$$$                                    SCL(jgrid,i,3,pt,0,alpha) = 
c$$$c$$$     1                              SCL(jgrid,i,3,pt,0,alpha)
c$$$c$$$     2                     + c0(j) * SCLmNC(jgrid,ixi(j),3,pt,0,alpha)
c$$$c$$$     3                     + c1(j) * SCLmNC(jgrid,ixi(j)+1,3,pt,0,alpha)
c$$$c$$$                                 endif
c$$$c$$$                              enddo
c$$$c$$$                           endif
c$$$c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$*
c$$$*     Heavy coefficient functions
c$$$*
c$$$                  if(Nf_FF.lt.6)then
c$$$                     do i=Nf_FF+1,6
c$$$                        do pt=0,ipt
c$$$                           do k=1,3
c$$$                              SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$*
c$$$*     Zero Mass Part
c$$$*
c$$$                              if(Q2.ge.m2th(i))then
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                           SC2(jgrid,i,k,pt,0,alpha)
c$$$     2                         + damp(i) * SC2zm(jgrid,nf,k,pt,0,alpha)
c$$$*
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                           SCL(jgrid,i,k,pt,0,alpha)
c$$$     2                         + damp(i) * SCLzm(jgrid,nf,k,pt,0,alpha)
c$$$*
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                           SC3(jgrid,i,k,pt,0,alpha)
c$$$     2                         + damp(i) * SC3zm(jgrid,nf,k,pt,0,alpha)
c$$$                              endif
c$$$                           enddo
c$$$*
c$$$*     Massive Parts
c$$$*
c$$$                           do k=1,2
c$$$                              if(W2.ge.4d0*m2th(i))then
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) =
c$$$     1                           SC2(jgrid,i,k,pt,0,alpha)
c$$$     2               + c0(i) * ( SC2mNC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     3               - damp(i) * SC2m0NC(jgrid,ixi(i),k,pt,0,alpha) )
c$$$     4               + c1(i) * ( SC2mNC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$     5               - damp(i) * SC2m0NC(jgrid,ixi(i)+1,k,pt,0,alpha) )
c$$$
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) =
c$$$     1                           SCL(jgrid,i,k,pt,0,alpha)
c$$$     2               + c0(i) * ( SCLmNC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     3               - damp(i) * SCLm0NC(jgrid,ixi(i),k,pt,0,alpha) )
c$$$     4               + c1(i) * ( SCLmNC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$     5               - damp(i) * SCLm0NC(jgrid,ixi(i)+1,k,pt,0,alpha) )
c$$$
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) =
c$$$     1                           SC3(jgrid,i,k,pt,0,alpha)
c$$$     2               + c0(i) * ( SC3mNC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     3               - damp(i) * SC3m0NC(jgrid,ixi(i),k,pt,0,alpha) )
c$$$     4               + c1(i) * ( SC3mNC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$     5               - damp(i) * SC3m0NC(jgrid,ixi(i)+1,k,pt,0,alpha) )
c$$$                              endif
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$         endif
c$$$*
c$$$*     Compute needed couplings
c$$$*
c$$$         call ComputeChargesDIS(Q2,bq,dq)
c$$$*
c$$$         do jgrid=1,ngrid
c$$$            do alpha=0,nin(jgrid)
c$$$*
c$$$*     Rearrange PDFs in case the target is not a proton.
c$$$*     (This assumes isospin symmetry)
c$$$*
c$$$               if(TargetDIS(1:7).eq."proton")then
c$$$                  frac = 1d0
c$$$               elseif(TargetDIS(1:7).eq."neutron")then
c$$$                  frac = 0d0
c$$$               elseif(TargetDIS.eq."isoscalar")then
c$$$                  frac = 0.5d0
c$$$               elseif(TargetDIS(1:4).eq."iron")then
c$$$                  frac = 0.47166350921037d0 !23.403d0 / 49.618d0
c$$$               endif
c$$$*     Reweight up and down by frac
c$$$               do beta=0,nin(jgrid)
c$$$                  fdw = fph(jgrid,1,beta)
c$$$                  fdb = fph(jgrid,-1,beta)
c$$$                  fup = fph(jgrid,2,beta)
c$$$                  fub = fph(jgrid,-2,beta)
c$$$*
c$$$                  fph(jgrid,1,beta)  = frac * fdw + ( 1d0 - frac ) * fup
c$$$                  fph(jgrid,-1,beta) = frac * fdb + ( 1d0 - frac ) * fub
c$$$                  fph(jgrid,2,beta)  = frac * fup + ( 1d0 - frac ) * fdw
c$$$                  fph(jgrid,-2,beta) = frac * fub + ( 1d0 - frac ) * fdb
c$$$               enddo
c$$$*
c$$$*     F2
c$$$*
c$$$               do i=1,7
c$$$                  F2(i,jgrid,alpha) = 0d0
c$$$               enddo
c$$$               do beta=0,nin(jgrid)-alpha
c$$$                  singlet = 0d0
c$$$                  do i=1,nf
c$$$                     singlet = singlet 
c$$$     1                       + fph(jgrid,i,alpha+beta)
c$$$     2                       + fph(jgrid,-i,alpha+beta)
c$$$                  enddo
c$$$*     F2 flavour by flavour
c$$$                  do i=1,6
c$$$                     do pt=0,ipt
c$$$                        F2t = bq(i) 
c$$$     1                      * ( SC2(jgrid,i,1,pt,0,beta) ! Gluon
c$$$     2                      * fph(jgrid,0,alpha+beta)
c$$$     3                      + SC2(jgrid,i,2,pt,0,beta)   ! Singlet
c$$$     4                      * singlet
c$$$     5                      + SC2(jgrid,i,3,pt,0,beta)   ! Non-singlet
c$$$     6                      * ( fph(jgrid,i,alpha+beta) 
c$$$     7                      + fph(jgrid,-i,alpha+beta) ) )
c$$$*
c$$$                        F2(i,jgrid,alpha) = F2(i,jgrid,alpha)
c$$$     1                                    + as**pt * F2t
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$*     F2 total (combine according to the massive scheme)
c$$$               do i=1,6
c$$$                  F2(7,jgrid,alpha) = F2(7,jgrid,alpha)
c$$$     1                              + F2(i,jgrid,alpha)
c$$$               enddo
c$$$*
c$$$*     FL
c$$$*
c$$$               do i=1,7
c$$$                  FL(i,jgrid,alpha) = 0d0
c$$$               enddo
c$$$               do beta=0,nin(jgrid)-alpha
c$$$                  singlet = 0d0
c$$$                  do i=1,nf
c$$$                     singlet = singlet 
c$$$     1                       + fph(jgrid,i,alpha+beta)
c$$$     2                       + fph(jgrid,-i,alpha+beta)
c$$$                  enddo
c$$$*     FL flavour by flavour
c$$$                  do i=1,6
c$$$                     do pt=0,ipt
c$$$                        FLt = bq(i) 
c$$$     1                      * ( SCL(jgrid,i,1,pt,0,beta) ! Gluon
c$$$     2                      * fph(jgrid,0,alpha+beta)
c$$$     3                      + SCL(jgrid,i,2,pt,0,beta)   ! Singlet
c$$$     4                      * singlet
c$$$     5                      + SCL(jgrid,i,3,pt,0,beta)   ! Non-singlet
c$$$     6                      * ( fph(jgrid,i,alpha+beta) 
c$$$     7                      + fph(jgrid,-i,alpha+beta) ) )
c$$$*
c$$$                        FL(i,jgrid,alpha) = FL(i,jgrid,alpha)
c$$$     1                                    + as**pt * FLt
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$*     FL total (combine according to the massive scheme)
c$$$               do i=1,6
c$$$                  FL(7,jgrid,alpha) = FL(7,jgrid,alpha)
c$$$     1                              + FL(i,jgrid,alpha)
c$$$               enddo
c$$$*
c$$$*     F3
c$$$*
c$$$               do i=1,7
c$$$                  F3(i,jgrid,alpha) = 0d0
c$$$               enddo
c$$$               do beta=0,nin(jgrid)-alpha
c$$$*     F3 flavour by flavour
c$$$                  do i=1,nf
c$$$                     do pt=0,ipt
c$$$                        F3t = dq(i) 
c$$$     1                      * SC3(jgrid,i,3,pt,0,beta)   ! Non-singlet
c$$$     2                      * ( fph(jgrid,i,alpha+beta) 
c$$$     3                      - fph(jgrid,-i,alpha+beta) )
c$$$*
c$$$                        F3(i,jgrid,alpha) = F3(i,jgrid,alpha)
c$$$     1                                    + as**pt * F3t
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$*     F3 total (it is always in the ZM scheme)
c$$$               do i=1,6
c$$$                  F3(7,jgrid,alpha) = F3(7,jgrid,alpha)
c$$$     1                              + F3(i,jgrid,alpha)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      elseif(ProcessDIS.eq."CC")then
c$$$*
c$$$*     Construct the coefficient functions according to the scheme chosen
c$$$*
c$$$         if(MassScheme.eq."ZM-VFNS")then
c$$$            do jgrid=1,ngrid
c$$$               do alpha=0,nin(jgrid)
c$$$                  do pt=0,ipt
c$$$                     do k=1,3
c$$$                        do i=3,nf
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,nf,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,nf,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC3zm(jgrid,nf,k,pt,0,alpha)
c$$$                        enddo
c$$$                        if(nf.lt.6)then
c$$$                           do i=nf+1,6
c$$$                              SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                           enddo
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         elseif(MassScheme(1:4).eq."FFNS")then
c$$$            do jgrid=1,ngrid
c$$$               do alpha=0,nin(jgrid)
c$$$                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
c$$$*
c$$$*     Light coefficient functions
c$$$*
c$$$                  do i=3,Nf_FF
c$$$                     do pt=0,ipt
c$$$                        do k=1,3
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                        enddo
c$$$
c$$$                     enddo
c$$$                  enddo
c$$$*
c$$$*     Heavy coefficient functions
c$$$*
c$$$                  if(Nf_FF.lt.6)then
c$$$                     do i=Nf_FF+1,6
c$$$                        do pt=0,ipt
c$$$                           do k=1,3
c$$$                              if(W2.ge.m2th(i))then
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                       c0(i) * SC2mCC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                     + c1(i) * SC2mCC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) =
c$$$     1                       c0(i) * SCLmCC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                     + c1(i) * SCLmCC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                       c0(i) * SC3mCC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                     + c1(i) * SC3mCC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                              else
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              endif
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$         elseif(MassScheme(1:4).eq."FFN0")then
c$$$            do jgrid=1,ngrid
c$$$               do alpha=0,nin(jgrid)
c$$$                  W2 = Q2 * ( 1d0 - xg(jgrid,alpha) ) / xg(jgrid,alpha)
c$$$*
c$$$*     Light coefficient functions
c$$$*
c$$$                  do i=1,Nf_FF
c$$$                     do pt=0,ipt
c$$$                        do k=1,3
c$$$                           SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SC2zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SCL(jgrid,i,k,pt,0,alpha) = 
c$$$     1                          SCLzm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                           SC3(jgrid,i,k,pt,0,alpha) =
c$$$     1                          SC3zm(jgrid,Nf_FF,k,pt,0,alpha)
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c$$$*
c$$$*     Heavy coefficient functions
c$$$*
c$$$                  if(Nf_FF.lt.6)then
c$$$                     do i=Nf_FF+1,6
c$$$                        do pt=0,ipt
c$$$                           do k=1,3
c$$$                              if(W2.ge.m2th(i))then
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 
c$$$     1                      c0(i) * SC2m0CC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                    + c1(i) * SC2m0CC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) =
c$$$     1                      c0(i) * SCLm0CC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                    + c1(i) * SCLm0CC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 
c$$$     1                      c0(i) * SC3m0CC(jgrid,ixi(i),k,pt,0,alpha)
c$$$     2                    + c1(i) * SC3m0CC(jgrid,ixi(i)+1,k,pt,0,alpha)
c$$$                              else
c$$$                                 SC2(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SCL(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                                 SC3(jgrid,i,k,pt,0,alpha) = 0d0
c$$$                              endif
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$         endif
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$*
c$$$         do jgrid=1,ngrid
c$$$            do alpha=0,nin(jgrid)
c$$$*
c$$$*     Rearrange PDFs in case the target is not a proton.
c$$$*     (This assumes isospin symmetry)
c$$$*
c$$$               if(TargetDIS(1:7).eq."proton")then
c$$$                  frac = 1d0
c$$$               elseif(TargetDIS(1:7).eq."neutron")then
c$$$                  frac = 0d0
c$$$               elseif(TargetDIS.eq."isoscalar")then
c$$$                  frac = 0.5d0
c$$$               elseif(TargetDIS(1:4).eq."iron")then
c$$$                  frac = 0.47166350921037d0 !23.403d0 / 49.618d0
c$$$               endif
c$$$*     Reweight up and down by frac
c$$$               do beta=0,nin(jgrid)
c$$$                  fdw = fph(jgrid,1,beta)
c$$$                  fdb = fph(jgrid,-1,beta)
c$$$                  fup = fph(jgrid,2,beta)
c$$$                  fub = fph(jgrid,-2,beta)
c$$$*
c$$$                  fph(jgrid,1,beta)  = frac * fdw + ( 1d0 - frac ) * fup
c$$$                  fph(jgrid,-1,beta) = frac * fdb + ( 1d0 - frac ) * fub
c$$$                  fph(jgrid,2,beta)  = frac * fup + ( 1d0 - frac ) * fdw
c$$$                  fph(jgrid,-2,beta) = frac * fub + ( 1d0 - frac ) * fdb
c$$$               enddo
c$$$*
c$$$               if(ProjectileDIS(1:8).eq."neutrino".or.
c$$$     1            ProjectileDIS(1:8).eq."positron")then
c$$$                  ipr = 1
c$$$               elseif(ProjectileDIS.eq."antineutrino".or.
c$$$     1            ProjectileDIS(1:8).eq."electron")then
c$$$                  ipr = - 1
c$$$               endif
c$$$*
c$$$*     F2
c$$$*
c$$$               do i=1,7
c$$$                  F2(i,jgrid,alpha) = 0d0
c$$$               enddo
c$$$               do beta=0,nin(jgrid)-alpha
c$$$                  singlet = 0d0
c$$$                  do i=1,nf
c$$$                     singlet = singlet 
c$$$     1                       + fph(jgrid,i,alpha+beta)
c$$$     2                       + fph(jgrid,-i,alpha+beta)
c$$$                  enddo
c$$$*     F2light (all in the component 3 of F2)
c$$$                  do pt=0,ipt
c$$$                     F2t = 2d0
c$$$     1                   * ( ( V_ud2 + V_us2 )        ! Gluon
c$$$     2                   * ( SC2(jgrid,3,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SC2(jgrid,3,2,pt,0,beta)   ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SC2(jgrid,3,3,pt,0,beta)   ! Non-singlet
c$$$     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
c$$$     8                   + ( V_ud2 + V_us2 )
c$$$     9                   * fph(jgrid,-2*ipr,alpha+beta) 
c$$$     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
c$$$*
c$$$                     F2(3,jgrid,alpha) = F2(3,jgrid,alpha)
c$$$     1                                 + as**pt * F2t
c$$$*     F2charm
c$$$                     F2t = 2d0
c$$$     1                   * ( ( V_cd2 + V_cs2 )        ! Gluon
c$$$     2                   * ( SC2(jgrid,4,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SC2(jgrid,4,2,pt,0,beta)   ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SC2(jgrid,4,3,pt,0,beta)   ! Non-singlet
c$$$     7                   * ( V_cd2
c$$$     8                   * ( fph(jgrid,1*ipr,alpha+beta)
c$$$     9                   + fph(jgrid,-4*ipr,alpha+beta) )
c$$$     1                   + V_cs2
c$$$     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
c$$$     3                   + fph(jgrid,-4*ipr,alpha+beta) ) ) )
c$$$*             
c$$$                     F2(4,jgrid,alpha) = F2(4,jgrid,alpha)
c$$$     1                                 + as**pt * F2t
c$$$*     F2bottom
c$$$                     F2t = 2d0
c$$$     1                   * ( ( V_ub2 + V_cb2 )        ! Gluon
c$$$     2                   * ( SC2(jgrid,5,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SC2(jgrid,5,2,pt,0,beta)   ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SC2(jgrid,5,3,pt,0,beta)   ! Non-singlet
c$$$     7                   * ( V_ub2
c$$$     8                   * ( fph(jgrid,-2*ipr,alpha+beta) 
c$$$     9                   + fph(jgrid,5*ipr,alpha+beta) )
c$$$     1                   + V_cb2
c$$$     2                   * ( fph(jgrid,-4*ipr,alpha+beta) 
c$$$     3                   + fph(jgrid,5*ipr,alpha+beta) ) ) )
c$$$*                    
c$$$                     F2(5,jgrid,alpha) = F2(5,jgrid,alpha)
c$$$     1                                 + as**pt * F2t
c$$$*     F2top
c$$$                     F2t = 2d0 
c$$$     1                   * ( ( V_td2 + V_ts2 + V_tb2 ) ! Gluon
c$$$     2                   * ( SC2(jgrid,6,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SC2(jgrid,6,2,pt,0,beta)    ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SC2(jgrid,6,3,pt,0,beta)    ! Non-singlet
c$$$     7                   * ( V_td2
c$$$     8                   * ( fph(jgrid,1*ipr,alpha+beta) 
c$$$     9                   + fph(jgrid,-6*ipr,alpha+beta) )
c$$$     1                   + V_ts2
c$$$     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
c$$$     3                   + fph(jgrid,-6*ipr,alpha+beta) )
c$$$     4                   + V_tb2
c$$$     5                   * ( fph(jgrid,5*ipr,alpha+beta) 
c$$$     6                   + fph(jgrid,-6*ipr,alpha+beta) ) ) )
c$$$*                  
c$$$                     F2(6,jgrid,alpha) = F2(6,jgrid,alpha)
c$$$     1                                 + as**pt * F2t
c$$$                  enddo
c$$$               enddo
c$$$*     F2 total
c$$$               do i=3,6
c$$$                  F2(7,jgrid,alpha) = F2(7,jgrid,alpha)
c$$$     1                              + F2(i,jgrid,alpha)
c$$$               enddo
c$$$*
c$$$*     FL
c$$$*
c$$$               do i=1,7
c$$$                  FL(i,jgrid,alpha) = 0d0
c$$$               enddo
c$$$               do beta=0,nin(jgrid)-alpha
c$$$                  singlet = 0d0
c$$$                  do i=1,nf
c$$$                     singlet = singlet 
c$$$     1                       + fph(jgrid,i,alpha+beta)
c$$$     2                       + fph(jgrid,-i,alpha+beta)
c$$$                  enddo
c$$$*     FLlight (all in the component 3 of FL)
c$$$                  do pt=0,ipt
c$$$                     FLt = 2d0
c$$$     1                   * ( ( V_ud2 + V_us2 )        ! Gluon
c$$$     2                   * ( SCL(jgrid,3,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SCL(jgrid,3,2,pt,0,beta)   ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SCL(jgrid,3,3,pt,0,beta)   ! Non-singlet
c$$$     7                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
c$$$     8                   + ( V_ud2 + V_us2 )
c$$$     9                   * fph(jgrid,-2*ipr,alpha+beta) 
c$$$     1                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
c$$$*
c$$$                     FL(3,jgrid,alpha) = FL(3,jgrid,alpha)
c$$$     1                                 + as**pt * FLt
c$$$*     FLcharm
c$$$                     FLt = 2d0
c$$$     1                   * ( ( V_cd2 + V_cs2 )        ! Gluon
c$$$     2                   * ( SCL(jgrid,4,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SCL(jgrid,4,2,pt,0,beta)   ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SCL(jgrid,4,3,pt,0,beta)   ! Non-singlet
c$$$     7                   * ( V_cd2
c$$$     8                   * ( fph(jgrid,1*ipr,alpha+beta) 
c$$$     9                   + fph(jgrid,-4*ipr,alpha+beta) )
c$$$     1                   + V_cs2
c$$$     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
c$$$     3                   + fph(jgrid,-4*ipr,alpha+beta) ) ) )
c$$$*                 
c$$$                     FL(4,jgrid,alpha) = FL(4,jgrid,alpha)
c$$$     1                                 + as**pt * FLt
c$$$*     FLbottom
c$$$                     FLt = 2d0
c$$$     1                   * ( ( V_ub2 + V_cb2 )        ! Gluon
c$$$     2                   * ( SCL(jgrid,5,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SCL(jgrid,5,2,pt,0,beta)   ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SCL(jgrid,5,3,pt,0,beta)   ! Non-singlet
c$$$     7                   * ( V_ub2
c$$$     8                   * ( fph(jgrid,-2*ipr,alpha+beta) 
c$$$     9                   + fph(jgrid,5*ipr,alpha+beta) )
c$$$     1                   + V_cb2
c$$$     2                   * ( fph(jgrid,-4*ipr,alpha+beta) 
c$$$     3                   + fph(jgrid,5*ipr,alpha+beta) ) ) )
c$$$*                   
c$$$                     FL(5,jgrid,alpha) = FL(5,jgrid,alpha)
c$$$     1                                 + as**pt * FLt
c$$$*     FLtop
c$$$                     FLt = 2d0 
c$$$     1                   * ( ( V_td2 + V_ts2 + V_tb2 ) ! Gluon
c$$$     2                   * ( SCL(jgrid,6,1,pt,0,beta)
c$$$     3                   * fph(jgrid,0,alpha+beta)
c$$$     4                   + SCL(jgrid,6,2,pt,0,beta)    ! Singlet
c$$$     5                   * singlet )
c$$$     6                   + SCL(jgrid,6,3,pt,0,beta)    ! Non-singlet
c$$$     7                   * ( V_td2
c$$$     8                   * ( fph(jgrid,1*ipr,alpha+beta) 
c$$$     9                   + fph(jgrid,-6*ipr,alpha+beta) )
c$$$     1                   + V_ts2
c$$$     2                   * ( fph(jgrid,3*ipr,alpha+beta) 
c$$$     3                   + fph(jgrid,-6*ipr,alpha+beta) )
c$$$     4                   + V_tb2
c$$$     5                   * ( fph(jgrid,5*ipr,alpha+beta) 
c$$$     6                   + fph(jgrid,-6*ipr,alpha+beta) ) ) )
c$$$*                 
c$$$                     FL(6,jgrid,alpha) = FL(6,jgrid,alpha)
c$$$     1                                 + as**pt * FLt
c$$$                  enddo
c$$$               enddo
c$$$*     FL total
c$$$               do i=3,6
c$$$                  FL(7,jgrid,alpha) = FL(7,jgrid,alpha)
c$$$     1                              + FL(i,jgrid,alpha)
c$$$               enddo
c$$$*
c$$$*     F3
c$$$*
c$$$               do i=1,7
c$$$                  F3(i,jgrid,alpha) = 0d0
c$$$               enddo
c$$$               do beta=0,nin(jgrid)-alpha
c$$$*     F3light (all in the component 3 of F3)
c$$$                  do pt=0,ipt
c$$$                     F3t = 2d0
c$$$     1                   * ( SC3(jgrid,3,3,pt,0,beta)        ! Non-singlet
c$$$     2                   * ( V_ud2 * fph(jgrid,1*ipr,alpha+beta) 
c$$$     3                   - ( V_ud2 + V_us2 )
c$$$     4                   * fph(jgrid,-2*ipr,alpha+beta) 
c$$$     5                   + V_us2 * fph(jgrid,3*ipr,alpha+beta) ) )
c$$$*
c$$$                     F3(3,jgrid,alpha) = F3(3,jgrid,alpha)
c$$$     1                                 + as**pt * F3t
c$$$*     F3charm
c$$$                     F3t = 2d0 * ( ( V_cd2 + V_cs2 )         ! Gluon
c$$$     1                   * SC3(jgrid,4,1,pt,0,beta)
c$$$     2                   * fph(jgrid,0,alpha+beta)
c$$$     3                   + SC3(jgrid,4,3,pt,0,beta)          ! Non-singlet
c$$$     4                   * ( V_cd2
c$$$     5                   * ( fph(jgrid,1*ipr,alpha+beta) 
c$$$     6                   - fph(jgrid,-4*ipr,alpha+beta) )
c$$$     7                   + V_cs2
c$$$     8                   * ( fph(jgrid,3*ipr,alpha+beta) 
c$$$     9                   - fph(jgrid,-4*ipr,alpha+beta) ) ) )
c$$$*                 
c$$$                     F3(4,jgrid,alpha) = F3(4,jgrid,alpha)
c$$$     1                                 + as**pt * F3t
c$$$*     F3bottom
c$$$                     F3t = 2d0
c$$$     1                   * ( SC3(jgrid,5,3,pt,0,beta)   ! Non-singlet
c$$$     2                   * ( V_ub2
c$$$     3                   * ( - fph(jgrid,-2*ipr,alpha+beta) 
c$$$     4                   + fph(jgrid,5*ipr,alpha+beta) )
c$$$     5                   + V_cb2
c$$$     6                   * ( - fph(jgrid,-4*ipr,alpha+beta) 
c$$$     7                   + fph(jgrid,5*ipr,alpha+beta) ) ) )
c$$$*                 
c$$$                     F3(5,jgrid,alpha) = F3(5,jgrid,alpha)
c$$$     1                                 + as**pt * F3t
c$$$*     F3top
c$$$                     F3t = 2d0 
c$$$     1                   * ( SC3(jgrid,6,3,pt,0,beta)   ! Non-singlet
c$$$     2                   * ( V_td2
c$$$     3                   * ( fph(jgrid,1*ipr,alpha+beta) 
c$$$     4                   - fph(jgrid,-6*ipr,alpha+beta) )
c$$$     5                   + V_ts2
c$$$     6                   * ( fph(jgrid,3*ipr,alpha+beta) 
c$$$     7                   - fph(jgrid,-6*ipr,alpha+beta) )
c$$$     8                   + V_tb2
c$$$     9                   * ( fph(jgrid,5*ipr,alpha+beta) 
c$$$     1                   - fph(jgrid,-6*ipr,alpha+beta) ) ) )
c$$$*                 
c$$$                     F3(6,jgrid,alpha) = F3(6,jgrid,alpha)
c$$$     1                                 + as**pt * F3t
c$$$                  enddo
c$$$               enddo
c$$$*     F3 total
c$$$               do i=3,6
c$$$                  F3(7,jgrid,alpha) = F3(7,jgrid,alpha)
c$$$     1                              + F3(i,jgrid,alpha)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      endif

      call cpu_time(t2)
*
      write(6,"(a,a,f7.3,a)") " Computation of the DIS operators",
     1                        " completed in",t2-t1," s"
      write(6,*) " "
*
      return
      end 
