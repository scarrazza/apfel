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
      integer ipr
      double precision Q2
      double precision as(0:2),a_QCD
      double precision bq(6),dq(6)
      double precision frac,fr3
      double precision C2g,C2ps,C2nsp,C2nsm
      double precision CLg,CLps,CLnsp,CLnsm
      double precision C3g,C3nsp,C3nsm
      double precision Kl,Kc,Kb,Kt
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
               ipr = 1
            elseif(ProjectileDIS.eq."antineutrino".or.
     1             ProjectileDIS(1:8).eq."electron")then
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
*     F2
*
*     Put together coefficient functions
*
                     C2g   = 0d0
                     C2ps  = 0d0
                     C2nsp = 0d0
                     C2nsm = 0d0
                     do pt=0,ipt
                        C2g   = C2g
     1                        + as(pt) * SC2zm(jgrid,nf,1,pt,alpha,beta)
                        C2ps  = C2ps
     1                        + as(pt) * SC2zm(jgrid,nf,2,pt,alpha,beta)
                        C2nsp = C2nsp
     1                        + as(pt) * SC2zm(jgrid,nf,3,pt,alpha,beta)
                        C2nsm = C2nsm
     1                        + as(pt) * SC2zm(jgrid,nf,4,pt,alpha,beta)
                     enddo
*     
*     Light Component
*     
*     Singlet
                     OpF2(jgrid,3,1,alpha,beta)  = 
     1                    2d0 * Kl * ( C2ps + C2nsp / 6d0 )
*     Gluon
                     OpF2(jgrid,3,2,alpha,beta)  = 2d0 * Kl * C2g
*     V3
                     OpF2(jgrid,3,4,alpha,beta)  = - ipr * fr3
     1                    * ( 2d0 * V_ud2 + V_us2 ) * C2nsm / 2d0
*     V8
                     OpF2(jgrid,3,5,alpha,beta)  = - ipr
     1                    * V_us2 * C2nsm / 2d0
*     T3
                     OpF2(jgrid,3,9,alpha,beta)  = fr3
     1                    * V_us2 * C2nsp / 2d0
*     T8
                     OpF2(jgrid,3,10,alpha,beta) = 
     1                    ( 2d0 * V_ud2 - V_us2 ) * C2nsp / 6d0
*     T15
                     OpF2(jgrid,3,11,alpha,beta) = Kl * C2nsp / 6d0
*     T24
                     OpF2(jgrid,3,12,alpha,beta) = Kl * C2nsp / 10d0
*     T35
                     OpF2(jgrid,3,13,alpha,beta) = Kl * C2nsp / 15d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Singlet
                        OpF2(jgrid,4,1,alpha,beta)  = 
     1                    2d0 * Kc * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,4,2,alpha,beta)  = 2d0 * Kc * C2g
*     V3
                        OpF2(jgrid,4,4,alpha,beta)  = 
     1                       - ipr * fr3 * V_cd2 * C2nsm / 2d0
*     V8
                        OpF2(jgrid,4,5,alpha,beta)  =
     1                       ipr * ( V_cd2 - 2d0 * V_cs2 ) * C2nsm / 6d0
*     V15
                        OpF2(jgrid,4,6,alpha,beta)  =
     1                       ipr * Kc * C2nsm / 3d0
*     T3
                        OpF2(jgrid,4,9,alpha,beta)  = 
     1                       - fr3 * V_cd2 * C2nsp / 2d0
*     T8
                        OpF2(jgrid,4,10,alpha,beta) = 
     1                       ( V_cd2 - 2d0 * V_cs2 ) * C2nsp / 6d0
*     T15
                        OpF2(jgrid,4,11,alpha,beta) = - Kc * C2nsp / 6d0
*     T24
                        OpF2(jgrid,4,12,alpha,beta) = Kc * C2nsp / 10d0
*     T35
                        OpF2(jgrid,4,13,alpha,beta) = Kc * C2nsp / 15d0
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Singlet
                        OpF2(jgrid,5,1,alpha,beta)  = 
     1                    2d0 * Kb * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,5,2,alpha,beta)  = 2d0 * Kb * C2g
*     V3
                        OpF2(jgrid,5,4,alpha,beta)  =
     1                       - ipr * fr3 * V_ub2 * C2nsm / 2d0
*     V8
                        OpF2(jgrid,5,5,alpha,beta)  = 
     1                       - ipr * V_ub2 * C2nsm / 6d0
*     V15
                        OpF2(jgrid,5,6,alpha,beta)  = - ipr
     1                       * ( V_ub2 - 3d0 * V_cb2 ) * C2nsm / 12d0
*     V24
                        OpF2(jgrid,5,7,alpha,beta)  = 
     1                       - ipr * Kb * C2nsm / 4d0
*     T3
                        OpF2(jgrid,5,9,alpha,beta)  =
     1                       fr3 * V_ub2 * C2nsp / 2d0 
*     T8
                        OpF2(jgrid,5,10,alpha,beta) = 
     1                       V_ub2 * C2nsp / 6d0
*     T15
                        OpF2(jgrid,5,11,alpha,beta) = 
     1                       ( V_ub2 - 3d0 * V_cb2 ) * C2nsp / 12d0
*     T24
                        OpF2(jgrid,5,12,alpha,beta) = 
     1                       - 3d0 * Kb * C2nsp / 20d0
*     T35
                        OpF2(jgrid,5,13,alpha,beta) = Kb * C2nsp / 15d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Singlet
                        OpF2(jgrid,6,1,alpha,beta)  = 
     1                    2d0 * Kt * ( C2ps + C2nsp / 6d0 )
*     Gluon
                        OpF2(jgrid,6,2,alpha,beta)  = 2d0 * Kt * C2g
*     V3
                        OpF2(jgrid,6,4,alpha,beta)  =
     1                       - ipr * fr3 * V_td2 * C2nsm / 2d0
*     V8
                        OpF2(jgrid,6,5,alpha,beta)  = 
     1                       ipr * ( V_td2 - 2d0 * V_ts2 ) * C2nsm / 6d0
*     V15
                        OpF2(jgrid,6,6,alpha,beta)  =
     1                       ipr * ( V_td2 + V_ts2 ) * C2nsm / 12d0
*     V24
                        OpF2(jgrid,6,7,alpha,beta)  = 
     1                       ipr * ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                       * C2nsm / 20d0
*     V35
                        OpF2(jgrid,6,8,alpha,beta)  =
     1                       ipr * Kt * C2nsm / 5d0
*     T3
                        OpF2(jgrid,6,9,alpha,beta)  =
     1                       - fr3 * V_td2 * C2nsp / 2d0 
*     T8
                        OpF2(jgrid,6,10,alpha,beta) =
     1                       ( V_td2 - 2d0 * V_ts2 ) * C2nsp / 6d0
*     T15
                        OpF2(jgrid,6,11,alpha,beta) = 
     1                       ( V_td2 + V_ts2 ) * C2nsp / 12d0
*     T24
                        OpF2(jgrid,6,12,alpha,beta) = 
     1                       ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                       * C2nsp / 20d0
*     T35
                        OpF2(jgrid,6,13,alpha,beta) =
     1                       - 2d0 * Kt * C2nsp / 15d0
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
                     CLnsm = 0d0
                     do pt=0,ipt
                        CLg   = CLg
     1                        + as(pt) * SCLzm(jgrid,nf,1,pt,alpha,beta)
                        CLps  = CLps
     1                        + as(pt) * SCLzm(jgrid,nf,2,pt,alpha,beta)
                        CLnsp = CLnsp
     1                        + as(pt) * SCLzm(jgrid,nf,3,pt,alpha,beta)
                        CLnsm = CLnsm
     1                        + as(pt) * SCLzm(jgrid,nf,4,pt,alpha,beta)
                     enddo
*     
*     Light Component
*     
*     Singlet
                     OpFL(jgrid,3,1,alpha,beta)  = 
     1                    2d0 * Kl * ( CLps + CLnsp / 6d0 )
*     Gluon
                     OpFL(jgrid,3,2,alpha,beta)  = 2d0 * Kl * CLg
*     V3
                     OpFL(jgrid,3,4,alpha,beta)  = - ipr * fr3
     1                    * ( 2d0 * V_ud2 + V_us2 ) * CLnsm / 2d0
*     V8
                     OpFL(jgrid,3,5,alpha,beta)  = - ipr
     1                    * V_us2 * CLnsm / 2d0
*     T3
                     OpFL(jgrid,3,9,alpha,beta)  = fr3
     1                    * V_us2 * CLnsp / 2d0
*     T8
                     OpFL(jgrid,3,10,alpha,beta) = 
     1                    ( 2d0 * V_ud2 - V_us2 ) * CLnsp / 6d0
*     T15
                     OpFL(jgrid,3,11,alpha,beta) = Kl * CLnsp / 6d0
*     T24
                     OpFL(jgrid,3,12,alpha,beta) = Kl * CLnsp / 10d0
*     T35
                     OpFL(jgrid,3,13,alpha,beta) = Kl * CLnsp / 15d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Singlet
                        OpFL(jgrid,4,1,alpha,beta)  = 
     1                    2d0 * Kc * ( CLps + CLnsp / 6d0 )
*     Gluon
                        OpFL(jgrid,4,2,alpha,beta)  = 2d0 * Kc * CLg
*     V3
                        OpFL(jgrid,4,4,alpha,beta)  = 
     1                       - ipr * fr3 * V_cd2 * CLnsm / 2d0
*     V8
                        OpFL(jgrid,4,5,alpha,beta)  =
     1                       ipr * ( V_cd2 - 2d0 * V_cs2 ) * CLnsm / 6d0
*     V15
                        OpFL(jgrid,4,6,alpha,beta)  =
     1                       ipr * Kc * CLnsm / 3d0
*     T3
                        OpFL(jgrid,4,9,alpha,beta)  = 
     1                       - fr3 * V_cd2 * CLnsp / 2d0
*     T8
                        OpFL(jgrid,4,10,alpha,beta) = 
     1                       ( V_cd2 - 2d0 * V_cs2 ) * CLnsp / 6d0
*     T15
                        OpFL(jgrid,4,11,alpha,beta) = - Kc * CLnsp / 6d0
*     T24
                        OpFL(jgrid,4,12,alpha,beta) = Kc * CLnsp / 10d0
*     T35
                        OpFL(jgrid,4,13,alpha,beta) = Kc * CLnsp / 15d0
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Singlet
                        OpFL(jgrid,5,1,alpha,beta)  = 
     1                    2d0 * Kb * ( CLps + CLnsp / 6d0 )
*     Gluon
                        OpFL(jgrid,5,2,alpha,beta)  = 2d0 * Kb * CLg
*     V3
                        OpFL(jgrid,5,4,alpha,beta)  =
     1                       - ipr * fr3 * V_ub2 * CLnsm / 2d0
*     V8
                        OpFL(jgrid,5,5,alpha,beta)  = 
     1                       - ipr * V_ub2 * CLnsm / 6d0
*     V15
                        OpFL(jgrid,5,6,alpha,beta)  = - ipr
     1                       * ( V_ub2 - 3d0 * V_cb2 ) * CLnsm / 12d0
*     V24
                        OpFL(jgrid,5,7,alpha,beta)  = 
     1                       - ipr * Kb * CLnsm / 4d0
*     T3
                        OpFL(jgrid,5,9,alpha,beta)  =
     1                       fr3 * V_ub2 * CLnsp / 2d0 
*     T8
                        OpFL(jgrid,5,10,alpha,beta) = 
     1                       V_ub2 * CLnsp / 6d0
*     T15
                        OpFL(jgrid,5,11,alpha,beta) = 
     1                       ( V_ub2 - 3d0 * V_cb2 ) * CLnsp / 12d0
*     T24
                        OpFL(jgrid,5,12,alpha,beta) = 
     1                       - 3d0 * Kb * CLnsp / 20d0
*     T35
                        OpFL(jgrid,5,13,alpha,beta) = Kb * CLnsp / 15d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Singlet
                        OpFL(jgrid,6,1,alpha,beta)  = 
     1                    2d0 * Kt * ( CLps + CLnsp / 6d0 )
*     Gluon
                        OpFL(jgrid,6,2,alpha,beta)  = 2d0 * Kt * CLg
*     V3
                        OpFL(jgrid,6,4,alpha,beta)  =
     1                       - ipr * fr3 * V_td2 * CLnsm / 2d0
*     V8
                        OpFL(jgrid,6,5,alpha,beta)  = 
     1                       ipr * ( V_td2 - 2d0 * V_ts2 ) * CLnsm / 6d0
*     V15
                        OpFL(jgrid,6,6,alpha,beta)  =
     1                       ipr * ( V_td2 + V_ts2 ) * CLnsm / 12d0
*     V24
                        OpFL(jgrid,6,7,alpha,beta)  = 
     1                       ipr * ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                       * CLnsm / 20d0
*     V35
                        OpFL(jgrid,6,8,alpha,beta)  =
     1                       ipr * Kt * CLnsm / 5d0
*     T3
                        OpFL(jgrid,6,9,alpha,beta)  =
     1                       - fr3 * V_td2 * CLnsp / 2d0 
*     T8
                        OpFL(jgrid,6,10,alpha,beta) =
     1                       ( V_td2 - 2d0 * V_ts2 ) * CLnsp / 6d0
*     T15
                        OpFL(jgrid,6,11,alpha,beta) = 
     1                       ( V_td2 + V_ts2 ) * CLnsp / 12d0
*     T24
                        OpFL(jgrid,6,12,alpha,beta) = 
     1                       ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                       * CLnsp / 20d0
*     T35
                        OpFL(jgrid,6,13,alpha,beta) =
     1                       - 2d0 * Kt * CLnsp / 15d0
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
*     F3
*
*     Put together coefficient functions
*
                     C3g   = 0d0
                     C3nsp = 0d0
                     C3nsm = 0d0
                     do pt=0,ipt
                        C3nsp = C3nsp
     1                        + as(pt) * SC3zm(jgrid,nf,3,pt,alpha,beta)
                        C3nsm = C3nsm
     1                        + as(pt) * SC3zm(jgrid,nf,4,pt,alpha,beta)
                     enddo
*     
*     Light Component
*     
*     Gluon
                     OpF3(jgrid,3,2,alpha,beta)  = 2d0 * Kl * C3g
*     Valence
                     OpF3(jgrid,3,3,alpha,beta)  = Kl * C3nsm / 3d0
*     V3
                     OpF3(jgrid,3,4,alpha,beta)  = fr3
     1                    * V_us2 * C3nsp / 2d0
*     V8
                     OpF3(jgrid,3,5,alpha,beta)  = 
     1                    ( 2d0 * V_ud2 - V_us2 ) * C3nsm / 6d0
*     V15
                     OpF3(jgrid,3,6,alpha,beta)  = Kl * C3nsm / 6d0
*     V24
                     OpF3(jgrid,3,7,alpha,beta)  = Kl * C3nsm / 10d0
*     V35
                     OpF3(jgrid,3,8,alpha,beta)  = Kl * C3nsm / 15d0
*     T3
                     OpF3(jgrid,3,9,alpha,beta)  = - ipr * fr3
     1                    * ( 2d0 * V_ud2 + V_us2 ) * C3nsp / 2d0
*     T8
                     OpF3(jgrid,3,10,alpha,beta) = - ipr
     1                    * V_us2 * C3nsp / 2d0
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Gluon
                        OpF3(jgrid,4,2,alpha,beta)  = 2d0 * Kc * C3g
*     Valence
                        OpF3(jgrid,4,3,alpha,beta)  = Kc * C3nsm / 3d0
*     V3
                        OpF3(jgrid,4,4,alpha,beta)  = 
     1                       - fr3 * V_cd2 * C3nsm / 2d0
*     V8
                        OpF3(jgrid,4,5,alpha,beta)  = 
     1                       ( V_cd2 - 2d0 * V_cs2 ) * C3nsm / 6d0
*     V15
                        OpF3(jgrid,4,6,alpha,beta)  = - Kc * C3nsm / 6d0
*     V24
                        OpF3(jgrid,4,7,alpha,beta)  = Kc * C3nsm / 10d0
*     V35
                        OpF3(jgrid,4,8,alpha,beta)  = Kc * C3nsm / 15d0
*     T3
                        OpF3(jgrid,4,9,alpha,beta)  = 
     1                       - ipr * fr3 * V_cd2 * C3nsp / 2d0
*     T8
                        OpF3(jgrid,4,10,alpha,beta) =
     1                       ipr * ( V_cd2 - 2d0 * V_cs2 ) * C3nsp / 6d0
*     T15
                        OpF3(jgrid,4,11,alpha,beta) =
     1                       ipr * Kc * C3nsp / 3d0
                     endif






















*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Gluon
                        OpF3(jgrid,5,2,alpha,beta)  = 2d0 * Kb * C3g
*     Valence
                        OpF3(jgrid,5,3,alpha,beta)  = Kb * C3nsm / 3d0
*     V3
                        OpF3(jgrid,5,4,alpha,beta)  =
     1                       fr3 * V_ub2 * C3nsm / 2d0 
*     V8
                        OpF3(jgrid,5,5,alpha,beta)  = 
     1                       V_ub2 * C3nsm / 6d0
*     V15
                        OpF3(jgrid,5,6,alpha,beta)  = 
     1                       ( V_ub2 - 3d0 * V_cb2 ) * C3nsm / 12d0
*     V24
                        OpF3(jgrid,5,7,alpha,beta)  = 
     1                       - 3d0 * Kb * C3nsm / 20d0
*     V35
                        OpF3(jgrid,5,8,alpha,beta)  = Kb * C3nsm / 15d0
*     T3
                        OpF3(jgrid,5,9,alpha,beta)  =
     1                       - ipr * fr3 * V_ub2 * C3nsp / 2d0
*     T8
                        OpF3(jgrid,5,10,alpha,beta) = 
     1                       - ipr * V_ub2 * C3nsp / 6d0
*     T15
                        OpF3(jgrid,5,11,alpha,beta) = - ipr
     1                       * ( V_ub2 - 3d0 * V_cb2 ) * C3nsp / 12d0
*     T24
                        OpF3(jgrid,5,12,alpha,beta) = 
     1                       - ipr * Kb * C3nsp / 4d0
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Gluon
                        OpF3(jgrid,6,2,alpha,beta)  = 2d0 * Kt * C3g
*     Valence
                        OpF3(jgrid,6,3,alpha,beta)  = Kt * C3nsm / 3d0
*     V3
                        OpF3(jgrid,6,4,alpha,beta)  =
     1                       - fr3 * V_td2 * C3nsp / 2d0 
*     V8
                        OpF3(jgrid,6,5,alpha,beta)  =
     1                       ( V_td2 - 2d0 * V_ts2 ) * C3nsp / 6d0
*     V15
                        OpF3(jgrid,6,6,alpha,beta)  = 
     1                       ( V_td2 + V_ts2 ) * C3nsp / 12d0
*     V24
                        OpF3(jgrid,6,7,alpha,beta)  = 
     1                       ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                       * C3nsp / 20d0
*     V35
                        OpF3(jgrid,6,8,alpha,beta)  =
     1                       - 2d0 * Kt * C3nsp / 15d0
*     T3
                        OpF3(jgrid,6,9,alpha,beta)  =
     1                       - ipr * fr3 * V_td2 * C3nsp / 2d0
*     T8
                        OpF3(jgrid,6,10,alpha,beta) = 
     1                       ipr * ( V_td2 - 2d0 * V_ts2 ) * C3nsp / 6d0
*     T15
                        OpF3(jgrid,6,11,alpha,beta) =
     1                       ipr * ( V_td2 + V_ts2 ) * C3nsp / 12d0
*     T24
                        OpF3(jgrid,6,12,alpha,beta) = 
     1                       ipr * ( V_td2 + V_ts2 - 4d0 * V_tb2 )
     2                       * C3nsp / 20d0
*     T35
                        OpF3(jgrid,6,13,alpha,beta) =
     1                       ipr * Kt * C3nsp / 5d0
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
                  enddo
               enddo
            enddo
         endif
      endif
*
      call cpu_time(t2)
*
      write(6,"(a,a,f7.3,a)") " Computation of the DIS operators",
     1                        " completed in",t2-t1," s"
      write(6,*) " "
*
      return
      end 
