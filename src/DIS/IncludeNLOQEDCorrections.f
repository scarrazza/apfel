************************************************************************
*
*     IncludeNLOQEDCorrections.f:
*
*     This routine include to the DIS operators the O(alpha) corrections
*     including that proportional to the photon PDF.
*
************************************************************************
      subroutine IncludeNLOQEDCorrections(Q2,muF2,nf,fr3,jxi,c0,c1,damp,
     1                                    bqi,dqi,ipr)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/DISOperators.h"
      include "../commons/m2th.h"
      include "../commons/ProcessDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/Nf_FF.h"
**
*     Internal Variables
*
      integer nf
      integer ipr
      integer jxi(4:6)
      double precision Q2,muF2
      double precision damp(4:6)
      double precision bqi(0:6),dqi(0:6)
      double precision c0(4:6),c1(4:6)
      double precision fr3
**
*     Internal Variables
*
      integer jgrid,ipdf,ihq
      integer alpha,beta,gbound
      integer ik
      double precision W2
      double precision bq(0:6),dq(0:6),e2q9(6),e2u,e2d
      double precision aem,a_QED
      double precision etap,etam
      double precision C2ph(3:6),C2ns(3:6)
      double precision CLph(3:6),CLns(3:6)
      double precision C3ph(3:6),C3ns(3:6)
      double precision C2gm(4:6),C2qm(4:6)
      double precision CLgm(4:6),CLqm(4:6)
      double precision C3gm(4:6),C3qm(4:6)
      data e2q9 / 1d0, 4d0, 1d0, 4d0, 1d0, 4d0 /
      parameter(e2u=4d0/9d0,e2d=1d0/9d0)
*
*     Compute alpha
*     (Remember that a_QCD takes as an argument the factorization scale
*     and converts it internally into the renormalization scale).
*
      aem  = a_QED(muF2)
*
*     Electromagnetic and Neutral current structure functions
*
      if(ProcessDIS.eq."EM".or.ProcessDIS.eq."NC")then
*
*     Rescale EW couplings
*
         do ihq=1,6
            bq(ihq) = e2q9(ihq) * bqi(ihq) / 9d0
            dq(ihq) = e2q9(ihq) * dqi(ihq) / 9d0
         enddo
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
*     Initialize coefficient functions
*     
                  do ihq=3,6
                     C2ph(ihq) = 0d0
                     C2ns(ihq) = 0d0
                     CLph(ihq) = 0d0
                     CLns(ihq) = 0d0
                     C3ns(ihq) = 0d0
                  enddo
*     
*     Construct coefficient functions according to the scheme chosen
*     
                  if(MassScheme.eq."ZM-VFNS")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,nf,1,1,alpha,beta) / TR
                     C2ns(3) = aem * SC2zm(jgrid,nf,3,1,alpha,beta) / CF
                     CLph(3) = aem * SCLzm(jgrid,nf,1,1,alpha,beta) / TR
                     CLns(3) = aem * SCLzm(jgrid,nf,3,1,alpha,beta) / CF
                     C3ns(3) = aem * SC3zm(jgrid,nf,4,1,alpha,beta) / CF
*     Heavy quark coefficient functions
                     if(nf.gt.3)then
                        do ihq=4,nf
                           C2ph(ihq) = C2ph(3)
                           C2ns(ihq) = C2ns(3)
                           CLph(ihq) = CLph(3)
                           CLns(ihq) = CLns(3)
                           C3ns(ihq) = C3ns(3)
                        enddo
                     endif
                  elseif(MassScheme(1:4).eq."FFN0")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     C2ns(3) = aem * SC2zm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     CLph(3) = aem * SCLzm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     CLns(3) = aem * SCLzm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     C3ns(3) = aem * SC3zm(jgrid,Nf_FF,4,1,alpha,beta)
     1                    / CF
*     Heavy quark coefficient functions
                     if(Nf_FF.gt.3)then
                        do ihq=4,Nf_FF
                           C2ph(ihq) = C2ph(3)
                           C2ns(ihq) = C2ns(3)
                           CLph(ihq) = CLph(3)
                           CLns(ihq) = CLns(3)
                           C3ns(ihq) = C3ns(3)
                        enddo
                     endif
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           if(Q2.ge.m2th(ihq))then
                              C2ph(ihq) = C2ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2m0NC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2m0NC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              CLph(ihq) = CLph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLm0NC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLm0NC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                           endif
                        enddo
                     endif
                  elseif(MassScheme(1:4).eq."FFNS")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     C2ns(3) = aem * SC2zm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     CLph(3) = aem * SCLzm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     CLns(3) = aem * SCLzm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     C3ns(3) = aem * SC3zm(jgrid,Nf_FF,4,1,alpha,beta)
     1                    / CF
*     Heavy quark coefficient functions
                     if(Nf_FF.gt.3)then
                        do ihq=4,Nf_FF
                           C2ph(ihq) = C2ph(3)
                           C2ns(ihq) = C2ns(3)
                           CLph(ihq) = CLph(3)
                           CLns(ihq) = CLns(3)
                           C3ns(ihq) = C3ns(3)
                        enddo
                     endif
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           if(W2.ge.4d0*m2th(ihq))then
                              C2ph(ihq) = C2ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2mNC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2mNC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              CLph(ihq) = CLph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLmNC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLmNC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                           endif
                        enddo
                     endif
                  elseif(MassScheme(1:5).eq."FONLL")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,nf,1,1,alpha,beta) / TR
                     C2ns(3) = aem * SC2zm(jgrid,nf,3,1,alpha,beta) / CF
                     CLph(3) = aem * SCLzm(jgrid,nf,1,1,alpha,beta) / TR
                     CLns(3) = aem * SCLzm(jgrid,nf,3,1,alpha,beta) / CF
                     C3ns(3) = aem * SC3zm(jgrid,nf,4,1,alpha,beta) / CF
*     Heavy quark coefficient functions
                     if(nf.gt.3)then
                        do ihq=4,nf
                           C2ph(ihq) = damp(ihq) * C2ph(3)
                           C2ns(ihq) = damp(ihq) * C2ns(3)
                           CLph(ihq) = damp(ihq) * CLph(3)
                           CLns(ihq) = damp(ihq) * CLns(3)
                           C3ns(ihq) = damp(ihq) * C3ns(3)
                        enddo
                     endif
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           if(W2.ge.4d0*m2th(ihq))then
                              C2ph(ihq) = C2ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2mNC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2mNC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              CLph(ihq) = CLph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLmNC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLmNC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                           endif
                           if(Q2.ge.m2th(ihq))then
                                 C2ph(ihq) = C2ph(ihq) + aem
     1                                * ( - damp(ihq) * ( c0(ihq)
     2                                * SC2m0NC(jgrid,jxi(ihq),
     3                                1,1,alpha,beta)
     4                                + c1(ihq)
     5                                * SC2m0NC(jgrid,jxi(ihq)+1,
     6                                1,1,alpha,beta) ) ) / TR
                                 CLph(ihq) = CLph(ihq) + aem
     1                                * ( - damp(ihq) * ( c0(ihq)
     2                                * SCLm0NC(jgrid,jxi(ihq),
     3                                1,1,alpha,beta)
     4                                + c1(ihq)
     5                                * SCLm0NC(jgrid,jxi(ihq)+1,
     6                                1,1,alpha,beta) ) ) / TR
                           endif
                        enddo
                     endif
                  endif
*     
*     F2
*     
*     Light Component
*     
*     Photon
                  OpF2(jgrid,3,0,alpha,beta)  =
     1                 OpF2(jgrid,3,0,alpha,beta)
     2                 + ( bq(1) + bq(2) + bq(3) ) * C2ph(3)
*     Singlet
                  OpF2(jgrid,3,1,alpha,beta)  =
     1                 OpF2(jgrid,3,1,alpha,beta)
     2                 + ( bq(1) + bq(2) + bq(3) ) * C2ns(3) / nf
*     T3
                  OpF2(jgrid,3,9,alpha,beta)  =
     1                 OpF2(jgrid,3,9,alpha,beta)
     2                 + fr3 * ( bq(2) - bq(1) ) * C2ns(3) / 2d0
*     T8
                  OpF2(jgrid,3,10,alpha,beta) =
     1                 OpF2(jgrid,3,10,alpha,beta)
     2                 + ( bq(1) + bq(2) - 2d0 * bq(3) ) * C2ns(3) / 6d0
*     T15, T24, T35
                  do ik=4,nf
                     OpF2(jgrid,3,7+ik,alpha,beta) =
     1                    OpF2(jgrid,3,7+ik,alpha,beta)
     2                    + ( bq(1) + bq(2) + bq(3) ) * C2ns(3)
     3                    / ik / ( ik -1 )
                  enddo
*     
*     Charm Component
*     
*     Photon
                  OpF2(jgrid,4,0,alpha,beta)  =
     1                 OpF2(jgrid,4,0,alpha,beta) + bq(4) * C2ph(4)
*     Singlet
                  OpF2(jgrid,4,1,alpha,beta)  =
     1                 OpF2(jgrid,4,1,alpha,beta)
     2                 + bq(4) * C2ns(4) / nf
*     T15
                  OpF2(jgrid,4,11,alpha,beta) =
     1                 OpF2(jgrid,4,11,alpha,beta)
     2                  - bq(4) * C2ns(4) / 4d0
*     T24, T35
                  do ik=5,nf
                     OpF2(jgrid,4,7+ik,alpha,beta) =
     1                    OpF2(jgrid,4,7+ik,alpha,beta)
     2                     + bq(4) * C2ns(4) / ik / ( ik -1 )
                  enddo
*     
*     Bottom Component
*     
*     Photon
                  OpF2(jgrid,5,0,alpha,beta)  =
     1                 OpF2(jgrid,5,0,alpha,beta) + bq(5) * C2ph(5)
*     Singlet
                  OpF2(jgrid,5,1,alpha,beta)  =
     1                 OpF2(jgrid,5,1,alpha,beta) + bq(5) * C2ns(5) / nf
*     T24
                  OpF2(jgrid,5,12,alpha,beta) =
     1                 OpF2(jgrid,5,12,alpha,beta)
     2                  - bq(5) * C2ns(5) / 5d0
*     T35
                  do ik=6,nf
                     OpF2(jgrid,5,7+ik,alpha,beta) =
     1                    OpF2(jgrid,5,7+ik,alpha,beta)
     2                     + bq(5) * C2ns(5)/ ik / ( ik - 1 )
                  enddo
*     
*     Top Component
*     
*     Photon
                  OpF2(jgrid,6,0,alpha,beta)  =
     1                 OpF2(jgrid,6,0,alpha,beta) + bq(6) * C2ph(6)
*     Singlet
                  OpF2(jgrid,6,1,alpha,beta)  =
     1                 OpF2(jgrid,6,1,alpha,beta) + bq(6) * C2ns(6) / nf
*     T35
                  OpF2(jgrid,6,13,alpha,beta) =
     1                 OpF2(jgrid,6,13,alpha,beta)
     2                  - bq(6) * C2ns(6) / 6d0
*
*     Total
*
                  do ipdf=0,13
                     OpF2(jgrid,7,ipdf,alpha,beta) = 0d0
                     do ihq=3,6
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
*     Photon
                  OpFL(jgrid,3,0,alpha,beta)  =
     1                 OpFL(jgrid,3,0,alpha,beta)
     2                 + ( bq(1) + bq(2) + bq(3) ) * CLph(3)
*     Singlet
                  OpFL(jgrid,3,1,alpha,beta)  =
     1                 OpFL(jgrid,3,1,alpha,beta)
     2                 + ( bq(1) + bq(2) + bq(3) ) * CLns(3) / nf
*     T3
                  OpFL(jgrid,3,9,alpha,beta)  =
     1                 OpFL(jgrid,3,9,alpha,beta)
     2                 + fr3 * ( bq(2) - bq(1) ) * CLns(3) / 2d0
*     T8
                  OpFL(jgrid,3,10,alpha,beta) =
     1                 OpFL(jgrid,3,10,alpha,beta)
     2                 + ( bq(1) + bq(2) - 2d0 * bq(3) ) * CLns(3) / 6d0
*     T15, T24, T35
                  do ik=4,nf
                     OpFL(jgrid,3,7+ik,alpha,beta) =
     1                    OpFL(jgrid,3,7+ik,alpha,beta)
     2                    + ( bq(1) + bq(2) + bq(3) ) * CLns(3)
     3                    / ik / ( ik -1 )
                  enddo
*     
*     Charm Component
*     
*     Photon
                  OpFL(jgrid,4,0,alpha,beta)  =
     1                 OpFL(jgrid,4,0,alpha,beta) + bq(4) * CLph(4)
*     Singlet
                  OpFL(jgrid,4,1,alpha,beta)  =
     1                 OpFL(jgrid,4,1,alpha,beta)
     2                 + bq(4) * CLns(4) / nf
*     T15
                  OpFL(jgrid,4,11,alpha,beta) =
     1                 OpFL(jgrid,4,11,alpha,beta)
     2                  - bq(4) * CLns(4) / 4d0
*     T24, T35
                  do ik=5,nf
                     OpFL(jgrid,4,7+ik,alpha,beta) =
     1                    OpFL(jgrid,4,7+ik,alpha,beta)
     2                     + bq(4) * CLns(4) / ik / ( ik -1 )
                  enddo
*     
*     Bottom Component
*     
*     Photon
                  OpFL(jgrid,5,0,alpha,beta)  =
     1                 OpFL(jgrid,5,0,alpha,beta) + bq(5) * CLph(5)
*     Singlet
                  OpFL(jgrid,5,1,alpha,beta)  =
     1                 OpFL(jgrid,5,1,alpha,beta) + bq(5) * CLns(5) / nf
*     T24
                  OpFL(jgrid,5,12,alpha,beta) =
     1                 OpFL(jgrid,5,12,alpha,beta)
     2                  - bq(5) * CLns(5) / 5d0
*     T35
                  do ik=6,nf
                     OpFL(jgrid,5,7+ik,alpha,beta) =
     1                    OpFL(jgrid,5,7+ik,alpha,beta)
     2                     + bq(5) * CLns(5)/ ik / ( ik - 1 )
                  enddo
*     
*     Top Component
*     
*     Photon
                  OpFL(jgrid,6,0,alpha,beta)  =
     1                 OpFL(jgrid,6,0,alpha,beta) + bq(6) * CLph(6)
*     Singlet
                  OpFL(jgrid,6,1,alpha,beta)  =
     1                 OpFL(jgrid,6,1,alpha,beta) + bq(6) * CLns(6) / nf
*     T35
                  OpFL(jgrid,6,13,alpha,beta) =
     1                 OpFL(jgrid,6,13,alpha,beta)
     2                  - bq(6) * CLns(6) / 6d0
*
*     Total
*
                  do ipdf=0,13
                     OpFL(jgrid,7,ipdf,alpha,beta) = 0d0
                     do ihq=3,6
                        OpFL(jgrid,7,ipdf,alpha,beta) = 
     1                       OpFL(jgrid,7,ipdf,alpha,beta)
     2                       + OpFL(jgrid,ihq,ipdf,alpha,beta)
                     enddo
                  enddo
*     
*     F3 (Only for neutral current processes and always ZM)
*     
                  if(ProcessDIS.eq."NC")then
*     
*     Light Component
*     
*     Valence
                     OpF3(jgrid,3,3,alpha,beta)  =
     1                    OpF3(jgrid,3,3,alpha,beta)
     2                    + ( dq(1) + dq(2) + dq(3) ) * C3ns(3) / nf
*     V3
                     OpF3(jgrid,3,4,alpha,beta)  =
     1                    OpF3(jgrid,3,4,alpha,beta)
     2                    + fr3 * ( dq(2) - dq(1) ) * C3ns(3) / 2d0
*     V8
                     OpF3(jgrid,3,5,alpha,beta)  =
     1                    OpF3(jgrid,3,5,alpha,beta) 
     2                    + ( dq(1) + dq(2) - 2d0 * dq(3) ) * C3ns(3)
     3                    / 6d0
*     V15, V24, V35
                     do ik=4,nf
                        OpF3(jgrid,3,2+ik,alpha,beta) =
     1                       OpF3(jgrid,3,2+ik,alpha,beta)
     2                       + ( dq(1) + dq(2) + dq(3) ) * C3ns(3)
     3                       / ik / ( ik -1 )
                     enddo
*     
*     Charm Component
*     
                     if(nf.ge.4)then
*     Valence
                        OpF3(jgrid,4,3,alpha,beta)  =
     1                       OpF3(jgrid,4,3,alpha,beta)
     2                        + dq(4) * C3ns(4) / nf
*     V15
                        OpF3(jgrid,4,6,alpha,beta) =
     1                       OpF3(jgrid,4,6,alpha,beta)
     2                       - dq(4) * C3ns(4) / 4d0
*     V24, V35
                        do ik=5,nf
                           OpF3(jgrid,4,2+ik,alpha,beta) =
     1                          OpF3(jgrid,4,2+ik,alpha,beta)
     2                           + dq(4) * C3ns(4) / ik / ( ik -1 )
                        enddo
                     endif
*     
*     Bottom Component
*     
                     if(nf.ge.5)then
*     Valence
                        OpF3(jgrid,5,3,alpha,beta)  =
     1                       OpF3(jgrid,5,3,alpha,beta)
     2                       + dq(5) * C3ns(5) / nf
*     V24
                        OpF3(jgrid,5,7,alpha,beta) =
     1                       OpF3(jgrid,5,7,alpha,beta)
     2                       - dq(5) * C3ns(5) / 5d0
*     V35
                        do ik=6,nf
                           OpF3(jgrid,5,2+ik,alpha,beta) =
     1                          OpF3(jgrid,5,2+ik,alpha,beta)
     2                          + dq(5) * C3ns(5) / ik / ( ik - 1 )
                        enddo
                     endif
*     
*     Top Component
*     
                     if(nf.ge.6)then
*     Valence
                        OpF3(jgrid,6,3,alpha,beta) =
     1                       OpF3(jgrid,6,3,alpha,beta)
     2                       + dq(6) * C3ns(6) / nf
*     V35
                        OpF3(jgrid,6,8,alpha,beta) =
     1                       OpF3(jgrid,6,8,alpha,beta)
     2                       - dq(6) * C3ns(6) / 6d0
                     endif
*     
*     Total
*     
                     do ipdf=0,13
                        OpF3(jgrid,7,ipdf,alpha,beta) = 0d0
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
*     Initialize coefficient functions
*
                  do ihq=3,6
                     C2ph(ihq) = 0d0
                     C2ns(ihq) = 0d0
                     CLph(ihq) = 0d0
                     CLns(ihq) = 0d0
                     C3ph(ihq) = 0d0
                     C3ns(ihq) = 0d0
                  enddo
*
*     Construct coefficient functions according to the scheme chosen
*
                  if(MassScheme.eq."ZM-VFNS")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,nf,1,1,alpha,beta) / TR
                     C2ns(3) = aem * SC2zm(jgrid,nf,3,1,alpha,beta) / CF
                     CLph(3) = aem * SCLzm(jgrid,nf,1,1,alpha,beta) / TR
                     CLns(3) = aem * SCLzm(jgrid,nf,3,1,alpha,beta) / CF
                     C3ns(3) = aem * SC3zm(jgrid,nf,4,1,alpha,beta) / CF
*     Heavy quark coefficient functions
                     if(nf.gt.3)then
                        do ihq=4,nf
                           C2ph(ihq) = C2ph(3)
                           C2ns(ihq) = C2ns(3)
                           CLph(ihq) = CLph(3)
                           CLns(ihq) = CLns(3)
                           C3ns(ihq) = C3ns(3)
                        enddo
                     endif
                  elseif(MassScheme.eq."FFN0")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     C2ns(3) = aem * SC2zm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     CLph(3) = aem * SCLzm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     CLns(3) = aem * SCLzm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     C3ns(3) = aem * SC3zm(jgrid,Nf_FF,4,1,alpha,beta)
     1                    / CF
*     Heavy quark coefficient functions
                     if(Nf_FF.gt.3)then
                        do ihq=4,Nf_FF
                           C2ph(ihq) = C2ph(3)
                           C2ns(ihq) = C2ns(3)
                           CLph(ihq) = CLph(3)
                           CLns(ihq) = CLns(3)
                           C3ns(ihq) = C3ns(3)
                        enddo
                     endif
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           if(Q2.ge.m2th(ihq))then
                              C2ph(ihq) = C2ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2m0CC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2m0CC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              C2ns(ihq) = C2ns(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2m0CC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2m0CC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) ) / CF
                              CLph(ihq) = CLph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLm0CC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLm0CC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              CLns(ihq) = CLns(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLm0CC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLm0CC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) ) / CF
                              C3ph(ihq) = C3ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC3m0CC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3m0CC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              C3ns(ihq) = C3ns(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC3m0CC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3m0CC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) ) / CF
                           endif
                        enddo
                     endif
                  elseif(MassScheme.eq."FFNS")then
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     C2ns(3) = aem * SC2zm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     CLph(3) = aem * SCLzm(jgrid,Nf_FF,1,1,alpha,beta)
     1                    / TR
                     CLns(3) = aem * SCLzm(jgrid,Nf_FF,3,1,alpha,beta)
     1                    / CF
                     C3ns(3) = aem * SC3zm(jgrid,Nf_FF,4,1,alpha,beta)
     1                    / CF
*     Heavy quark coefficient functions
                     if(Nf_FF.gt.3)then
                        do ihq=4,Nf_FF
                           C2ph(ihq) = C2ph(3)
                           C2ns(ihq) = C2ns(3)
                           CLph(ihq) = CLph(3)
                           CLns(ihq) = CLns(3)
                           C3ns(ihq) = C3ns(3)
                        enddo
                     endif
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           if(W2.ge.m2th(ihq))then
                              C2ph(ihq) = C2ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2mCC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2mCC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              C2ns(ihq) = C2ns(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2mCC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2mCC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) ) / CF
                              CLph(ihq) = CLph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLmCC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLmCC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              CLns(ihq) = CLns(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLmCC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLmCC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) ) / CF
                              C3ph(ihq) = C3ph(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC3mCC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3mCC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) ) / TR
                              C3ns(ihq) = C3ns(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC3mCC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3mCC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) ) / CF
                           endif
                        enddo
                     endif
                  elseif(MassScheme(1:5).eq."FONLL")then
*     Heavy quark coefficient functions
                     do ihq=4,6
                        C2gm(ihq) = 0d0
                        C2qm(ihq) = 0d0
                        CLgm(ihq) = 0d0
                        CLqm(ihq) = 0d0
                        C3gm(ihq) = 0d0
                        C3qm(ihq) = 0d0
                     enddo
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           if(W2.ge.m2th(ihq))then
                              C2gm(ihq) = C2gm(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2mCC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2mCC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) )
                              C2qm(ihq) = C2qm(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC2mCC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2mCC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) )
                              CLgm(ihq) = CLgm(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLmCC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLmCC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) )
                              CLqm(ihq) = CLqm(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SCLmCC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLmCC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) )
                              C3gm(ihq) = C3gm(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC3mCC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3mCC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) )
                              C3qm(ihq) = C3qm(ihq) + aem
     1                             * ( c0(ihq)
     2                             * SC3mCC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3mCC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) )
                           endif
                           if(Q2.ge.m2th(ihq))then
                              C2gm(ihq) = C2gm(ihq) - damp(ihq)
     1                             * aem * ( c0(ihq)
     2                             * SC2m0CC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2m0CC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) )
                              C2qm(ihq) = C2qm(ihq) - damp(ihq)
     1                             * aem * ( c0(ihq)
     2                             * SC2m0CC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC2m0CC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) )
                              CLgm(ihq) = CLgm(ihq) - damp(ihq)
     1                             * aem * ( c0(ihq)
     2                             * SCLm0CC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLm0CC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) )
                              CLqm(ihq) = CLqm(ihq) - damp(ihq)
     1                             * aem * ( c0(ihq)
     2                             * SCLm0CC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SCLm0CC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) )
                              C3gm(ihq) = C3gm(ihq) - damp(ihq)
     1                             * aem * ( c0(ihq)
     2                             * SC3m0CC(jgrid,jxi(ihq),
     3                             1,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3m0CC(jgrid,jxi(ihq)+1,
     6                             1,1,alpha,beta) )
                              C3qm(ihq) = C3qm(ihq) - damp(ihq)
     1                             * aem * ( c0(ihq)
     2                             * SC3m0CC(jgrid,jxi(ihq),
     3                             3,1,alpha,beta)
     4                             + c1(ihq)
     5                             * SC3m0CC(jgrid,jxi(ihq)+1,
     6                             3,1,alpha,beta) )
                           endif
                        enddo
                     endif
*
*     ZM contribution
*
*     Light coefficient functions
                     C2ph(3) = aem * SC2zm(jgrid,nf,1,1,alpha,beta) / TR
                     C2ns(3) = aem * SC2zm(jgrid,nf,3,1,alpha,beta) / CF
                     CLph(3) = aem * SCLzm(jgrid,nf,1,1,alpha,beta) / TR
                     CLns(3) = aem * SCLzm(jgrid,nf,3,1,alpha,beta) / CF
                     C3ns(3) = aem * SC3zm(jgrid,nf,4,1,alpha,beta) / CF
*     Heavy quark coefficient functions
                     if(nf.gt.3)then
                        do ihq=4,nf
                           C2ph(ihq) = damp(ihq) * C2ph(3)
                           C2ns(ihq) = damp(ihq) * C2ns(3)
                           CLph(ihq) = damp(ihq) * CLph(3)
                           CLns(ihq) = damp(ihq) * CLns(3)
                           C3ns(ihq) = damp(ihq) * C3ns(3)
                        enddo
                     endif
*
*     Include massive contributions
*
                     if(Nf_FF.lt.6)then
                        do ihq=Nf_FF+1,6
                           C2ph(ihq) = C2ph(ihq) + C2gm(ihq)
                           C2ns(ihq) = C2ns(ihq) + C2qm(ihq)
                           CLph(ihq) = CLph(ihq) + CLgm(ihq)
                           CLns(ihq) = CLns(ihq) + CLqm(ihq)
                           C3ph(ihq) = C3ph(ihq) + C3gm(ihq)
                           C3ns(ihq) = C3ns(ihq) + C3qm(ihq)
                        enddo
                     endif
                  endif
*
*     The following corrections are computed assuming a diagonal CKM matrix.
*
*     Coefficients
*
                  etap = ( e2u + e2d ) / 2d0
                  etam = ( e2u - e2d ) / 2d0
*
*     F2
*
*     Light Component
*
*     Photon
                  OpF2(jgrid,3,0,alpha,beta)  =
     1                 OpF2(jgrid,3,0,alpha,beta)
     2                 + 2d0 * etap * C2ph(3)
*     Singlet
                  OpF2(jgrid,3,1,alpha,beta)  =
     1                 OpF2(jgrid,3,1,alpha,beta)
     2                 + etap * C2ns(3) / 3d0
*     Valence
                  OpF2(jgrid,3,3,alpha,beta)  =
     1                 OpF2(jgrid,3,3,alpha,beta)
     2                 - ipr * etam * C2ns(3) / 3d0
*     V3
                  OpF2(jgrid,3,4,alpha,beta)  =
     1                 OpF2(jgrid,3,4,alpha,beta)
     2                 - ipr * fr3 * etap * C2ns(3)
*     T3
                  OpF2(jgrid,3,9,alpha,beta)  =
     1                 OpF2(jgrid,3,9,alpha,beta)
     2                 + fr3 * etam * C2ns(3)
*     V8
                  OpF2(jgrid,3,5,alpha,beta)  =
     1                 OpF2(jgrid,3,5,alpha,beta)
     2                 - ipr * etam * C2ns(3) / 3d0
*     T8
                  OpF2(jgrid,3,10,alpha,beta)  =
     1                 OpF2(jgrid,3,10,alpha,beta)
     2                 + etap * C2ns(3) / 3d0
*     V15
                  OpF2(jgrid,3,6,alpha,beta)  =
     1                 OpF2(jgrid,3,6,alpha,beta)
     2                 - ipr * etam * C2ns(3) / 6d0
*     T15
                  OpF2(jgrid,3,11,alpha,beta)  =
     1                 OpF2(jgrid,3,11,alpha,beta)
     2                 + etap * C2ns(3) / 6d0
*     V24
                  OpF2(jgrid,3,7,alpha,beta)  =
     1                 OpF2(jgrid,3,7,alpha,beta)
     2                 - ipr * etam * C2ns(3) / 10d0
*     T24
                  OpF2(jgrid,3,12,alpha,beta)  =
     1                 OpF2(jgrid,3,12,alpha,beta)
     2                 + etap * C2ns(3) / 10d0
*     V35
                  OpF2(jgrid,3,8,alpha,beta)  =
     1                 OpF2(jgrid,3,8,alpha,beta)
     2                 - ipr * etam * C2ns(3) / 15d0
*     T35
                  OpF2(jgrid,3,13,alpha,beta)  =
     1                 OpF2(jgrid,3,13,alpha,beta)
     2                 + etap * C2ns(3) / 15d0
*
*     Charm Component
*
*     Photon
                  OpF2(jgrid,4,0,alpha,beta)  =
     1                 OpF2(jgrid,4,0,alpha,beta)
     2                 + 2d0 * etap * C2ph(4)
*     Singlet
                  OpF2(jgrid,4,1,alpha,beta)  =
     1                 OpF2(jgrid,4,1,alpha,beta)
     2                 + etap * C2ns(4) / 3d0
*     Valence
                  OpF2(jgrid,4,3,alpha,beta)  =
     1                 OpF2(jgrid,4,3,alpha,beta)
     2                 - ipr * etam * C2ns(4) / 3d0
*     V8
                  OpF2(jgrid,4,5,alpha,beta)  =
     1                 OpF2(jgrid,4,5,alpha,beta)
     2                 - ipr * e2d * C2ns(4) / 3d0
*     T8
                  OpF2(jgrid,4,10,alpha,beta)  =
     1                 OpF2(jgrid,4,10,alpha,beta)
     2                 - e2d * C2ns(4) / 3d0
*     V15
                  OpF2(jgrid,4,6,alpha,beta)  =
     1                 OpF2(jgrid,4,6,alpha,beta)
     2                 + ipr * ( e2d + 3d0 * e2u ) * C2ns(4) / 12d0
*     T15
                  OpF2(jgrid,4,11,alpha,beta)  =
     1                 OpF2(jgrid,4,11,alpha,beta)
     2                 + ( e2d - 3d0 * e2u ) * C2ns(4) / 12d0
*     V24
                  OpF2(jgrid,4,7,alpha,beta)  =
     1                 OpF2(jgrid,4,7,alpha,beta)
     2                 - ipr * etam * C2ns(4) / 10d0
*     T24
                  OpF2(jgrid,4,12,alpha,beta)  =
     1                 OpF2(jgrid,4,12,alpha,beta)
     2                 + etap * C2ns(4) / 10d0
*     V35
                  OpF2(jgrid,4,8,alpha,beta)  =
     1                 OpF2(jgrid,4,8,alpha,beta)
     2                 - ipr * etam * C2ns(4) / 15d0
*     T35
                  OpF2(jgrid,4,13,alpha,beta)  =
     1                 OpF2(jgrid,4,13,alpha,beta)
     2                 + etap * C2ns(4) / 15d0
*
                  if(MassScheme(1:5).eq."FONLL")then
*     Subtract the spurious charm contribution
                     if(Nf_FF.le.3)then
*     Singlet
                        OpF2(jgrid,4,1,alpha,beta)  =
     1                       OpF2(jgrid,4,1,alpha,beta)
     2                       - etap * C2qm(4) / 6d0
*     Valence
                        OpF2(jgrid,4,3,alpha,beta)  =
     1                       OpF2(jgrid,4,3,alpha,beta)
     2                       + ipr * etam * C2qm(4) / 6d0
*     V15
                        OpF2(jgrid,4,6,alpha,beta)  =
     1                       OpF2(jgrid,4,6,alpha,beta)
     2                       - ipr * etam * C2qm(4) / 4d0
*     T15
                        OpF2(jgrid,4,11,alpha,beta) =
     1                       OpF2(jgrid,4,11,alpha,beta)
     2                       + etap * C2qm(4) / 4d0
*     V24
                        OpF2(jgrid,4,7,alpha,beta)  =
     1                       OpF2(jgrid,4,7,alpha,beta)
     2                       + ipr * etam * C2qm(4) / 20d0
*     T24
                        OpF2(jgrid,4,12,alpha,beta) =
     1                       OpF2(jgrid,4,12,alpha,beta)
     2                       - etap * C2qm(4) / 20d0
*     V35
                        OpF2(jgrid,4,8,alpha,beta)  =
     1                       OpF2(jgrid,4,8,alpha,beta)
     2                       + ipr * etam * C2qm(4) / 30d0
*     T35
                        OpF2(jgrid,4,13,alpha,beta) =
     1                       OpF2(jgrid,4,13,alpha,beta)
     2                       - etap * C2qm(4) / 30d0
                     endif
                  endif
*
*     Bottom component (no corrections)
*
*
*     Top component (no corrections)
*
*     Photon
                  OpF2(jgrid,6,0,alpha,beta)  =
     1                 OpF2(jgrid,6,0,alpha,beta)
     2                 + 2d0 * etap * C2ph(6)
*     Singlet
                  OpF2(jgrid,6,1,alpha,beta)  =
     1                 OpF2(jgrid,6,1,alpha,beta)
     2                 + etap * C2ns(6) / 3d0
*     Valence
                  OpF2(jgrid,6,3,alpha,beta)  =
     1                 OpF2(jgrid,6,3,alpha,beta)
     2                 - ipr * etam * C2ns(6) / 3d0
*     V24
                  OpF2(jgrid,6,7,alpha,beta)  =
     1                 OpF2(jgrid,6,7,alpha,beta)
     2                 - ipr * e2d * C2ns(6) / 5d0
*     T24
                  OpF2(jgrid,6,12,alpha,beta)  =
     1                 OpF2(jgrid,6,12,alpha,beta)
     2                 - e2d * C2ns(6) / 5d0
*     V35
                  OpF2(jgrid,6,8,alpha,beta)  =
     1                 OpF2(jgrid,6,8,alpha,beta)
     2                 + ipr * ( e2d + 5d0 * e2u ) * C2ns(6) / 30d0
*     T35
                  OpF2(jgrid,6,13,alpha,beta)  =
     1                 OpF2(jgrid,6,13,alpha,beta)
     2                 + ( e2d - 5d0 * e2u ) * C2ns(6) / 30d0
*
                  if(MassScheme(1:5).eq."FONLL")then
*     Subtract the spurious top contribution
                     if(Nf_FF.le.5)then
*     Singlet
                        OpF2(jgrid,6,1,alpha,beta)  =
     1                       OpF2(jgrid,6,1,alpha,beta)
     2                       - etap * C2qm(6) / 6d0
*     Valence
                        OpF2(jgrid,6,3,alpha,beta)  =
     1                       OpF2(jgrid,6,3,alpha,beta)
     2                       + ipr * etam * C2qm(6) / 6d0
*     V35
                        OpF2(jgrid,6,8,alpha,beta)  =
     1                       OpF2(jgrid,6,8,alpha,beta)
     2                       - ipr * etam * C2qm(6) / 6d0
*     T35
                        OpF2(jgrid,6,13,alpha,beta) =
     1                       OpF2(jgrid,6,13,alpha,beta)
     2                       + etap * C2qm(6) / 6d0
                     endif
*     Subtract the spurious bottom contribution
                     if(Nf_FF.le.4)then
*     Singlet
                        OpF2(jgrid,6,1,alpha,beta)  =
     1                       OpF2(jgrid,6,1,alpha,beta)
     2                       - etap * C2qm(6) / 6d0
*     Valence
                        OpF2(jgrid,6,3,alpha,beta)  =
     1                       OpF2(jgrid,6,3,alpha,beta)
     2                       - ipr * etam * C2qm(6) / 6d0
*     V24
                        OpF2(jgrid,6,7,alpha,beta)  =
     1                       OpF2(jgrid,6,7,alpha,beta)
     2                       + ipr * etam * C2qm(6) / 5d0
*     T24
                        OpF2(jgrid,6,12,alpha,beta) =
     1                       OpF2(jgrid,6,12,alpha,beta)
     2                       + etap * C2qm(6) / 5d0
*     V35
                        OpF2(jgrid,6,8,alpha,beta)  =
     1                       OpF2(jgrid,6,8,alpha,beta)
     2                       - ipr * etam * C2qm(6) / 30d0
*     T35
                        OpF2(jgrid,6,13,alpha,beta) =
     1                       OpF2(jgrid,6,13,alpha,beta)
     2                       - etap * C2qm(6) / 30d0
                     endif
                  endif
*
*     Total
*
                  do ipdf=0,13
                     OpF2(jgrid,7,ipdf,alpha,beta) = 0d0
                     do ihq=3,6
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
*     Photon
                  OpFL(jgrid,3,0,alpha,beta)  =
     1                 OpFL(jgrid,3,0,alpha,beta)
     2                 + 2d0 * etap * CLph(3)
*     Singlet
                  OpFL(jgrid,3,1,alpha,beta)  =
     1                 OpFL(jgrid,3,1,alpha,beta)
     2                 + etap * CLns(3) / 3d0
*     Valence
                  OpFL(jgrid,3,3,alpha,beta)  =
     1                 OpFL(jgrid,3,3,alpha,beta)
     2                 - ipr * etam * CLns(3) / 3d0
*     V3
                  OpFL(jgrid,3,4,alpha,beta)  =
     1                 OpFL(jgrid,3,4,alpha,beta)
     2                 - ipr * fr3 * etap * CLns(3)
*     T3
                  OpFL(jgrid,3,9,alpha,beta)  =
     1                 OpFL(jgrid,3,9,alpha,beta)
     2                 + fr3 * etam * CLns(3)
*     V8
                  OpFL(jgrid,3,5,alpha,beta)  =
     1                 OpFL(jgrid,3,5,alpha,beta)
     2                 - ipr * etam * CLns(3) / 3d0
*     T8
                  OpFL(jgrid,3,10,alpha,beta)  =
     1                 OpFL(jgrid,3,10,alpha,beta)
     2                 + etap * CLns(3) / 3d0
*     V15
                  OpFL(jgrid,3,6,alpha,beta)  =
     1                 OpFL(jgrid,3,6,alpha,beta)
     2                 - ipr * etam * CLns(3) / 6d0
*     T15
                  OpFL(jgrid,3,11,alpha,beta)  =
     1                 OpFL(jgrid,3,11,alpha,beta)
     2                 + etap * CLns(3) / 6d0
*     V24
                  OpFL(jgrid,3,7,alpha,beta)  =
     1                 OpFL(jgrid,3,7,alpha,beta)
     2                 - ipr * etam * CLns(3) / 10d0
*     T24
                  OpFL(jgrid,3,12,alpha,beta)  =
     1                 OpFL(jgrid,3,12,alpha,beta)
     2                 + etap * CLns(3) / 10d0
*     V35
                  OpFL(jgrid,3,8,alpha,beta)  =
     1                 OpFL(jgrid,3,8,alpha,beta)
     2                 - ipr * etam * CLns(3) / 15d0
*     T35
                  OpFL(jgrid,3,13,alpha,beta)  =
     1                 OpFL(jgrid,3,13,alpha,beta)
     2                 + etap * CLns(3) / 15d0
*
*     Charm Component
*
*     Photon
                  OpFL(jgrid,4,0,alpha,beta)  =
     1                 OpFL(jgrid,4,0,alpha,beta)
     2                 + 2d0 * etap * CLph(4)
*     Singlet
                  OpFL(jgrid,4,1,alpha,beta)  =
     1                 OpFL(jgrid,4,1,alpha,beta)
     2                 + etap * CLns(4) / 3d0
*     Valence
                  OpFL(jgrid,4,3,alpha,beta)  =
     1                 OpFL(jgrid,4,3,alpha,beta)
     2                 - ipr * etam * CLns(4) / 3d0
*     V8
                  OpFL(jgrid,4,5,alpha,beta)  =
     1                 OpFL(jgrid,4,5,alpha,beta)
     2                 - ipr * e2d * CLns(4) / 3d0
*     T8
                  OpFL(jgrid,4,10,alpha,beta)  =
     1                 OpFL(jgrid,4,10,alpha,beta)
     2                 - e2d * CLns(4) / 3d0
*     V15
                  OpFL(jgrid,4,6,alpha,beta)  =
     1                 OpFL(jgrid,4,6,alpha,beta)
     2                 + ipr * ( e2d + 3d0 * e2u ) * CLns(4) / 12d0
*     T15
                  OpFL(jgrid,4,11,alpha,beta)  =
     1                 OpFL(jgrid,4,11,alpha,beta)
     2                 + ( e2d - 3d0 * e2u ) * CLns(4) / 12d0
*     V24
                  OpFL(jgrid,4,7,alpha,beta)  =
     1                 OpFL(jgrid,4,7,alpha,beta)
     2                 - ipr * etam * CLns(4) / 10d0
*     T24
                  OpFL(jgrid,4,12,alpha,beta)  =
     1                 OpFL(jgrid,4,12,alpha,beta)
     2                 + etap * CLns(4) / 10d0
*     V35
                  OpFL(jgrid,4,8,alpha,beta)  =
     1                 OpFL(jgrid,4,8,alpha,beta)
     2                 - ipr * etam * CLns(4) / 15d0
*     T35
                  OpFL(jgrid,4,13,alpha,beta)  =
     1                 OpFL(jgrid,4,13,alpha,beta)
     2                 + etap * CLns(4) / 15d0
*
                  if(MassScheme(1:5).eq."FONLL")then
*     Subtract the spurious charm contribution
                     if(Nf_FF.le.3)then
*     Singlet
                        OpFL(jgrid,4,1,alpha,beta)  =
     1                       OpFL(jgrid,4,1,alpha,beta)
     2                       - etap * CLqm(4) / 6d0
*     Valence
                        OpFL(jgrid,4,3,alpha,beta)  =
     1                       OpFL(jgrid,4,3,alpha,beta)
     2                       + ipr * etam * CLqm(4) / 6d0
*     V15
                        OpFL(jgrid,4,6,alpha,beta)  =
     1                       OpFL(jgrid,4,6,alpha,beta)
     2                       - ipr * etam * CLqm(4) / 4d0
*     T15
                        OpFL(jgrid,4,11,alpha,beta) =
     1                       OpFL(jgrid,4,11,alpha,beta)
     2                       + etap * CLqm(4) / 4d0
*     V24
                        OpFL(jgrid,4,7,alpha,beta)  =
     1                       OpFL(jgrid,4,7,alpha,beta)
     2                       + ipr * etam * CLqm(4) / 20d0
*     T24
                        OpFL(jgrid,4,12,alpha,beta) =
     1                       OpFL(jgrid,4,12,alpha,beta)
     2                       - etap * CLqm(4) / 20d0
*     V35
                        OpFL(jgrid,4,8,alpha,beta)  =
     1                       OpFL(jgrid,4,8,alpha,beta)
     2                       + ipr * etam * CLqm(4) / 30d0
*     T35
                        OpFL(jgrid,4,13,alpha,beta) =
     1                       OpFL(jgrid,4,13,alpha,beta)
     2                       - etap * CLqm(4) / 30d0
                     endif
                  endif
*
*     Bottom component (no corrections)
*
*
*     Top component (no corrections)
*
*     Photon
                  OpFL(jgrid,6,0,alpha,beta)  =
     1                 OpFL(jgrid,6,0,alpha,beta)
     2                 + 2d0 * etap * CLph(6)
*     Singlet
                  OpFL(jgrid,6,1,alpha,beta)  =
     1                 OpFL(jgrid,6,1,alpha,beta)
     2                 + etap * CLns(6) / 3d0
*     Valence
                  OpFL(jgrid,6,3,alpha,beta)  =
     1                 OpFL(jgrid,6,3,alpha,beta)
     2                 - ipr * etam * CLns(6) / 3d0
*     V24
                  OpFL(jgrid,6,7,alpha,beta)  =
     1                 OpFL(jgrid,6,7,alpha,beta)
     2                 - ipr * e2d * CLns(6) / 5d0
*     T24
                  OpFL(jgrid,6,12,alpha,beta)  =
     1                 OpFL(jgrid,6,12,alpha,beta)
     2                 - e2d * CLns(6) / 5d0
*     V35
                  OpFL(jgrid,6,8,alpha,beta)  =
     1                 OpFL(jgrid,6,8,alpha,beta)
     2                 + ipr * ( e2d + 5d0 * e2u ) * CLns(6) / 30d0
*     T35
                  OpFL(jgrid,6,13,alpha,beta)  =
     1                 OpFL(jgrid,6,13,alpha,beta)
     2                 + ( e2d - 5d0 * e2u ) * CLns(6) / 30d0
*
                  if(MassScheme(1:5).eq."FONLL")then
*     Subtract the spurious top contribution
                     if(Nf_FF.le.5)then
*     Singlet
                        OpFL(jgrid,6,1,alpha,beta)  =
     1                       OpFL(jgrid,6,1,alpha,beta)
     2                       - etap * CLqm(6) / 6d0
*     Valence
                        OpFL(jgrid,6,3,alpha,beta)  =
     1                       OpFL(jgrid,6,3,alpha,beta)
     2                       + ipr * etam * CLqm(6) / 6d0
*     V35
                        OpFL(jgrid,6,8,alpha,beta)  =
     1                       OpFL(jgrid,6,8,alpha,beta)
     2                       - ipr * etam * CLqm(6) / 6d0
*     T35
                        OpFL(jgrid,6,13,alpha,beta) =
     1                       OpFL(jgrid,6,13,alpha,beta)
     2                       + etap * CLqm(6) / 6d0
                     endif
*     Subtract the spurious bottom contribution
                     if(Nf_FF.le.4)then
*     Singlet
                        OpFL(jgrid,6,1,alpha,beta)  =
     1                       OpFL(jgrid,6,1,alpha,beta)
     2                       - etap * CLqm(6) / 6d0
*     Valence
                        OpFL(jgrid,6,3,alpha,beta)  =
     1                       OpFL(jgrid,6,3,alpha,beta)
     2                       - ipr * etam * CLqm(6) / 6d0
*     V24
                        OpFL(jgrid,6,7,alpha,beta)  =
     1                       OpFL(jgrid,6,7,alpha,beta)
     2                       + ipr * etam * CLqm(6) / 5d0
*     T24
                        OpFL(jgrid,6,12,alpha,beta) =
     1                       OpFL(jgrid,6,12,alpha,beta)
     2                       + etap * CLqm(6) / 5d0
*     V35
                        OpFL(jgrid,6,8,alpha,beta)  =
     1                       OpFL(jgrid,6,8,alpha,beta)
     2                       - ipr * etam * CLqm(6) / 30d0
*     T35
                        OpFL(jgrid,6,13,alpha,beta) =
     1                       OpFL(jgrid,6,13,alpha,beta)
     2                       - etap * CLqm(6) / 30d0
                     endif
                  endif
*
*     Total
*
                  do ipdf=0,13
                     OpFL(jgrid,7,ipdf,alpha,beta) = 0d0
                     do ihq=3,6
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
*     Photon
                  OpF3(jgrid,3,0,alpha,beta)  =
     1                 OpF3(jgrid,3,0,alpha,beta)
     2                 + 2d0 * etap * C3ph(3)
*     Singlet
                  OpF3(jgrid,3,1,alpha,beta)  =
     1                 OpF3(jgrid,3,1,alpha,beta)
     2                  - ipr * etam * C3ns(3) / 3d0
*     Valence
                  OpF3(jgrid,3,3,alpha,beta)  =
     1                 OpF3(jgrid,3,3,alpha,beta)
     2                 + etap * C3ns(3) / 3d0
*     V3
                  OpF3(jgrid,3,4,alpha,beta)  =
     1                 OpF3(jgrid,3,4,alpha,beta)
     2                 + fr3 * etam * C3ns(3)
*     T3
                  OpF3(jgrid,3,9,alpha,beta)  =
     1                 OpF3(jgrid,3,9,alpha,beta)
     2                 - ipr * fr3 * etap * C3ns(3)
*     V8
                  OpF3(jgrid,3,5,alpha,beta)  =
     1                 OpF3(jgrid,3,5,alpha,beta)
     2                 + etap * C3ns(3) / 3d0
*     T8
                  OpF3(jgrid,3,10,alpha,beta)  =
     1                 OpF3(jgrid,3,10,alpha,beta)
     2                 - ipr * etam * C3ns(3) / 3d0
*     V15
                  OpF3(jgrid,3,6,alpha,beta)  =
     1                 OpF3(jgrid,3,6,alpha,beta)
     2                  + etap * C3ns(3) / 6d0
*     T15
                  OpF3(jgrid,3,11,alpha,beta)  =
     1                 OpF3(jgrid,3,11,alpha,beta)
     2                 - ipr * etam* C3ns(3) / 6d0
*     V24
                  OpF3(jgrid,3,7,alpha,beta)  =
     1                 OpF3(jgrid,3,7,alpha,beta)
     2                  + etap * C3ns(3) / 10d0
*     T24
                  OpF3(jgrid,3,12,alpha,beta)  =
     1                 OpF3(jgrid,3,12,alpha,beta)
     2                 - ipr * etam * C3ns(3) / 10d0
*     V35
                  OpF3(jgrid,3,8,alpha,beta)  =
     1                 OpF3(jgrid,3,8,alpha,beta)
     2                  + etap * C3ns(3) / 15d0
*     T35
                  OpF3(jgrid,3,13,alpha,beta)  =
     1                 OpF3(jgrid,3,13,alpha,beta)
     2                 - ipr * etam * C3ns(3) / 15d0
*
*     Charm Component
*
*     Photon
                  OpF3(jgrid,4,0,alpha,beta)  =
     1                 OpF3(jgrid,4,0,alpha,beta)
     2                 + 2d0 * etap * C3ph(4)
*     Singlet
                  OpF3(jgrid,4,1,alpha,beta)  =
     1                 OpF3(jgrid,4,1,alpha,beta)
     2                 - ipr * etam * C3ns(4) / 3d0
*     Valence
                  OpF3(jgrid,4,3,alpha,beta)  =
     1                 OpF3(jgrid,4,3,alpha,beta)
     2                 + etap * C3ns(4) / 3d0
*     V8
                  OpF3(jgrid,4,5,alpha,beta)  =
     1                 OpF3(jgrid,4,5,alpha,beta)
     2                 - e2d * C3ns(4) / 3d0
*     T8
                  OpF3(jgrid,4,10,alpha,beta)  =
     1                 OpF3(jgrid,4,10,alpha,beta)
     2                 - ipr * e2d * C3ns(4) / 3d0
*     V15
                  OpF3(jgrid,4,6,alpha,beta)  =
     1                 OpF3(jgrid,4,6,alpha,beta)
     2                 + ( e2d - 3d0 * e2u ) * C3ns(4) / 12d0
*     T15
                  OpF3(jgrid,4,11,alpha,beta)  =
     1                 OpF3(jgrid,4,11,alpha,beta)
     2                 + ipr * ( e2d + 3d0 * e2u ) * C3ns(4) / 12d0
*     V24
                  OpF3(jgrid,4,7,alpha,beta)  =
     1                 OpF3(jgrid,4,7,alpha,beta)
     2                 + etap * C3ns(4) / 10d0
*     T24
                  OpF3(jgrid,4,12,alpha,beta)  =
     1                 OpF3(jgrid,4,12,alpha,beta)
     2                 - ipr * etam * C3ns(4) / 10d0
*     V35
                  OpF3(jgrid,4,8,alpha,beta)  =
     1                 OpF3(jgrid,4,8,alpha,beta)
     2                 + etap * C3ns(4) / 15d0
*     T35
                  OpF3(jgrid,4,13,alpha,beta)  =
     1                 OpF3(jgrid,4,13,alpha,beta)
     2                 - ipr * etam * C3ns(4) / 15d0
*
                  if(MassScheme(1:5).eq."FONLL")then
*     Subtract the spurious charm contribution
                     if(Nf_FF.le.3)then
*     Singlet
                        OpF3(jgrid,4,1,alpha,beta)  =
     1                       OpF3(jgrid,4,1,alpha,beta)
     2                       + ipr * etam * C3qm(4) / 6d0
*     Valence
                        OpF3(jgrid,4,3,alpha,beta)  =
     1                       OpF3(jgrid,4,3,alpha,beta)
     2                       - etap * C3qm(4) / 6d0
*     V15
                        OpF3(jgrid,4,6,alpha,beta)  =
     1                       OpF3(jgrid,4,6,alpha,beta)
     2                       + etap * C3qm(4) / 4d0
*     T15
                        OpF3(jgrid,4,11,alpha,beta) =
     1                       OpF3(jgrid,4,11,alpha,beta)
     2                       - ipr * etam * C3qm(4) / 4d0
*     V24
                        OpF3(jgrid,4,7,alpha,beta)  =
     1                       OpF3(jgrid,4,7,alpha,beta)
     2                       - etap * C3qm(4) / 20d0
*     T24
                        OpF3(jgrid,4,12,alpha,beta) =
     1                       OpF3(jgrid,4,12,alpha,beta)
     2                       + ipr * etam * C3qm(4) / 20d0
*     V35
                        OpF3(jgrid,4,8,alpha,beta)  =
     1                       OpF3(jgrid,4,8,alpha,beta)
     2                       - etap * C3qm(4) / 30d0
*     T35
                        OpF3(jgrid,4,13,alpha,beta) =
     1                       OpF3(jgrid,4,13,alpha,beta)
     2                       + ipr * etam * C3qm(4) / 30d0
                     endif
                  endif
*
*     Bottom component (no corrections)
*
*
*     Top component (no corrections)
*
*     Photon
                  OpF3(jgrid,6,0,alpha,beta)  =
     1                 OpF3(jgrid,6,0,alpha,beta)
     2                 + 2d0 * etap * C3ph(6)
*     Singlet
                  OpF3(jgrid,6,1,alpha,beta)  =
     1                 OpF3(jgrid,6,1,alpha,beta)
     2                 - ipr * etam * C3ns(6) / 3d0
*     Valence
                  OpF3(jgrid,6,3,alpha,beta)  =
     1                 OpF3(jgrid,6,3,alpha,beta)
     2                 + etap * C3ns(6) / 3d0
*     V24
                  OpF3(jgrid,6,7,alpha,beta)  =
     1                 OpF3(jgrid,6,7,alpha,beta)
     2                 - e2d * C3ns(6) / 5d0
*     T24
                  OpF3(jgrid,6,12,alpha,beta)  =
     1                 OpF3(jgrid,6,12,alpha,beta)
     2                 - ipr * e2d * C3ns(6) / 5d0
*     V35
                  OpF3(jgrid,6,8,alpha,beta)  =
     1                 OpF3(jgrid,6,8,alpha,beta)
     2                 + ( e2d - 5d0 * e2u ) * C3ns(6) / 30d0
*     T35
                  OpF3(jgrid,6,13,alpha,beta)  =
     1                 OpF3(jgrid,6,13,alpha,beta)
     2                 + ipr * ( e2d + 5d0 * e2u ) * C3ns(6) / 30d0
*
                  if(MassScheme(1:5).eq."FONLL")then
*     Subtract the spurious top contribution
                     if(Nf_FF.le.5)then
*     Singlet
                        OpF3(jgrid,6,1,alpha,beta)  =
     1                       OpF3(jgrid,6,1,alpha,beta)
     2                       + ipr * etam * C3qm(6) / 6d0
*     Valence
                        OpF3(jgrid,6,3,alpha,beta)  =
     1                       OpF3(jgrid,6,3,alpha,beta)
     2                       - etap * C3qm(6) / 6d0
*     V35
                        OpF3(jgrid,6,8,alpha,beta)  =
     1                       OpF3(jgrid,6,8,alpha,beta)
     2                       + etap * C3qm(6) / 6d0
*     T35
                        OpF3(jgrid,6,13,alpha,beta) =
     1                       OpF3(jgrid,6,13,alpha,beta)
     2                       - ipr * etam * C3qm(6) / 6d0
                     endif
*     Subtract the spurious bottom contribution
                     if(Nf_FF.le.4)then
*     Singlet
                        OpF3(jgrid,6,1,alpha,beta)  =
     1                       OpF3(jgrid,6,1,alpha,beta)
     2                       - ipr * etam * C3qm(6) / 6d0
*     Valence
                        OpF3(jgrid,6,3,alpha,beta)  =
     1                       OpF3(jgrid,6,3,alpha,beta)
     2                       - etap * C3qm(6) / 6d0
*     V24
                        OpF3(jgrid,6,7,alpha,beta)  =
     1                       OpF3(jgrid,6,7,alpha,beta)
     2                       + etap * C3qm(6) / 5d0
*     T24
                        OpF3(jgrid,6,12,alpha,beta) =
     1                       OpF3(jgrid,6,12,alpha,beta)
     2                       + ipr * etam * C3qm(6) / 5d0
*     V35
                        OpF3(jgrid,6,8,alpha,beta)  =
     1                       OpF3(jgrid,6,8,alpha,beta)
     2                       - etap * C3qm(6) / 30d0
*     T35
                        OpF3(jgrid,6,13,alpha,beta) =
     1                       OpF3(jgrid,6,13,alpha,beta)
     2                       - ipr * etam * C3qm(6) / 30d0
                     endif
                  endif
*
*     Total
*
                  do ipdf=0,13
                     OpF3(jgrid,7,ipdf,alpha,beta) = 0d0
                     do ihq=3,6
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
      return
      end 
