************************************************************************
*
*     JoinOperatorsQCD.f:
*
*     This routine joins the QCD evolution operators.
*
*     QCD evolution basis:
*     0   1   2   3   4   5   6   7   8   9  10  11  12  13
*     gm  Sg   g   V  V3  V8 V15 V24 V35  T3  T8 T15 T24 T35
*
************************************************************************
      subroutine JoinOperatorsQCD(fevQCD)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/m2th.h"
      include "../commons/EvolutionMatrices.h"
**
*     Internal Variables
*
      integer i,j,k,l
      integer nf,nfm
      integer alpha,beta,gamma,delta
      double precision coup,a_QCD!,integralsMatching
      double precision integns(0:nint_max,0:nint_max)
      double precision integsg(2,2,0:nint_max,0:nint_max)
      double precision MQCD(0:13,0:13,0:nint_max,0:nint_max)
      double precision EQCD(0:13,0:13,0:nint_max,0:nint_max)
      double precision fevQCDb(0:13,0:nint_max)
**
*     Input and Output Variables
*
      double precision fevQCD(0:13,0:nint_max)
*
*     Set evolution operators to zero
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            do i=0,13
               do j=0,13
                  EQCD(i,j,alpha,beta) = 0d0
               enddo
            enddo
         enddo
      enddo
*
      if(nfi.ne.nff)then
         do nf=nfi,nff-sgn,sgn
*
*     Construct matching conditions Operator
*
            if(sgn.eq.-1)then
               nfm = nf
            elseif(sgn.eq.1)then
               nfm = nf + 1
            endif
*     Get alphas value at the heavy quark threshold (with nfm active flavours)
            coup = a_QCD(m2th(nfm))
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  integns(alpha,beta) = 0d0
                  do i=1,2
                     do j=1,2
                        integsg(1,1,alpha,beta) = 0d0
                        integsg(1,2,alpha,beta) = 0d0
                        integsg(2,1,alpha,beta) = 0d0
                        integsg(2,2,alpha,beta) = 0d0
                     enddo
                  enddo
               enddo
            enddo
*     Contruct matching conditions at this threshod
            do alpha=0,nin(igrid)
               do beta=alpha,nin(igrid)
                  integns(alpha,beta)     = 1d0
c     1                 integralsMatching(0,beta-alpha,coup,1)
                  integsg(1,1,alpha,beta) = 1d0
c     1                 integralsMatching(0,beta-alpha,coup,2)
                  integsg(1,2,alpha,beta) = 0d0
c     1                 integralsMatching(0,beta-alpha,coup,3)
                  integsg(2,1,alpha,beta) = 0d0
c     1                 integralsMatching(0,beta-alpha,coup,4)
                  integsg(2,2,alpha,beta) = 1d0
c     1                 integralsMatching(0,beta-alpha,coup,5)
               enddo
            enddo
*     
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)-alpha
*     Singlet and Gluon
                  do i=1,2
                     do j=1,2
                        MQCD(i,j,alpha,beta) = 
     1                       integsg(i,j,alpha,beta)
                     enddo
                  enddo
*     Total Valence
                  MQCD(3,3,alpha,beta)   = integns(alpha,beta)
*     V3
                  MQCD(4,4,alpha,beta)   = integns(alpha,beta)
*     V8
                  MQCD(5,5,alpha,beta)   = integns(alpha,beta)
*     V15
                  MQCD(6,6,alpha,beta)   = integns(alpha,beta)
*     V24
                  MQCD(7,7,alpha,beta)   = integns(alpha,beta)
*     V35
                  MQCD(8,8,alpha,beta)   = integns(alpha,beta)
*     T3
                  MQCD(9,9,alpha,beta)   = integns(alpha,beta)
*     T8
                  MQCD(10,10,alpha,beta) = integns(alpha,beta)
*     T15
                  if(nfm.eq.4)then
                     MQCD(11,1,alpha,beta) = integns(alpha,beta) 
     1                    - 3d0 * ( integsg(1,1,alpha,beta) 
     2                    - integns(alpha,beta) )
                     MQCD(11,2,alpha,beta) = 
     1                    - 3d0 * integsg(1,2,alpha,beta) 
                  else
                     MQCD(11,11,alpha,beta) = integns(alpha,beta)
                  endif
*     T24
                  if(nfm.eq.4)then
                     do j=1,2
                        MQCD(12,j,alpha,beta) = 
     1                       integsg(1,j,alpha,beta)
                     enddo
                  elseif(nfm.eq.5)then
                     MQCD(12,1,alpha,beta) = integns(alpha,beta) 
     1                    - 4d0 * ( integsg(1,1,alpha,beta) 
     2                    - integns(alpha,beta) )
                     MQCD(12,2,alpha,beta) = 
     1                    - 4d0 * integsg(1,2,alpha,beta)
                  else
                     MQCD(12,12,alpha,beta) = integns(alpha,beta)
                  endif
*     T35
                  if(nfm.eq.6)then
                     MQCD(12,1,alpha,beta) = integns(alpha,beta) 
     1                    - 5d0 * ( integsg(1,1,alpha,beta) 
     2                    - integns(alpha,beta) )
                     MQCD(12,2,alpha,beta) = 
     1                    - 5d0 * integsg(1,2,alpha,beta)
                  else
                     do j=1,2
                        MQCD(13,j,alpha,beta) = 
     1                       integsg(1,j,alpha,beta)
                     enddo
                  endif
               enddo
            enddo
*
*     Now combine evolution tables for different numbers of active
*     flavours.
*
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
*     Singlet and Gluon
                  do i=1,2
                     do j=1,2
                        do gamma=0,nin(igrid)
                           do delta=0,nin(igrid)
                              do k=1,2
                                 do l=1,2
                                    EQCD(i,j,alpha,beta) = 
     1                              EQCD(i,j,alpha,beta)
     2                              + MQCDsg(nf+sgn,i,k,alpha,gamma)
     3                              * MQCD(k,l,gamma,delta)
     4                              * MQCDsg(nf,l,j,delta,beta)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
*     Total Valence
                  do gamma=0,nin(igrid)
                     do delta=0,nin(igrid)
                        EQCD(3,3,alpha,beta) = 
     1                       EQCD(3,3,alpha,beta)
     2                       + MQCDnsv(nf+sgn,alpha,gamma)
     3                       * MQCD(3,3,gamma,delta)
     4                       * MQCDnsv(nf,delta,beta)
                     enddo
                  enddo
*     V3
                  do gamma=0,nin(igrid)
                     do delta=0,nin(igrid)
                        EQCD(4,4,alpha,beta) = 
     1                       EQCD(4,4,alpha,beta)
     2                       + MQCDnsm(nf+sgn,alpha,gamma)
     3                       * MQCD(4,4,gamma,delta)
     4                       * MQCDnsm(nf,delta,beta)
                     enddo
                  enddo
*     V8
                  do gamma=0,nin(igrid)
                     do delta=0,nin(igrid)
                        EQCD(5,5,alpha,beta) = 
     1                       EQCD(5,5,alpha,beta)
     2                       + MQCDnsm(nf+sgn,alpha,gamma)
     3                       * MQCD(5,5,gamma,delta)
     4                       * MQCDnsm(nf,delta,beta)
                     enddo
                  enddo
*     V15
                  if(nfm.lt.4)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(6,3,alpha,beta) = 
     1                          EQCD(6,3,alpha,beta)
     2                          + MQCDnsv(nf+sgn,alpha,gamma)
     3                          * MQCD(6,6,gamma,delta)
     4                          * MQCDnsv(nf,delta,beta)
                        enddo
                     enddo
                  elseif(nfm.eq.4)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(6,3,alpha,beta) = 
     1                          EQCD(6,3,alpha,beta)
     2                          + MQCDnsm(nf+sgn,alpha,gamma)
     3                          * MQCD(6,6,gamma,delta)
     4                          * MQCDnsv(nf,delta,beta)
                        enddo
                     enddo
                  else
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(6,6,alpha,beta) = 
     1                          EQCD(6,6,alpha,beta)
     2                          + MQCDnsm(nf+sgn,alpha,gamma)
     3                          * MQCD(6,6,gamma,delta)
     4                          * MQCDnsm(nf,delta,beta)
                        enddo
                     enddo
                  endif
*     V24
                  if(nfm.lt.5)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(7,3,alpha,beta) = 
     1                          EQCD(7,3,alpha,beta)
     2                          + MQCDnsv(nf+sgn,alpha,gamma)
     3                          * MQCD(7,7,gamma,delta)
     4                          * MQCDnsv(nf,delta,beta)
                        enddo
                     enddo
                  elseif(nfm.eq.5)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(7,3,alpha,beta) = 
     1                          EQCD(7,3,alpha,beta)
     2                          + MQCDnsm(nf+sgn,alpha,gamma)
     3                          * MQCD(7,7,gamma,delta)
     4                          * MQCDnsv(nf,delta,beta)
                        enddo
                     enddo
                  else
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(7,7,alpha,beta) = 
     1                          EQCD(7,7,alpha,beta)
     2                          + MQCDnsm(nf+sgn,alpha,gamma)
     3                          * MQCD(7,7,gamma,delta)
     4                          * MQCDnsm(nf,delta,beta)
                        enddo
                     enddo
                  endif
*     V35
                  if(nfm.lt.6)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(8,3,alpha,beta) = 
     1                          EQCD(8,3,alpha,beta)
     2                          + MQCDnsv(nf+sgn,alpha,gamma)
     3                          * MQCD(8,8,gamma,delta)
     4                          * MQCDnsv(nf,delta,beta)
                        enddo
                     enddo
                  elseif(nfm.eq.6)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(8,3,alpha,beta) = 
     1                          EQCD(8,3,alpha,beta)
     2                          + MQCDnsm(nf+sgn,alpha,gamma)
     3                          * MQCD(8,8,gamma,delta)
     4                          * MQCDnsv(nf,delta,beta)
                        enddo
                     enddo
                  else
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(8,8,alpha,beta) = 
     1                          EQCD(8,8,alpha,beta)
     2                          + MQCDnsm(nf+sgn,alpha,gamma)
     3                          * MQCD(8,8,gamma,delta)
     4                          * MQCDnsm(nf,delta,beta)
                        enddo
                     enddo
                  endif
*     T3
                  do gamma=0,nin(igrid)
                     do delta=0,nin(igrid)
                        EQCD(9,9,alpha,beta) = 
     1                       EQCD(9,9,alpha,beta)
     2                       + MQCDnsp(nf+sgn,alpha,gamma)
     3                       * MQCD(9,9,gamma,delta)
     4                       * MQCDnsp(nf,delta,beta)
                     enddo
                  enddo
*     T8
                  do gamma=0,nin(igrid)
                     do delta=0,nin(igrid)
                        EQCD(10,10,alpha,beta) = 
     1                       EQCD(10,10,alpha,beta)
     2                       + MQCDnsp(nf+sgn,alpha,gamma)
     3                       * MQCD(10,10,gamma,delta)
     4                       * MQCDnsp(nf,delta,beta)
                     enddo
                  enddo
*     T15
                  if(nfm.lt.4)then
                     do j=1,2
                        do gamma=0,nin(igrid)
                           do k=1,2
                              EQCD(11,j,alpha,beta) = 
     1                             EQCD(11,j,alpha,beta)
     2                             + MQCDsg(nf+sgn,1,k,alpha,gamma)
     3                             * MQCDsg(nf,k,j,gamma,beta)
                           enddo
                        enddo
                     enddo
                  elseif(nfm.eq.4)then
                     do j=1,2
                        do gamma=0,nin(igrid)
                           EQCD(11,j,alpha,beta) = 
     1                          EQCD(11,j,alpha,beta)
     2                          + MQCDnsp(nf+sgn,alpha,gamma)
     3                          * MQCDsg(nf,1,j,gamma,beta)
                        enddo
                     enddo
                  else
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(11,11,alpha,beta) = 
     1                          EQCD(11,11,alpha,beta)
     2                          + MQCDnsp(nf+sgn,alpha,gamma)
     3                          * MQCD(11,11,gamma,delta)
     4                          * MQCDnsp(nf,delta,beta)
                        enddo
                     enddo
                  endif
*     T24
                  if(nfm.lt.5)then
                     do j=1,2
                        do gamma=0,nin(igrid)
                           do k=1,2
                              EQCD(12,j,alpha,beta) = 
     1                             EQCD(12,j,alpha,beta)
     2                             + MQCDsg(nf+sgn,1,k,alpha,gamma)
     3                             * MQCDsg(nf,k,j,gamma,beta)
                           enddo
                        enddo
                     enddo
                  elseif(nfm.eq.5)then
                     do j=1,2
                        do gamma=0,nin(igrid)
                           EQCD(12,j,alpha,beta) = 
     1                          EQCD(12,j,alpha,beta)
     2                          + MQCDnsp(nf+sgn,alpha,gamma)
     3                          * MQCDsg(nf,1,j,gamma,beta)
                        enddo
                     enddo
                  else
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(12,12,alpha,beta) = 
     1                          EQCD(12,12,alpha,beta)
     2                          + MQCDnsp(nf+sgn,alpha,gamma)
     3                          * MQCD(12,12,gamma,delta)
     4                          * MQCDnsp(nf,delta,beta)
                        enddo
                     enddo
                  endif
*     T35
                  if(nfm.lt.6)then
                     do j=1,2
                        do gamma=0,nin(igrid)
                           do k=1,2
                              EQCD(13,j,alpha,beta) = 
     1                             EQCD(13,j,alpha,beta)
     2                             + MQCDsg(nf+sgn,1,k,alpha,gamma)
     3                             * MQCDsg(nf,k,j,gamma,beta)
                           enddo
                        enddo
                     enddo
                  elseif(nfm.eq.6)then
                     do j=1,2
                        do gamma=0,nin(igrid)
                           EQCD(13,j,alpha,beta) = 
     1                          EQCD(13,j,alpha,beta)
     2                          + MQCDnsp(nf+sgn,alpha,gamma)
     3                          * MQCDsg(nf,1,j,gamma,beta)
                        enddo
                     enddo
                  else
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
                           EQCD(13,13,alpha,beta) = 
     1                          EQCD(13,13,alpha,beta)
     2                          + MQCDnsp(nf+sgn,alpha,gamma)
     3                          * MQCD(13,13,gamma,delta)
     4                          * MQCDnsp(nf,delta,beta)
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      else
         do alpha=0,nin(igrid)
            do beta=0,nin(igrid)
*     Singlet and Gluon
               do i=1,2
                  do j=1,2
                     EQCD(i,j,alpha,beta) = 
     1                    MQCDsg(nfi,i,j,alpha,beta)
                  enddo
*     Total Valence
                  EQCD(3,3,alpha,beta) = 
     1                 MQCDnsv(nfi,alpha,beta)
*     V3
                  EQCD(4,4,alpha,beta) = 
     1                 MQCDnsm(nfi,alpha,beta)
*     V8
                  EQCD(5,5,alpha,beta) = 
     1                 MQCDnsm(nfi,alpha,beta)
*     V15
                  if(nfi.lt.4)then
                     EQCD(6,3,alpha,beta) = 
     1                    MQCDnsv(nfi,alpha,beta)
                  else
                     EQCD(6,6,alpha,beta) = 
     1                    MQCDnsm(nfi,alpha,beta)
                  endif
*     V24
                  if(nfi.lt.5)then
                     EQCD(7,3,alpha,beta) = 
     1                    MQCDnsv(nfi,alpha,beta)
                  else
                     EQCD(7,7,alpha,beta) = 
     1                    MQCDnsm(nfi,alpha,beta)
                  endif
*     V35
                  if(nfi.lt.6)then
                     EQCD(8,3,alpha,beta) = 
     1                    MQCDnsv(nfi,alpha,beta)
                  else
                     EQCD(8,8,alpha,beta) = 
     1                    MQCDnsm(nfi,alpha,beta)
                  endif
*     T3
                  EQCD(9,9,alpha,beta) = 
     1                 MQCDnsp(nfi,alpha,beta)
*     T8
                  EQCD(10,10,alpha,beta) = 
     1                 MQCDnsp(nfi,alpha,beta)
               enddo
*     T15
               if(nfi.lt.4)then
                  do j=1,2
                     EQCD(11,j,alpha,beta) = 
     1                    MQCDsg(nfi,1,j,alpha,beta)
                  enddo
               else
                  EQCD(11,11,alpha,beta) = 
     1                 MQCDnsp(nfi,alpha,beta)
               endif
*     T24
               if(nfi.lt.5)then
                  do j=1,2
                     EQCD(12,j,alpha,beta) = 
     1                    MQCDsg(nfi,1,j,alpha,beta)
                  enddo
               else
                  EQCD(12,12,alpha,beta) = 
     1                 MQCDnsp(nfi,alpha,beta)
               endif
*     T25
               if(nfi.lt.6)then
                  do j=1,2
                     EQCD(13,j,alpha,beta) = 
     1                    MQCDsg(nfi,1,j,alpha,beta)
                  enddo
               else
                  EQCD(13,13,alpha,beta) = 
     1                 MQCDnsp(nfi,alpha,beta)
               endif
            enddo
         enddo
      endif
*
*     Covolute input PDFs with the joint evolution operator
*
      do alpha=0,nin(igrid)
         do i=0,13
            fevQCDb(i,alpha) = 0d0
            do beta=0,nin(igrid)
               do j=0,13
                  fevQCDb(i,alpha) = fevQCDb(i,alpha)
     1                 + EQCD(i,j,alpha,beta) * fevQCD(j,beta)
               enddo
            enddo
         enddo
      enddo
*
      do alpha=0,nin(igrid)
         do i=0,13
            fevQCD(i,alpha) = fevQCDb(i,alpha)
         enddo
      enddo
*
      return
      end
