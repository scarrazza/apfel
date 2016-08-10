************************************************************************
*
*     JoinOperatorsUni.f:
*
*     This routine joins the Unified evolution operators computed with 
*     different numbers of active flavours.
*
*     Unified evolution basis:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     g   gm  Sig Dsg Tu1 Tu2 Td1 Td2 V   DV  Vu1 Vu2 Vd1 Vd2
*     
************************************************************************
      subroutine JoinOperatorsUni(jgrid)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/ThresholdAlphaQCD.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/transQCD.h"
      include "../commons/EvolutionOperator.h"
**
*     Input Variables
*
      integer jgrid
**
*     Internal Variables
*
      integer i,j,k,l
      integer nf,nfm
      integer alpha,beta,gamma,delta
      double precision coup,integralsMatching
      double precision MatQCDns(0:nint_max,0:nint_max),Match(2)
      double precision MatQCDsg(2,2,0:nint_max,0:nint_max)
      double precision EvUnib(0:13,0:13,0:nint_max,0:nint_max)
      double precision EvUni(0:13,0:13,0:nint_max,0:nint_max)
*
      if(sgn.ne.1)then
         write(6,*) "In JoinOperatorsUni.f:"
         write(6,*) "Backward evolution not allowed yet."
         call exit(-10)
      endif
*
      do alpha=0,nin(igrid)
         do beta=alpha,nin(igrid)
*     Set evolution operators to zero
            do i=0,13
               do j=0,13
                  EvUni(i,j,alpha,beta) = 0d0
               enddo
            enddo
*     Gluon, photon, Sigma, DeltaSigma
            do i=0,3
               do j=0,3
                  EvUni(i,j,alpha,beta) =
     1                 MUnisg1(nfi,nli,i,j,alpha,beta)
               enddo
            enddo
*     Tu1
            if(nfi.lt.4)then
               do j=0,3
                  EvUni(4,j,alpha,beta) =
     1                 ( MUnisg1(nfi,nli,2,j,alpha,beta) 
     2                 + MUnisg1(nfi,nli,3,j,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(4,4,alpha,beta) = MUninspu(nfi,nli,alpha,beta)
            endif
*     Tu2
            if(nfi.lt.6)then
               do j=0,3
                  EvUni(5,j,alpha,beta) =
     1                 ( MUnisg1(nfi,nli,2,j,alpha,beta) 
     2                 + MUnisg1(nfi,nli,3,j,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(5,5,alpha,beta) = MUninspu(nfi,nli,alpha,beta)
            endif
*     Td1
            EvUni(6,6,alpha,beta) = MUninspd(nfi,nli,alpha,beta)
*     Td2
            if(nfi.lt.5)then
               do j=0,3
                  EvUni(7,j,alpha,beta) =
     1                 ( MUnisg1(nfi,nli,2,j,alpha,beta) 
     2                 - MUnisg1(nfi,nli,3,j,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(7,7,alpha,beta) = MUninspd(nfi,nli,alpha,beta)
            endif
*     Total Valence, DeltaValence
            do i=1,2
               do j=1,2
                  EvUni(7+i,7+j,alpha,beta) =
     1                 MUnisg2(nfi,nli,i,j,alpha,beta)
               enddo
            enddo
*     Vu1
            if(nfi.lt.4)then
               do j=1,2
                  EvUni(10,j,alpha,beta) =
     1                 ( MUnisg2(nfi,nli,1,j,alpha,beta) 
     2                 + MUnisg2(nfi,nli,2,j,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(10,10,alpha,beta) = MUninsmu(nfi,nli,alpha,beta)
            endif
*     Vu2
            if(nfi.lt.6)then
               do j=1,2
                  EvUni(11,j,alpha,beta) =
     1                 ( MUnisg2(nfi,nli,1,j,alpha,beta) 
     2                 + MUnisg2(nfi,nli,2,j,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(11,11,alpha,beta) = MUninsmu(nfi,nli,alpha,beta)
            endif
*     Vd1
            EvUni(12,12,alpha,beta) = MUninsmd(nfi,nli,alpha,beta)
*     Vd2
            if(nfi.lt.5)then
               do j=1,2
                  EvUni(13,j,alpha,beta) =
     1                 ( MUnisg2(nfi,nli,1,j,alpha,beta) 
     2                 - MUnisg2(nfi,nli,2,j,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(13,13,alpha,beta) = MUninsmd(nfi,nli,alpha,beta)
            endif
         enddo
      enddo
*
*     If the initial and the final numbers of flavours are different ...
*
      if(nfi.ne.nff)then
         do nf=nfi,nff-1
*
*     Set temporary evolution operators to zero
*
            do alpha=0,nin(igrid)
               do beta=alpha,nin(igrid)
                  do i=0,13
                     do j=0,13
                        EvUnib(i,j,alpha,beta) = 0d0
                     enddo
                  enddo
               enddo
            enddo
*     
            nfm = nf + 1
*     Get alphas value at the heavy quark threshold (with nfm active flavours)
c            coup = asthUp(nfm)
            coup = asthDown(nfm)
*     Contruct matching conditions at this threshod in the QCD basis
            if(IsExt(igrid))then
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     MatQCDns(alpha,beta)     =
     1                    integralsMatching(nf+1,alpha,beta,coup,1,sgn)
                     MatQCDsg(1,1,alpha,beta) =
     1                    integralsMatching(nf+1,alpha,beta,coup,2,sgn)
                     MatQCDsg(1,2,alpha,beta) =
     1                    integralsMatching(nf+1,alpha,beta,coup,3,sgn)
                     MatQCDsg(2,1,alpha,beta) =
     1                    integralsMatching(nf+1,alpha,beta,coup,4,sgn)
                     MatQCDsg(2,2,alpha,beta) =
     1                    integralsMatching(nf+1,alpha,beta,coup,5,sgn)
                  enddo
               enddo
            else
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     MatQCDns(alpha,beta)     =
     1                   integralsMatching(nf+1,0,beta-alpha,coup,1,sgn)
                     MatQCDsg(1,1,alpha,beta) =
     1                   integralsMatching(nf+1,0,beta-alpha,coup,2,sgn)
                     MatQCDsg(1,2,alpha,beta) =
     1                   integralsMatching(nf+1,0,beta-alpha,coup,3,sgn)
                     MatQCDsg(2,1,alpha,beta) =
     1                   integralsMatching(nf+1,0,beta-alpha,coup,4,sgn)
                     MatQCDsg(2,2,alpha,beta) =
     1                   integralsMatching(nf+1,0,beta-alpha,coup,5,sgn)
                  enddo
               enddo
            endif
*     
*     Now combine evolution tables for different numbers of active
*     flavours.
*     
            do alpha=0,nin(igrid)
               do beta=alpha,nin(igrid)
*     Singlet and Gluon
                  do i=1,2
                     do j=1,2
                        do gamma=0,nin(igrid)
                           do delta=0,nin(igrid)
                              do k=1,2
                                 do l=1,2
                                    EvUnib(i,j,alpha,beta) = 
     1                                   EvUnib(i,j,alpha,beta)
     2                                   + MQCDsg(nf+1,i,k,alpha,gamma)
     3                                   * MatQCDsg(k,l,gamma,delta)
     4                                   * EvUni(l,j,delta,beta)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
*     
                  do gamma=0,nin(igrid)
                     do delta=0,nin(igrid)
*     Total Valence
                        EvUnib(3,3,alpha,beta) = 
     1                       EvUnib(3,3,alpha,beta)
     2                       + MQCDnsv(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvUni(3,3,delta,beta)
*     V3
                        EvUnib(4,4,alpha,beta) = 
     1                       EvUnib(4,4,alpha,beta)
     2                       + MQCDnsm(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvUni(4,4,delta,beta)
*     V8
                        EvUnib(5,5,alpha,beta) = 
     1                       EvUnib(5,5,alpha,beta)
     2                       + MQCDnsm(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvUni(5,5,delta,beta)
*     T3
                        EvUnib(9,9,alpha,beta) = 
     1                       EvUnib(9,9,alpha,beta)
     2                       + MQCDnsp(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvUni(9,9,delta,beta)
*     T8
                        EvUnib(10,10,alpha,beta) = 
     1                       EvUnib(10,10,alpha,beta)
     2                       + MQCDnsp(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvUni(10,10,delta,beta)
                     enddo
                  enddo
*     Charm threshold
                  if(nfm.eq.4)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
*     V15
                           EvUnib(6,3,alpha,beta) = 
     1                          EvUnib(6,3,alpha,beta)
     2                          + MQCDnsm(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvUni(6,3,delta,beta)
*     V24
                           EvUnib(7,3,alpha,beta) = 
     1                          EvUnib(7,3,alpha,beta)
     2                          + MQCDnsv(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvUni(7,3,delta,beta)
*     V35
                           EvUnib(8,3,alpha,beta) = 
     1                          EvUnib(8,3,alpha,beta)
     2                          + MQCDnsv(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvUni(8,3,delta,beta)
*     T15
                           Match(1) = MatQCDns(gamma,delta) 
     1                          - 3d0 * ( MatQCDsg(1,1,gamma,delta) 
     2                          - MatQCDns(gamma,delta) )
                           Match(2) = - 3d0 * MatQCDsg(1,2,gamma,delta)
                           do j=1,2
                              do k=1,2
                                 EvUnib(11,j,alpha,beta) = 
     1                                EvUnib(11,j,alpha,beta)
     2                                + MQCDnsp(nf+sgn,alpha,gamma)
     3                                * Match(k)
     4                                * EvUni(k,j,delta,beta)
                              enddo
                           enddo
*     
                           do j=1,2
*     T24
                              do k=1,2
                                 do l=1,2
                                    EvUnib(12,j,alpha,beta) = 
     1                                   EvUnib(12,j,alpha,beta)
     2                                   + MQCDsg(nf+1,1,k,alpha,gamma)
     3                                   * MatQCDsg(k,l,gamma,delta)
     4                                   * EvUni(l,j,delta,beta)
                                 enddo
                              enddo
*     T35
                              EvUnib(13,j,alpha,beta) = 
     1                             EvUnib(12,j,alpha,beta)
                           enddo
                        enddo
                     enddo
*     Bottom threshold
                  elseif(nfm.eq.5)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
*     V15
                           if(nfi.ge.4)then
                              EvUnib(6,6,alpha,beta) = 
     1                             EvUnib(6,6,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(6,6,delta,beta)
                           else
                              EvUnib(6,3,alpha,beta) = 
     1                             EvUnib(6,3,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(6,3,delta,beta)
                           endif
*     V24
                           EvUnib(7,3,alpha,beta) = 
     1                          EvUnib(7,3,alpha,beta)
     2                          + MQCDnsm(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvUni(7,3,delta,beta)
*     V35
                           EvUnib(8,3,alpha,beta) = 
     1                          EvUnib(8,3,alpha,beta)
     2                          + MQCDnsv(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvUni(8,3,delta,beta)
*     T15
                           if(nfi.ge.4)then
                              EvUnib(11,11,alpha,beta) = 
     1                             EvUnib(11,11,alpha,beta)
     2                             + MQCDnsp(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(11,11,delta,beta)
                           else
                              do j=1,2
                                 EvUnib(11,j,alpha,beta) = 
     1                                EvUnib(11,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * MatQCDns(gamma,delta)
     4                                * EvUni(11,j,delta,beta)
                              enddo
                           endif
*     T24
                           Match(1) = MatQCDns(gamma,delta) 
     1                          - 4d0 * ( MatQCDsg(1,1,gamma,delta) 
     2                          - MatQCDns(gamma,delta) )
                           Match(2) = - 4d0 * MatQCDsg(1,2,gamma,delta)
                           do j=1,2
                              do k=1,2
                                 EvUnib(12,j,alpha,beta) = 
     1                                EvUnib(12,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * Match(k)
     4                                * EvUni(k,j,delta,beta)
                              enddo
                           enddo
*     T35
                           do j=1,2
                              do k=1,2
                                 do l=1,2
                                    EvUnib(13,j,alpha,beta) = 
     1                                   EvUnib(13,j,alpha,beta)
     2                                   + MQCDsg(nf+1,1,k,alpha,gamma)
     3                                   * MatQCDsg(k,l,gamma,delta)
     4                                   * EvUni(l,j,delta,beta)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
*     Top threshold
                  elseif(nfm.eq.6)then
                     do gamma=0,nin(igrid)
                        do delta=0,nin(igrid)
*     V15
                           if(nfi.ge.4)then
                              EvUnib(6,6,alpha,beta) = 
     1                             EvUnib(6,6,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(6,6,delta,beta)
                           else
                              EvUnib(6,3,alpha,beta) = 
     1                             EvUnib(6,3,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(6,3,delta,beta)
                           endif
*     V24
                           if(nfi.ge.5)then
                              EvUnib(7,7,alpha,beta) = 
     1                             EvUnib(7,7,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(7,7,delta,beta)
                           else
                              EvUnib(7,3,alpha,beta) = 
     1                             EvUnib(7,3,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(7,3,delta,beta)
                           endif
*     V35
                           EvUnib(8,3,alpha,beta) = 
     1                          EvUnib(8,3,alpha,beta)
     2                          + MQCDnsm(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvUni(8,3,delta,beta)
*     T15
                           if(nfi.ge.4)then
                              EvUnib(11,11,alpha,beta) = 
     1                             EvUnib(11,11,alpha,beta)
     2                             + MQCDnsp(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(11,11,delta,beta)
                           else
                              do j=1,2
                                 EvUnib(11,j,alpha,beta) = 
     1                                EvUnib(11,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * MatQCDns(gamma,delta)
     4                                * EvUni(11,j,delta,beta)
                              enddo
                           endif
*     T24
                           if(nfi.ge.5)then
                              EvUnib(12,12,alpha,beta) = 
     1                             EvUnib(12,12,alpha,beta)
     2                             + MQCDnsp(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvUni(12,12,delta,beta)
                           else
                              do j=1,2
                                 EvUnib(12,j,alpha,beta) = 
     1                                EvUnib(12,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * MatQCDns(gamma,delta)
     4                                * EvUni(12,j,delta,beta)
                              enddo
                           endif
*     T35
                           Match(1) = MatQCDns(gamma,delta) 
     1                          - 5d0 * ( MatQCDsg(1,1,gamma,delta) 
     2                          - MatQCDns(gamma,delta) )
                           Match(2) = - 5d0 * MatQCDsg(1,2,gamma,delta)
                           do j=1,2
                              do k=1,2
                                 EvUnib(13,j,alpha,beta) = 
     1                                EvUnib(13,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * Match(k)
     4                                * EvUni(k,j,delta,beta)
                              enddo
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
*     
*     Copy the backup evolution operators into the main ones
*     
            do alpha=0,nin(igrid)
               do beta=alpha,nin(igrid)
                  do i=0,13
                     do j=0,13
                        EvUni(i,j,alpha,beta) = EvUnib(i,j,alpha,beta)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
*
*     Tranform the evolution operators from the evolution to the physical basis
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            do i=0,13
               do j=0,13
                  Ev2EvUni(jgrid,i,j,alpha,beta) = EvUni(i,j,alpha,beta)
                  Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 0d0
                  Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 0d0
                  do k=0,13
                     Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 
     1                    Ev2PhQCD(jgrid,i-7,j,alpha,beta)
     2                       + Tev2phQCD(nff,i,k)
     3                       * EvUni(k,j,alpha,beta)
                     do l=0,13
                        Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 
     1                       Ph2PhQCD(jgrid,i-7,j-7,alpha,beta)
     2                       + Tev2phQCD(nff,i,k)
     3                       * EvUni(K,alpha,beta)
     4                       * Tph2evQCD(nfi,l,j)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c*
c*     Print rotation matrices
c*
c      do k=3,6
c         write(56,*) "Tph2evQCD, nf =",k
c         do i=0,13
c            write(56,"(a,14(2x,i3),a)")
c     1           "{",(int(Tph2evQCD(k,i,j)),j=0,13)," }"
c         enddo
c         write(56,*) "   "
c         write(57,*) "120d0*Tev2phQCD, nf =",k
c         do i=0,13
c            write(57,"(a,14(2x,i4),a)") 
c     1           "{",(int(120d0*Tev2phQCD(k,i,j)),j=0,13)," }"
c         enddo
c         write(57,*) "   "
c      enddo
c      stop
c$$$*
c$$$*     Rotation matrix factory
c$$$*
c$$$      do j=0,13
c$$$         do i=0,13
c$$$            Tev2phQCD(5,i,j) = Tev2phQCD(6,i,j)
c$$$         enddo
c$$$      enddo
c$$$      do j=0,13
c$$$         Tev2phQCD(5,j,1)  = Tev2phQCD(6,j,1) + Tev2phQCD(6,j,13)
c$$$         Tev2phQCD(5,j,13) = 0d0
c$$$         Tev2phQCD(5,j,3) = Tev2phQCD(6,j,3) + Tev2phQCD(6,j,8)
c$$$         Tev2phQCD(5,j,8) = 0d0
c$$$      enddo
c$$$      do j=0,13
c$$$         do i=0,13
c$$$            Tev2phQCD(4,i,j) = Tev2phQCD(5,i,j)
c$$$         enddo
c$$$      enddo
c$$$      do j=0,13
c$$$         Tev2phQCD(4,j,1)  = Tev2phQCD(5,j,1) + Tev2phQCD(5,j,12)
c$$$         Tev2phQCD(4,j,12) = 0d0
c$$$
c$$$         Tev2phQCD(4,j,3) = Tev2phQCD(5,j,3) + Tev2phQCD(5,j,7)
c$$$         Tev2phQCD(4,j,7) = 0d0
c$$$      enddo
c$$$      do j=0,13
c$$$         do i=0,13
c$$$            Tev2phQCD(3,i,j) = Tev2phQCD(4,i,j)
c$$$         enddo
c$$$      enddo
c$$$      do j=0,13
c$$$         Tev2phQCD(3,j,1)  = Tev2phQCD(4,j,1) + Tev2phQCD(4,j,11)
c$$$         Tev2phQCD(3,j,11) = 0d0
c$$$         Tev2phQCD(3,j,3) = Tev2phQCD(4,j,3) + Tev2phQCD(4,j,6)
c$$$         Tev2phQCD(3,j,6) = 0d0
c$$$      enddo
c$$$      do k=3,5
c$$$         do j=0,13
c$$$            do i=0,13
c$$$               Tph2evQCD(k,i,j) = 0d0
c$$$            enddo
c$$$         enddo
c$$$         Tph2evQCD(k,0,0) = 1d0
c$$$      enddo
c$$$      do j=2,12
c$$$         do i=0,13
c$$$            Tph2evQCD(5,i,j) = Tph2evQCD(6,i,j)
c$$$         enddo
c$$$      enddo
c$$$      do j=3,11
c$$$         do i=0,13
c$$$            Tph2evQCD(4,i,j) = Tph2evQCD(5,i,j)
c$$$         enddo
c$$$      enddo
c$$$      do j=4,10
c$$$         do i=0,13
c$$$            Tph2evQCD(3,i,j) = Tph2evQCD(4,i,j)
c$$$         enddo
c$$$      enddo
c$$$      do k=3,6
c$$$         do i=0,13
c$$$            do j=0,13
c$$$               if(dabs(Tev2phQCD(k,i,j)).lt.1d-8) Tev2phQCD(k,i,j) = 0d0 
c$$$               write(56,"(a,i1,a,i2,a,i2,a,f18.14,a)")
c$$$     2              "      data Tev2phQCD(",k,",",i,",",j,
c$$$     1              ")   /",Tev2phQCD(k,i,j),"d0 /"
c$$$            enddo
c$$$         enddo
c$$$         write(56,*) "  "
c$$$      enddo
c$$$      do k=3,6
c$$$         do i=0,13
c$$$            do j=0,13
c$$$               if(dabs(Tph2evQCD(k,i,j)).lt.1d-8) Tph2evQCD(k,i,j) = 0d0 
c$$$               write(56,"(a,i1,a,i2,a,i2,a,f18.14,a)")
c$$$     2              "      data Tph2evQCD(",k,",",i,",",j,
c$$$     1              ")   /",Tph2evQCD(k,i,j),"d0 /"
c$$$            enddo
c$$$         enddo
c$$$         write(56,*) "  "
c$$$      enddo
c$$$      stop
*
      return
      end
