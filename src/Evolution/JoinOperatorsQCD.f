************************************************************************
*
*     JoinOperatorsQCD.f:
*
*     This routine joins the QCD evolution operators computed with 
*     different numbers of active flavours.
*
*     QCD evolution basis:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     gm  Sg  g   V   V3  V8  V15 V24 V35 T3  T8  T15 T24 T35
*     
************************************************************************
      subroutine JoinOperatorsQCD(jgrid)
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
      double precision EvQCDb(0:13,0:13,0:nint_max,0:nint_max)
      double precision EvQCD(0:13,0:13,0:nint_max,0:nint_max)
*
      if(sgn.ne.1)then
         write(6,*) "In JoinOperatorsQCD.f:"
         write(6,*) "Backward evolution not allowed yet."
         call exit(-10)
      endif
*
      do alpha=0,nin(jgrid)
         do beta=alpha,nin(jgrid)
*     Set evolution operators to zero
            do i=0,13
               do j=0,13
                  EvQCD(i,j,alpha,beta) = 0d0
               enddo
            enddo
*     Singlet and Gluon
            do i=1,2
               do j=1,2
                  EvQCD(i,j,alpha,beta) =
     1                 MQCDsg(nfi,i,j,alpha,beta)
               enddo
            enddo
*     Total Valence
            EvQCD(3,3,alpha,beta) =
     1           MQCDnsv(nfi,alpha,beta)
*     V3
            EvQCD(4,4,alpha,beta) =
     1           MQCDnsm(nfi,alpha,beta)
*     V8
            EvQCD(5,5,alpha,beta) =
     1           MQCDnsm(nfi,alpha,beta)
*     V15
            if(nfi.lt.4)then
               EvQCD(6,3,alpha,beta) =
     1              MQCDnsv(nfi,alpha,beta)
            else
               EvQCD(6,6,alpha,beta) =
     1              MQCDnsm(nfi,alpha,beta)
            endif
*     V24
            if(nfi.lt.5)then
               EvQCD(7,3,alpha,beta) =
     1              MQCDnsv(nfi,alpha,beta)
            else
               EvQCD(7,7,alpha,beta) =
     1              MQCDnsm(nfi,alpha,beta)
            endif
*     V35
            if(nfi.lt.6)then
               EvQCD(8,3,alpha,beta) =
     1              MQCDnsv(nfi,alpha,beta)
            else
               EvQCD(8,8,alpha,beta) =
     1              MQCDnsm(nfi,alpha,beta)
            endif
*     T3
            EvQCD(9,9,alpha,beta) =
     1           MQCDnsp(nfi,alpha,beta)
*     T8
            EvQCD(10,10,alpha,beta) =
     1           MQCDnsp(nfi,alpha,beta)
*     T15
            if(nfi.lt.4)then
               do j=1,2
                  EvQCD(11,j,alpha,beta) =
     1                 MQCDsg(nfi,1,j,alpha,beta)
               enddo
            else
               EvQCD(11,11,alpha,beta) =
     1              MQCDnsp(nfi,alpha,beta)
            endif
*     T24
            if(nfi.lt.5)then
               do j=1,2
                  EvQCD(12,j,alpha,beta) =
     1                 MQCDsg(nfi,1,j,alpha,beta)
               enddo
            else
               EvQCD(12,12,alpha,beta) =
     1              MQCDnsp(nfi,alpha,beta)
            endif
*     T35
            if(nfi.lt.6)then
               do j=1,2
                  EvQCD(13,j,alpha,beta) =
     1                 MQCDsg(nfi,1,j,alpha,beta)
               enddo
            else
               EvQCD(13,13,alpha,beta) =
     1              MQCDnsp(nfi,alpha,beta)
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
            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
                  do i=0,13
                     do j=0,13
                        EvQCDb(i,j,alpha,beta) = 0d0
                     enddo
                  enddo
               enddo
            enddo
*     
            nfm = nf + 1
*     Get alphas value at the heavy quark threshold (with nfm active flavours)
            coup = asthUp(nfm)
c            coup = asthDown(nfm)
*     Contruct matching conditions at this threshod
            if(IsExt(jgrid))then
               do alpha=0,nin(jgrid)
                  do beta=alpha,nin(jgrid)
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
               do alpha=0,nin(jgrid)
                  do beta=alpha,nin(jgrid)
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
            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
*     Singlet and Gluon
                  do i=1,2
                     do j=1,2
                        do gamma=0,nin(jgrid)
                           do delta=0,nin(jgrid)
                              do k=1,2
                                 do l=1,2
                                    EvQCDb(i,j,alpha,beta) = 
     1                                   EvQCDb(i,j,alpha,beta)
     2                                   + MQCDsg(nf+1,i,k,alpha,gamma)
     3                                   * MatQCDsg(k,l,gamma,delta)
     4                                   * EvQCD(l,j,delta,beta)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
*     
                  do gamma=0,nin(jgrid)
                     do delta=0,nin(jgrid)
*     Total Valence
                        EvQCDb(3,3,alpha,beta) = 
     1                       EvQCDb(3,3,alpha,beta)
     2                       + MQCDnsv(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvQCD(3,3,delta,beta)
*     V3
                        EvQCDb(4,4,alpha,beta) = 
     1                       EvQCDb(4,4,alpha,beta)
     2                       + MQCDnsm(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvQCD(4,4,delta,beta)
*     V8
                        EvQCDb(5,5,alpha,beta) = 
     1                       EvQCDb(5,5,alpha,beta)
     2                       + MQCDnsm(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvQCD(5,5,delta,beta)
*     T3
                        EvQCDb(9,9,alpha,beta) = 
     1                       EvQCDb(9,9,alpha,beta)
     2                       + MQCDnsp(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvQCD(9,9,delta,beta)
*     T8
                        EvQCDb(10,10,alpha,beta) = 
     1                       EvQCDb(10,10,alpha,beta)
     2                       + MQCDnsp(nf+1,alpha,gamma)
     3                       * MatQCDns(gamma,delta)
     4                       * EvQCD(10,10,delta,beta)
                     enddo
                  enddo
*     Charm threshold
                  if(nfm.eq.4)then
                     do gamma=0,nin(jgrid)
                        do delta=0,nin(jgrid)
*     V15
                           EvQCDb(6,3,alpha,beta) = 
     1                          EvQCDb(6,3,alpha,beta)
     2                          + MQCDnsm(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvQCD(6,3,delta,beta)
*     V24
                           EvQCDb(7,3,alpha,beta) = 
     1                          EvQCDb(7,3,alpha,beta)
     2                          + MQCDnsv(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvQCD(7,3,delta,beta)
*     V35
                           EvQCDb(8,3,alpha,beta) = 
     1                          EvQCDb(8,3,alpha,beta)
     2                          + MQCDnsv(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvQCD(8,3,delta,beta)
*     T15
                           Match(1) = MatQCDns(gamma,delta) 
     1                          - 3d0 * ( MatQCDsg(1,1,gamma,delta) 
     2                          - MatQCDns(gamma,delta) )
                           Match(2) = - 3d0 * MatQCDsg(1,2,gamma,delta)
                           do j=1,2
                              do k=1,2
                                 EvQCDb(11,j,alpha,beta) = 
     1                                EvQCDb(11,j,alpha,beta)
     2                                + MQCDnsp(nf+sgn,alpha,gamma)
     3                                * Match(k)
     4                                * EvQCD(k,j,delta,beta)
                              enddo
                           enddo
*     
                           do j=1,2
*     T24
                              do k=1,2
                                 do l=1,2
                                    EvQCDb(12,j,alpha,beta) = 
     1                                   EvQCDb(12,j,alpha,beta)
     2                                   + MQCDsg(nf+1,1,k,alpha,gamma)
     3                                   * MatQCDsg(k,l,gamma,delta)
     4                                   * EvQCD(l,j,delta,beta)
                                 enddo
                              enddo
*     T35
                              EvQCDb(13,j,alpha,beta) = 
     1                             EvQCDb(12,j,alpha,beta)
                           enddo
                        enddo
                     enddo
*     Bottom threshold
                  elseif(nfm.eq.5)then
                     do gamma=0,nin(jgrid)
                        do delta=0,nin(jgrid)
*     V15
                           if(nfi.ge.4)then
                              EvQCDb(6,6,alpha,beta) = 
     1                             EvQCDb(6,6,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(6,6,delta,beta)
                           else
                              EvQCDb(6,3,alpha,beta) = 
     1                             EvQCDb(6,3,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(6,3,delta,beta)
                           endif
*     V24
                           EvQCDb(7,3,alpha,beta) = 
     1                          EvQCDb(7,3,alpha,beta)
     2                          + MQCDnsm(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvQCD(7,3,delta,beta)
*     V35
                           EvQCDb(8,3,alpha,beta) = 
     1                          EvQCDb(8,3,alpha,beta)
     2                          + MQCDnsv(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvQCD(8,3,delta,beta)
*     T15
                           if(nfi.ge.4)then
                              EvQCDb(11,11,alpha,beta) = 
     1                             EvQCDb(11,11,alpha,beta)
     2                             + MQCDnsp(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(11,11,delta,beta)
                           else
                              do j=1,2
                                 EvQCDb(11,j,alpha,beta) = 
     1                                EvQCDb(11,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * MatQCDns(gamma,delta)
     4                                * EvQCD(11,j,delta,beta)
                              enddo
                           endif
*     T24
                           Match(1) = MatQCDns(gamma,delta) 
     1                          - 4d0 * ( MatQCDsg(1,1,gamma,delta) 
     2                          - MatQCDns(gamma,delta) )
                           Match(2) = - 4d0 * MatQCDsg(1,2,gamma,delta)
                           do j=1,2
                              do k=1,2
                                 EvQCDb(12,j,alpha,beta) = 
     1                                EvQCDb(12,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * Match(k)
     4                                * EvQCD(k,j,delta,beta)
                              enddo
                           enddo
*     T35
                           do j=1,2
                              do k=1,2
                                 do l=1,2
                                    EvQCDb(13,j,alpha,beta) = 
     1                                   EvQCDb(13,j,alpha,beta)
     2                                   + MQCDsg(nf+1,1,k,alpha,gamma)
     3                                   * MatQCDsg(k,l,gamma,delta)
     4                                   * EvQCD(l,j,delta,beta)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
*     Top threshold
                  elseif(nfm.eq.6)then
                     do gamma=0,nin(jgrid)
                        do delta=0,nin(jgrid)
*     V15
                           if(nfi.ge.4)then
                              EvQCDb(6,6,alpha,beta) = 
     1                             EvQCDb(6,6,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(6,6,delta,beta)
                           else
                              EvQCDb(6,3,alpha,beta) = 
     1                             EvQCDb(6,3,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(6,3,delta,beta)
                           endif
*     V24
                           if(nfi.ge.5)then
                              EvQCDb(7,7,alpha,beta) = 
     1                             EvQCDb(7,7,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(7,7,delta,beta)
                           else
                              EvQCDb(7,3,alpha,beta) = 
     1                             EvQCDb(7,3,alpha,beta)
     2                             + MQCDnsm(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(7,3,delta,beta)
                           endif
*     V35
                           EvQCDb(8,3,alpha,beta) = 
     1                          EvQCDb(8,3,alpha,beta)
     2                          + MQCDnsm(nf+1,alpha,gamma)
     3                          * MatQCDns(gamma,delta)
     4                          * EvQCD(8,3,delta,beta)
*     T15
                           if(nfi.ge.4)then
                              EvQCDb(11,11,alpha,beta) = 
     1                             EvQCDb(11,11,alpha,beta)
     2                             + MQCDnsp(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(11,11,delta,beta)
                           else
                              do j=1,2
                                 EvQCDb(11,j,alpha,beta) = 
     1                                EvQCDb(11,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * MatQCDns(gamma,delta)
     4                                * EvQCD(11,j,delta,beta)
                              enddo
                           endif
*     T24
                           if(nfi.ge.5)then
                              EvQCDb(12,12,alpha,beta) = 
     1                             EvQCDb(12,12,alpha,beta)
     2                             + MQCDnsp(nf+1,alpha,gamma)
     3                             * MatQCDns(gamma,delta)
     4                             * EvQCD(12,12,delta,beta)
                           else
                              do j=1,2
                                 EvQCDb(12,j,alpha,beta) = 
     1                                EvQCDb(12,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * MatQCDns(gamma,delta)
     4                                * EvQCD(12,j,delta,beta)
                              enddo
                           endif
*     T35
                           Match(1) = MatQCDns(gamma,delta) 
     1                          - 5d0 * ( MatQCDsg(1,1,gamma,delta) 
     2                          - MatQCDns(gamma,delta) )
                           Match(2) = - 5d0 * MatQCDsg(1,2,gamma,delta)
                           do j=1,2
                              do k=1,2
                                 EvQCDb(13,j,alpha,beta) = 
     1                                EvQCDb(13,j,alpha,beta)
     2                                + MQCDnsp(nf+1,alpha,gamma)
     3                                * Match(k)
     4                                * EvQCD(k,j,delta,beta)
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
            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
                  do i=0,13
                     do j=0,13
                        EvQCD(i,j,alpha,beta) = EvQCDb(i,j,alpha,beta)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
*
*     Tranform the evolution operators from the evolution to the physical basis
*
      do alpha=0,nin(jgrid)
         do beta=0,nin(jgrid)
            do i=0,13
               do j=0,13
                  Ev2EvQCD(jgrid,i,j,alpha,beta) = EvQCD(i,j,alpha,beta)
                  Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 0d0
                  Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 0d0
                  do k=0,13
                     Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 
     1                    Ev2PhQCD(jgrid,i-7,j,alpha,beta)
     2                       + Tev2phQCD(nff,i,k)
     3                       * EvQCD(k,j,alpha,beta)
                     do l=0,13
                        Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 
     1                       Ph2PhQCD(jgrid,i-7,j-7,alpha,beta)
     2                       + Tev2phQCD(nff,i,k)
     3                       * EvQCD(k,l,alpha,beta)
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
