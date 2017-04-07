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
      include "../commons/transUni.h"
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
      integer nf,nfm,mfi,mff,nl
      integer sgnnf(3:5)
      integer alpha,beta,gamma
      integer a,b
      double precision coup,integralsMatching
      double precision MatchUni(0:13,0:13,0:nint_max,0:nint_max)
      double precision Mns,Msg(2,2),MatchQCD(0:13,0:13)
      double precision EvUni(0:13,0:13,0:nint_max,0:nint_max)
      double precision EvUniM(0:13,0:13,0:nint_max,0:nint_max)
      double precision EvUniN(0:13,0:13,0:nint_max,0:nint_max)
      double precision deltaab
*
*     No backward evolution allowed (yet)
*
      if(sgn.ne.1)then
         write(6,*) "In JoinOperatorsUni.f:"
         write(6,*) "Backward evolution not allowed yet."
         call exit(-10)
      endif
*
*     Check that the tau threshold is never crossed
*
      if(nli.ne.nlf)then
         write(6,*) "In JoinOperatorsUni.f:"
         write(6,*) "Tau mass threshold crossing not allowed."
         write(6,*) "Set the tau mass either to zero or to infinity."
         call exit(-10)
      endif
*
      nl = nli
*
      mfi = nfli(nl)
      mff = nflf(nl)
      do alpha=0,nin(jgrid)
         do beta=alpha,nin(jgrid)
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
     1                 MUnisg1(mfi,nl,i+1,j+1,alpha,beta)
               enddo
            enddo
*     Tu1
            if(mfi.lt.4)then
               do j=0,3
                  EvUni(4,j,alpha,beta) =
     1                 ( MUnisg1(mfi,nl,3,j+1,alpha,beta) 
     2                 + MUnisg1(mfi,nl,4,j+1,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(4,4,alpha,beta) = MUninspu(mfi,nl,alpha,beta)
            endif
*     Tu2
            if(mfi.lt.6)then
               do j=0,3
                  EvUni(5,j,alpha,beta) =
     1                 ( MUnisg1(mfi,nl,3,j+1,alpha,beta) 
     2                 + MUnisg1(mfi,nl,4,j+1,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(5,5,alpha,beta) = MUninspu(mfi,nl,alpha,beta)
            endif
*     Td1
            EvUni(6,6,alpha,beta) = MUninspd(mfi,nl,alpha,beta)
*     Td2
            if(mfi.lt.5)then
               do j=0,3
                  EvUni(7,j,alpha,beta) =
     1                 ( MUnisg1(mfi,nl,3,j+1,alpha,beta) 
     2                 - MUnisg1(mfi,nl,4,j+1,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(7,7,alpha,beta) = MUninspd(mfi,nl,alpha,beta)
            endif
*     Total Valence, DeltaValence
            do i=8,9
               do j=8,9
                  EvUni(i,j,alpha,beta) =
     1                 MUnisg2(mfi,nl,i-7,j-7,alpha,beta)
               enddo
            enddo
*     Vu1
            if(mfi.lt.4)then
               do j=8,9
                  EvUni(10,j,alpha,beta) =
     1                 ( MUnisg2(mfi,nl,1,j-7,alpha,beta) 
     2                 + MUnisg2(mfi,nl,2,j-7,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(10,10,alpha,beta) = MUninsmu(mfi,nl,alpha,beta)
            endif
*     Vu2
            if(mfi.lt.6)then
               do j=8,9
                  EvUni(11,j,alpha,beta) =
     1                 ( MUnisg2(mfi,nl,1,j-7,alpha,beta) 
     2                 + MUnisg2(mfi,nl,2,j-7,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(11,11,alpha,beta) = MUninsmu(mfi,nl,alpha,beta)
            endif
*     Vd1
            EvUni(12,12,alpha,beta) = MUninsmd(mfi,nl,alpha,beta)
*     Vd2
            if(mfi.lt.5)then
               do j=8,9
                  EvUni(13,j,alpha,beta) =
     1                 ( MUnisg2(mfi,nl,1,j-7,alpha,beta) 
     2                 - MUnisg2(mfi,nl,2,j-7,alpha,beta) ) / 2d0
               enddo
            else
               EvUni(13,13,alpha,beta) = MUninsmd(mfi,nl,alpha,beta)
            endif
         enddo
      enddo
*
*     If the initial and the final numbers of flavours are different ...
*
      if(mfi.ne.mff)then
*     Sign to be applied to the matching conditions of DeltaSigma
         sgnnf(3) =  1
         sgnnf(4) = -1
         sgnnf(5) =  1
         do nf=mfi,mff-1
            nfm = nf + 1
*     Get alphas value at the heavy quark threshold (with nfm active flavours)
            coup = asthUp(nfm)
c            coup = asthDown(nfm)
*     Contruct matching conditions at this threshod in the QCD basis
            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
                  if(IsExt(jgrid))then
                     a = alpha
                     b = beta
                  else
                     a = 0
                     b = beta - alpha
                  endif
*     Define delta function
                  deltaab = 0d0
                  if(beta.eq.alpha) deltaab = 1d0
*     Call QCD matching conditions
                  Mns      = integralsMatching(nfm,a,b,coup,1,sgn)
                  Msg(1,1) = integralsMatching(nfm,a,b,coup,2,sgn)
                  Msg(1,2) = integralsMatching(nfm,a,b,coup,3,sgn)
                  Msg(2,1) = integralsMatching(nfm,a,b,coup,4,sgn)
                  Msg(2,2) = integralsMatching(nfm,a,b,coup,5,sgn)
*
*     Matching condition in the QCD evolution basis
*
                  do i=0,13
                     do j=0,13
                        MatchQCD(i,j) = 0d0
                     enddo
                  enddo
*     Gamma
                  MatchQCD(0,0) = deltaab
*     Singlet and gluon
                  do i=1,2
                     do j=1,2
                        MatchQCD(i,j) = Msg(i,j)
                     enddo
                  enddo
*     Total Valence
                  MatchQCD(3,3)   = Mns
*     V3
                  MatchQCD(4,4)   = Mns
*     V8
                  MatchQCD(5,5)   = Mns
*     T3
                  MatchQCD(9,9)   = Mns
*     T8
                  MatchQCD(10,10) = Mns
*     Charm threshold
                  if(nfm.eq.4)then
*     V15
                     MatchQCD(6,3) = Mns
*     V24
                     MatchQCD(7,3) = Mns
*     V35
                     MatchQCD(8,3) = Mns
*     T15
                     MatchQCD(11,1) = Mns - 3d0 * ( Msg(1,1) - Mns )
                     MatchQCD(11,2) = - 3d0 * Msg(1,2)
*     T24
                     MatchQCD(12,1) = Msg(1,1)
                     MatchQCD(12,2) = Msg(1,2)
*     T35
                     MatchQCD(13,1) = Msg(1,1)
                     MatchQCD(13,2) = Msg(1,2)
*     Bottom threshold
                  elseif(nfm.eq.5)then
*     V15
                     MatchQCD(6,6) = Mns
*     V24
                     MatchQCD(7,3) = Mns
*     V35
                     MatchQCD(8,3) = Mns
*     T15
                     MatchQCD(11,11) = Mns
*     T24
                     MatchQCD(12,1) = Mns - 4d0 * ( Msg(1,1) - Mns )
                     MatchQCD(12,2) = - 4d0 * Msg(1,2)
*     T35
                     MatchQCD(13,1) = Msg(1,1)
                     MatchQCD(13,2) = Msg(1,2)
*     Top threshold
                  elseif(nfm.eq.6)then
*     V15
                     MatchQCD(6,6) = Mns
*     V24
                     MatchQCD(7,7) = Mns
*     V35
                     MatchQCD(8,3) = Mns
*     T15
                     MatchQCD(11,11) = Mns
*     T24
                     MatchQCD(12,12) = Mns
*     T35
                     MatchQCD(13,1) = Mns - 5d0 * ( Msg(1,1) - Mns )
                     MatchQCD(13,2) = - 5d0 * Msg(1,2)
                  endif
*
*     Rotate mathcing conditions from QCD to the Uni basis
*
                  do i=0,13
                     do j=0,13
                        MatchUni(i,j,alpha,beta) = 0d0
                        do k=0,13
                           do l=0,13
                              MatchUni(i,j,alpha,beta) =
     1                             MatchUni(i,j,alpha,beta)
     2                             + TevQCD2evUni(i,k)
     3                             * MatchQCD(k,l)
     4                             * TevUni2evQCD(l,j)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
*
*     Apply matching conditions
*
            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
                  do i=0,13
                     do j=0,13
                        EvUniM(i,j,alpha,beta) = 0d0
                        do gamma=alpha,beta
                           do k=0,13
                              EvUniM(i,j,alpha,beta) = 
     1                             EvUniM(i,j,alpha,beta)
     2                             + MatchUni(i,k,alpha,gamma)
     3                             * EvUni(k,j,gamma,beta)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
*
*     Now construct evolution matrix for the next step
*

            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
*     Set evolution operators to zero
                  do i=0,13
                     do j=0,13
                        EvUniN(i,j,alpha,beta) = 0d0
                     enddo
                  enddo
*     Gluon, photon, Sigma, DeltaSigma
                  do i=0,3
                     do j=0,3
                        EvUniN(i,j,alpha,beta) =
     1                       MUnisg1(nfm,nl,i+1,j+1,alpha,beta)
                     enddo
                  enddo
*     Tu1
                  if(nfm.lt.4)then
                     do j=0,3
                        EvUniN(4,j,alpha,beta) =
     1                       ( MUnisg1(nfm,nl,3,j+1,alpha,beta) 
     2                       + MUnisg1(nfm,nl,4,j+1,alpha,beta) ) / 2d0
                     enddo
                  else
                     EvUniN(4,4,alpha,beta) =
     1                    MUninspu(nfm,nl,alpha,beta)
                  endif
*     Tu2
                  if(nfm.lt.6)then
                     do j=0,3
                        EvUniN(5,j,alpha,beta) =
     1                       ( MUnisg1(nfm,nl,3,j+1,alpha,beta) 
     2                       + MUnisg1(nfm,nl,4,j+1,alpha,beta) ) / 2d0
                     enddo
                  else
                     EvUniN(5,5,alpha,beta) =
     1                    MUninspu(nfm,nl,alpha,beta)
                  endif
*     Td1
                  EvUniN(6,6,alpha,beta) = MUninspd(nfm,nl,alpha,beta)
*     Td2
                  if(nfm.lt.5)then
                     do j=0,3
                        EvUniN(7,j,alpha,beta) =
     1                       ( MUnisg1(nfm,nl,3,j+1,alpha,beta) 
     2                       - MUnisg1(nfm,nl,4,j+1,alpha,beta) ) / 2d0
                     enddo
                  else
                     EvUniN(7,7,alpha,beta) =
     1                    MUninspd(nfm,nl,alpha,beta)
                  endif
*     Total Valence, DeltaValence
                  do i=8,9
                     do j=8,9
                        EvUniN(i,j,alpha,beta) =
     1                       MUnisg2(nfm,nl,i-7,j-7,alpha,beta)
                     enddo
                  enddo
*     Vu1
                  if(nfm.lt.4)then
                     do j=8,9
                        EvUniN(10,j,alpha,beta) =
     1                       ( MUnisg2(nfm,nl,1,j-7,alpha,beta) 
     2                       + MUnisg2(nfm,nl,2,j-7,alpha,beta) ) / 2d0
                     enddo
                  else
                     EvUniN(10,10,alpha,beta) =
     1                    MUninsmu(nfm,nl,alpha,beta)
                  endif
*     Vu2
                  if(nfm.lt.6)then
                     do j=8,9
                        EvUniN(11,j,alpha,beta) =
     1                       ( MUnisg2(nfm,nl,1,j-7,alpha,beta) 
     2                       + MUnisg2(nfm,nl,2,j-7,alpha,beta) ) / 2d0
                     enddo
                  else
                     EvUniN(11,11,alpha,beta) =
     1                    MUninsmu(nfm,nl,alpha,beta)
                  endif
*     Vd1
                  EvUniN(12,12,alpha,beta) = MUninsmd(nfm,nl,alpha,beta)
*     Vd2
                  if(nfm.lt.5)then
                     do j=8,9
                        EvUniN(13,j,alpha,beta) =
     1                       ( MUnisg2(nfm,nl,1,j-7,alpha,beta) 
     2                       - MUnisg2(nfm,nl,2,j-7,alpha,beta) ) / 2d0
                     enddo
                  else
                     EvUniN(13,13,alpha,beta) =
     1                    MUninsmd(nfm,nl,alpha,beta)
                  endif
               enddo
            enddo
*
*     Apply next step
*
            do alpha=0,nin(jgrid)
               do beta=alpha,nin(jgrid)
                  do i=0,13
                     do j=0,13
                        EvUni(i,j,alpha,beta) = 0d0
                        do gamma=alpha,beta
                           do k=0,13
                              EvUni(i,j,alpha,beta) = 
     1                             EvUni(i,j,alpha,beta)
     3                             + EvUniN(i,k,alpha,gamma)
     4                             * EvUniM(k,j,gamma,beta)
                           enddo
                        enddo
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
*     Transform the evolution operator from the Uni to the QCD evolution basis
                  Ev2EvQCD(jgrid,i,j,alpha,beta) = 0d0
                  do k=0,13
                     do l=0,13
                        Ev2EvQCD(jgrid,i,j,alpha,beta) =
     1                       Ev2EvQCD(jgrid,i,j,alpha,beta)
     2                       + TevUni2evQCD(i,k)
     3                       * EvUni(k,l,alpha,beta)
     4                       * TevQCD2evUni(l,j)
                     enddo
                  enddo
*
*     Adjust the evolution operator to reabsorb distribution that reduce
*     to the singlet or the total valence.
*
                  if(mfi.le.5)then
*     V35 = V
                     do k=0,13
                        Ev2EvQCD(jgrid,k,3,alpha,beta) =
     1                       Ev2EvQCD(jgrid,k,3,alpha,beta)
     2                       + Ev2EvQCD(jgrid,k,8,alpha,beta)
                        Ev2EvQCD(jgrid,k,8,alpha,beta) = 0d0
*     T35 = Sigma
                        Ev2EvQCD(jgrid,k,1,alpha,beta) =
     1                       Ev2EvQCD(jgrid,k,1,alpha,beta)
     2                       + Ev2EvQCD(jgrid,k,13,alpha,beta)
                        Ev2EvQCD(jgrid,k,13,alpha,beta) = 0d0
                     enddo
                  endif
*
                  if(mfi.le.4)then
*     V24 = V
                     do k=0,13
                        Ev2EvQCD(jgrid,k,3,alpha,beta) =
     1                       Ev2EvQCD(jgrid,k,3,alpha,beta)
     2                       + Ev2EvQCD(jgrid,k,7,alpha,beta)
                        Ev2EvQCD(jgrid,k,7,alpha,beta) = 0d0
*     T24 = Sigma
                        Ev2EvQCD(jgrid,k,1,alpha,beta) =
     1                       Ev2EvQCD(jgrid,k,1,alpha,beta)
     2                       + Ev2EvQCD(jgrid,k,12,alpha,beta)
                        Ev2EvQCD(jgrid,k,12,alpha,beta) = 0d0
                     enddo
                  endif
*
                  if(mfi.le.3)then
*     V15 = V
                     do k=0,13
                        Ev2EvQCD(jgrid,k,3,alpha,beta) =
     1                       Ev2EvQCD(jgrid,k,3,alpha,beta)
     2                       + Ev2EvQCD(jgrid,k,6,alpha,beta)
                        Ev2EvQCD(jgrid,k,6,alpha,beta) = 0d0
*     T15 = Sigma
                        Ev2EvQCD(jgrid,k,1,alpha,beta) =
     1                       Ev2EvQCD(jgrid,k,1,alpha,beta)
     2                       + Ev2EvQCD(jgrid,k,11,alpha,beta)
                        Ev2EvQCD(jgrid,k,11,alpha,beta) = 0d0
                     enddo
                  endif
               enddo
            enddo
*
            do i=0,13
               do j=0,13
                  Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 0d0
                  Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 0d0
                  do k=0,13
                     Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 
     1                    Ev2PhQCD(jgrid,i-7,j,alpha,beta)
     2                    + Tev2phQCD(mff,i,k)
     3                    * Ev2EvQCD(jgrid,k,j,alpha,beta)
                     do l=0,13
                        Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 
     1                       Ph2PhQCD(jgrid,i-7,j-7,alpha,beta)
     2                       + Tev2phQCD(mff,i,k)
     3                       * Ev2EvQCD(jgrid,k,l,alpha,beta)
     4                       * Tph2evQCD(mfi,l,j)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
*
      return
      end
