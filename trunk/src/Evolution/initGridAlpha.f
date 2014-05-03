************************************************************************
*
*     initGridAlpha.f:
*
*     This routine constructs a grid for the strong coupling constant 
*     alphas to be used for the precomputation of the resummed splitting
*     functions.
*
************************************************************************
      subroutine initGridAlpha
*
      implicit none
*
      include "../commons/gridAlpha.h"
      include "../commons/scales.h"
      include "../commons/m2th.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/Smallx.h"
      include "../commons/ipt.h"
*
*     Variables
*
      integer tau,inf,ig
      integer nfi,nff
      integer nsubgrids,nas(4),ncount
      double precision a_QCD,amax,amin,step
      double precision lbound(4),ubound(4)
      double precision asmh(4:6)
      double precision eps
      parameter(eps=1d-10)
*
*     Initialize grid
*
      do tau=0,na
         ag(tau) = 0d0
      enddo
*
*     Values of alphas at the energy bounds
*
      amin = a_QCD(Q2min)
      amax = a_QCD(Q2max)
*
*     Fixed Flavour Number Scheme
*
      if(Evs.eq."FF")then
*     Linear distribution
         ag(0) = amin
         step = ( amax - amin ) / dble(na)
         do tau=1,na
            ag(tau) = ag(tau-1) + step
         enddo
*     Value of the strong coupling at the heavy quark threshold
         if(Nf_FF.eq.3)then
            asmh(4) = 0d0
            asmh(5) = 0d0
            asmh(6) = 0d0
         elseif(Nf_FF.eq.4)then
            asmh(4) = 1d20
            asmh(5) = 0d0
            asmh(6) = 0d0
         elseif(Nf_FF.eq.5)then
            asmh(4) = 1d20
            asmh(5) = 1d20
            asmh(6) = 0d0
         elseif(Nf_FF.eq.6)then
            asmh(4) = 1d20
            asmh(5) = 1d20
            asmh(6) = 1d20
         endif
*
*     Variable Flavour Number Scheme
*
      elseif(Evs.eq."VF")then
*     Find number of flavours at the energy bounds     
         if(Q2max.gt.m2th(6))then
            nff = 6
         elseif(Q2max.gt.m2th(5))then
            nff = 5
         elseif(Q2max.gt.m2th(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfMaxAlpha) nff = nfMaxAlpha
*
         if(Q2min.gt.m2th(6))then
            nfi = 6
         elseif(Q2min.gt.m2th(5))then
            nfi = 5
         elseif(Q2min.gt.m2th(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfMaxAlpha) nfi = nfMaxAlpha
*     If the final and the initial number of flavours are equal...
         if(nfi.eq.nff)then
            ag(0) = amin
            step = ( amax - amin ) / dble(na)
            do tau=1,na
               ag(tau) = ag(tau-1) + step
            enddo
*     Value of the strong coupling at the heavy quark threshold
            if(nfi.eq.3)then
               asmh(4) = 0d0
               asmh(5) = 0d0
               asmh(6) = 0d0
            elseif(nfi.eq.4)then
               asmh(4) = 1d20
               asmh(5) = 0d0
               asmh(6) = 0d0
            elseif(nfi.eq.5)then
               asmh(4) = 1d20
               asmh(5) = 1d20
               asmh(6) = 0d0
            elseif(nfi.eq.6)then
               asmh(4) = 1d20
               asmh(5) = 1d20
               asmh(6) = 1d20
            endif
         else
*     ...else find the values of alphas at and just below the the mass thresholds
            if(nfi.lt.nff)then
               nsubgrids = 1
               lbound(nsubgrids) = amin
               do inf=nfi+1,nff
                  ubound(nsubgrids) = a_QCD(m2th(inf)-eps)
                  nsubgrids = nsubgrids + 1
                  lbound(nsubgrids) = a_QCD(m2th(inf))
               enddo
               ubound(nsubgrids) = amax
*     Value of the strong coupling at the heavy quark threshold
               if(nfi.gt.3)then
                  do inf=4,nfi
                     asmh(inf) = 1d20
                  enddo
               endif
               if(nff.lt.6)then
                  do inf=nff+1,6
                     asmh(inf) = 0d0
                  enddo
               endif
               do inf=nfi+1,nff
                  asmh(inf) = a_QCD(m2th(inf))
               enddo
            elseif(nfi.gt.nff)then
               nsubgrids = 1
               lbound(nsubgrids) = amin
               do inf=nff+1,nfi
                  ubound(nsubgrids) = a_QCD(m2th(inf)-eps)
                  nsubgrids = nsubgrids + 1
                  lbound(nsubgrids) = a_QCD(m2th(inf))
               enddo
               ubound(nsubgrids) = amax
*     Value of the strong coupling at the heavy quark threshold
               if(nff.gt.3)then
                  do inf=4,nff
                     asmh(inf) = 1d20
                  enddo
               endif
               if(nfi.lt.6)then
                  do inf=nfi+1,6
                     asmh(inf) = 0d0
                  enddo
               endif
               do inf=nff+1,nfi
                  asmh(inf) = a_QCD(m2th(inf))
               enddo
            endif
*     Evaluate the number of nodes of the subgrids according
*     to the coverage of alphas.
            ncount = 0
            do ig=1,nsubgrids
               nas(ig) = max(floor( ( na + 1 )
     1                 * ( ubound(ig) - lbound(ig) )
     2                 / ( amax - amin )) - 1,1)
               ncount = ncount + nas(ig) + 1
            enddo
*     Check that the sum of the number of nodes of the subgrids is
*     equal to "na+1", if not adjust.
            if(ncount.ne.(na+1))then
               if(ncount-na-1.eq.1)then
                  nas(1) = nas(1) - 1
               elseif(ncount-na-1.eq.-1)then
                  nas(nsubgrids) = nas(nsubgrids) + 1
               else
                  write(6,*) "In initGridAlpha.f:"
                  write(6,*) "The number of nodes of the subgrids does",
     1                       " not match:",ncount,na+1
                  call exit(-10)
               endif
            endif
*     Compute the subgrids 
            ncount = 0
            do ig=1,nsubgrids
               ag(ncount) = lbound(ig)
c               write(45,*) ncount,ag(ncount)
               step = ( ubound(ig) - lbound(ig) ) / dble(nas(ig))
               do tau=1,nas(ig)
                  ncount = ncount + 1
                  ag(ncount)  = ag(ncount-1) + step
c                  write(45,*) ncount,ag(ncount)
               enddo
               ncount = ncount + 1
            enddo
         endif
      endif
*
*     Initialize HELL
*
      call initHELL(LogAcc,ipt,asmh(4),asmh(5),asmh(6))
*
      return
      end
