************************************************************************
*
*     GetSIATotalCrossSection.f:
*
*     This function returns SIA total cross section as a function of Q
*     in GeV and a given perturbative order "pto" and for a given
*     component.
*     Reference: eq. (2.19) of arXiv:hep-ph/9609377.
*
************************************************************************
      function GetSIATotalCrossSection(pto,Q,comp)
*
      implicit none
*
      include "../commons/TimeLike.h"
      include "../commons/m2th.h"
      include "../commons/consts.h"
      include "../commons/kfacQ.h"
      include "../commons/ColorFactors.h"
      include "../commons/MaxFlavourPDFs.h"
**
*     Input Variables
*
      integer pto
      double precision Q
      character*6 comp
**
*     Internal Variables
*
      integer i,nf,nfi,nff
      integer NC
      double precision Q2,a_QCD,as(2),a_QED,alpha2
      double precision bq(0:6),dq(0:6),bqt(0:6),sumq
      double precision sigma0tot
      double precision lnQ2M2
      double precision Ree
      double precision hbarc2
      parameter(NC=3)
      parameter(hbarc2=0.389379338d6)
**
*     Output Variables
*
      double precision GetSIATotalCrossSection
*
      if(.not.TimeLike)then
         write(6,*) "GetSIATotalCrossSection: this function can be"
         write(6,*) "   called only if the time-like evolution has"
         write(6,*) "   been set."
         write(6,*) "   Set it with 'SetTimeLikeEvolution(true)'."
         write(6,*) "  "
         call exit(-10)
      endif
*
      Q2 = Q * Q
*
      as(1) = a_QCD(Q2)
      as(2) = as(1) * as(1)
*
      alpha2 = ( 4d0 * pi * a_QED(Q2) )**2
*
      call ComputeChargesDIS(Q2,bq,dq,bqt)
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
      if(nf.gt.nfMaxPDFs) nf = nfMaxPDFs
*
      nfi = 1
      nff = nf
      if(comp(1:5).eq."total")then
         nfi = 1
         nff = nf
      elseif(comp(1:5).eq."light")then
         nfi = 1
         nff = 3
      elseif(comp(1:5).eq."charm")then
         nfi = 4
         nff = 4
         if(nf.lt.4)then
            GetSIATotalCrossSection = 0d0
            return
         endif
      elseif(comp(1:6).eq."bottom")then
         nfi = 5
         nff = 5
         if(nf.lt.5)then
            GetSIATotalCrossSection = 0d0
            return
         endif
      elseif(comp(1:3).eq."top")then
         nfi = 6
         nff = 6
         if(nf.lt.6)then
            GetSIATotalCrossSection = 0d0
            return
         endif
      endif
*
      sumq = 0d0
      do i=nfi,nff
         sumq = sumq + bq(0) * bq(i)
      enddo
*
      sigma0tot = 4d0 * pi * alpha2 * NC * sumq / 3d0 / Q2
*
      lnQ2M2 = - dlog(kfacQ)
*
      Ree = 1d0
      if(pto.ge.1)then
         Ree = Ree + as(1) * CF * 3d0
      endif
      if(pto.ge.2)then
         Ree = Ree + as(2) * ( CF**2 * ( - 3d0 / 2d0 )
     1       + CA * CF * ( - 11d0 * lnQ2M2 - 44d0 * zeta3
     2       + 123d0 / 2d0 ) + nf * CF * TR * ( 4d0 * lnQ2M2
     3       + 16d0 * zeta3 - 22d0 ) )
      endif
*
      GetSIATotalCrossSection = Ree * sigma0tot
*
*     Covert cross section in nbarn
*
      GetSIATotalCrossSection = hbarc2 * GetSIATotalCrossSection
*
      return
      end
