************************************************************************
*
*     TruncatedEvolveAPFEL.f:
*
*     This ruotine computes the evolved PDFs on the grids.
*
************************************************************************
      subroutine TruncatedEvolveAPFEL(Q0,Q)
*
      implicit none
*
      include "../commons/InAPFEL.h"
      include "../commons/PDFEvolution.h"
      include "../commons/ipt.h"
      include "../commons/scales.h"
      include "../commons/grid.h"
      include "../commons/Th.h"
      include "../commons/FastEvol.h"
      include "../commons/Nf_FF.h"
      include "../commons/Evs.h"
      include "../commons/m2th.h"
      include "../commons/fph.h"
      include "../commons/EpsTrunc.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/pdfset.h"
**
*     Input Variables
*
      double precision Q0,Q
**
*     Internal Variables
*
      integer inf,nfi,nff
      integer jgrid,ipdf,alpha
      integer sgn
      double precision mu2i(3:7),mu2f(3:7)
      double precision Q20,Q2
      double precision t1,t2
      double precision fpheps(-1:1,0:ngrid_max,-6:6,0:nint_max)
      double precision fgammaeps(-1:1,0:ngrid_max,0:nint_max)
      double precision fleptoneps(-1:1,0:ngrid_max,-3:3,0:nint_max)
      double precision tiny
      parameter(tiny=1d-10)
      double precision fqpre(0:ngrid_max,-6:6,0:nint_max)
      double precision flpre(0:ngrid_max,-3:3,0:nint_max)
      common / pretabAPFEL / fqpre,flpre
      character*100 pdfsetbkp
*
*     Check that truncated solution has been set
*
      if(PDFevol(1:9).ne."truncated")then
         write(6,*) "TruncatedEvolveAPFEL: this function can be used",
     1              " only if the 'truncated' solution has been set."
         write(6,*) "call SetPDFEvolution('truncated') before calling",
     1              " 'TruncatedEvolveAPFEL'"
         write(6,*) "   "
         call exit(-10)
      endif
*
*     For the LO solution use the standard 'EvolveAPFEL'
*
      if(ipt.eq.0)then
         call EvolveAPFEL(Q0,Q)
         return
      endif
*
*     Define initial and final number of flavours
*
      Q20 = Q0 * Q0
      Q2  = Q * Q
      if(Evs.eq."FF")then
         nfi = Nf_FF
         nff = Nf_FF
*
         sgn = 1
         mu2i(nfi) = Q20
         mu2f(nff) = Q2
      elseif(Evs.eq."VF")then
*
*     Find initial and final number of flavours
*
         if(Q2.gt.m2th(6))then
            nff = 6
         elseif(Q2.gt.m2th(5))then
            nff = 5
         elseif(Q2.gt.m2th(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfMaxPDFs) nff = nfMaxPDFs
*
         if(Q20.gt.m2th(6))then
            nfi = 6
         elseif(Q20.gt.m2th(5))then
            nfi = 5
         elseif(Q20.gt.m2th(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfMaxPDFs) nfi = nfMaxPDFs
*
         sgn = 1
         if(Q2.lt.Q20) sgn = - 1
*
         mu2i(nfi) = Q20
         if(sgn.eq.1)then
            do inf=nfi+1,nff
               mu2i(inf) = m2th(inf) + tiny
            enddo
            do inf=nfi,nff-1
               mu2f(inf) = m2th(inf+1) - tiny
            enddo
         elseif(sgn.eq.-1)then
            do inf=nfi-1,nff,sgn
               mu2i(inf) = m2th(inf+1) + tiny
            enddo
            do inf=nfi,nff+1,sgn
               mu2f(inf) = m2th(inf) - tiny
            enddo
         endif
         mu2f(nff) = Q2
      endif
*
      pdfsetbkp = pdfset
      do inf=nfi,nff,sgn
         EpsEff = 0d0
         call EvolveAPFEL(dsqrt(mu2i(inf)),dsqrt(mu2f(inf)))
         do jgrid=0,ngrid
            do alpha=0,nin(jgrid)
               do ipdf=-6,6
                  fpheps(0,jgrid,ipdf,alpha) = fph(jgrid,ipdf,alpha)
               enddo
               fgammaeps(0,jgrid,alpha) = fgamma(jgrid,alpha)
               do ipdf=-3,3
                  fleptoneps(0,jgrid,ipdf,alpha) = 
     1                 flepton(jgrid,ipdf,alpha)
               enddo
            enddo
         enddo
         EpsEff = EpsTrunc
         call EvolveAPFEL(dsqrt(mu2i(inf)),dsqrt(mu2f(inf)))
         do jgrid=0,ngrid
            do alpha=0,nin(jgrid)
               do ipdf=-6,6
                  fpheps(1,jgrid,ipdf,alpha) = fph(jgrid,ipdf,alpha)
               enddo
               fgammaeps(1,jgrid,alpha) = fgamma(jgrid,alpha)
               do ipdf=-3,3
                  fleptoneps(1,jgrid,ipdf,alpha) = 
     1                 flepton(jgrid,ipdf,alpha)
               enddo
            enddo
         enddo
         EpsEff = - EpsTrunc
         call EvolveAPFEL(dsqrt(mu2i(inf)),dsqrt(mu2f(inf)))
         do jgrid=0,ngrid
            do alpha=0,nin(jgrid)
               do ipdf=-6,6
                  fpheps(-1,jgrid,ipdf,alpha) = fph(jgrid,ipdf,alpha)
               enddo
               fgammaeps(-1,jgrid,alpha) = fgamma(jgrid,alpha)
               do ipdf=-3,3
                  fleptoneps(-1,jgrid,ipdf,alpha) = 
     1                 flepton(jgrid,ipdf,alpha)
               enddo
            enddo
         enddo
*
*     Combine PDFs
*
*     NLO
*
         if(ipt.ge.1)then
            do jgrid=0,ngrid
               do alpha=0,nin(jgrid)
                  do ipdf=-6,6
                     fph(jgrid,ipdf,alpha) = fpheps(0,jgrid,ipdf,alpha)
     1                    + ( fpheps(1,jgrid,ipdf,alpha)
     2                    - fpheps(-1,jgrid,ipdf,alpha) )
     3                    / 2d0 / EpsTrunc
                  enddo
                  fgamma(jgrid,alpha) = fgammaeps(0,jgrid,alpha)
     1                 + ( fgammaeps(1,jgrid,alpha)
     2                 - fgammaeps(-1,jgrid,alpha) )
     3                 / 2d0 / EpsTrunc
                  do ipdf=-3,3
                     flepton(jgrid,ipdf,alpha) = 
     1                    fleptoneps(0,jgrid,ipdf,alpha)
     2                    + ( fleptoneps(1,jgrid,ipdf,alpha)
     3                    - fleptoneps(-1,jgrid,ipdf,alpha) )
     4                    / 2d0 / EpsTrunc
                  enddo
               enddo
            enddo
         endif
*
*     NNLO
*
         if(ipt.ge.2)then
            do jgrid=0,ngrid
               do alpha=0,nin(jgrid)
                  do ipdf=-6,6
                     fph(jgrid,ipdf,alpha) = fph(jgrid,ipdf,alpha)
     1                    + ( fpheps(1,jgrid,ipdf,alpha)
     2                    - 2d0 * fpheps(0,jgrid,ipdf,alpha)
     3                    + fpheps(-1,jgrid,ipdf,alpha) )
     4                    / 2d0 / EpsTrunc / EpsTrunc
                  enddo
                  fgamma(jgrid,alpha) = fgamma(jgrid,alpha)
     1                 + ( fgammaeps(1,jgrid,alpha)
     2                 - 2d0 * fgammaeps(0,jgrid,alpha)
     2                 + fgammaeps(-1,jgrid,alpha) )
     3                 / 2d0 / EpsTrunc / EpsTrunc
                  do ipdf=-3,3
                     flepton(jgrid,ipdf,alpha) = 
     1                    flepton(jgrid,ipdf,alpha)
     2                    + ( fleptoneps(1,jgrid,ipdf,alpha)
     3                    - 2d0 * fleptoneps(0,jgrid,ipdf,alpha)
     3                    + fleptoneps(-1,jgrid,ipdf,alpha) )
     4                    / 2d0 / EpsTrunc / EpsTrunc
                  enddo
               enddo
            enddo
         endif
         call SetPDFSet("apfel")
*     Matching at the thresholds
         if(inf.lt.nff) call EvolveAPFEL(dsqrt(mu2f(inf)),
     1                                   dsqrt(mu2i(inf+1)))
*
*     Save pretabulated PDFs
*
         do jgrid=0,ngrid
            do alpha=0,nin(jgrid)
               do ipdf=-6,6
                  fqpre(jgrid,ipdf,alpha) = fph(jgrid,ipdf,alpha)
               enddo
               do ipdf=-3,3
                  flpre(jgrid,ipdf,alpha) = flepton(jgrid,ipdf,alpha)
               enddo
            enddo
         enddo
         call SetPDFSet("pretabulated")
      enddo
      pdfset = pdfsetbkp
*
      return
      end
*
************************************************************************
      subroutine pretabulatedPDFs(jgrid,alpha,fq,fl)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer jgrid,alpha
**
*     Internal Variables
*
      integer ifl,ilept
      double precision fqpre(0:ngrid_max,-6:6,0:nint_max)
      double precision flpre(0:ngrid_max,-3:3,0:nint_max)
      common / pretabAPFEL / fqpre,flpre
**
*     Output Variables
*
      double precision fq(-6:6)
      double precision fl(-3:3)
*
      do ifl=-6,6
         fq(ifl) = fqpre(jgrid,ifl,alpha)
      enddo
      do ilept=-3,3
         fl(ilept) = flpre(jgrid,ilept,alpha)
      enddo
*
      return
      end
