************************************************************************
*
*     TruncatedEvolveAPFEL.f:
*
*     This ruotine computes the evolved PDFs on the grids.
*
************************************************************************
      subroutine TruncatedEvolveAPFEL(Q20,Q2)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/grid.h"
      include "../commons/Th.h"
      include "../commons/FastEvol.h"
      include "../commons/Nf_FF.h"
      include "../commons/Evs.h"
      include "../commons/m2th.h"
      include "../commons/fph.h"
      include "../commons/f0ph.h"
      include "../commons/EpsTrunc.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/pdfset.h"
      include "../commons/EvolOp.h"
      include "../commons/EvolutionMatrices.h"
**
*     Input Variables
*
      double precision Q20,Q2
**
*     Internal Variables
*
      integer inf,mfi,mff,mfmax
      integer ipdf,jpdf,alpha,beta
      integer sign
      integer ieps
      double precision mu2i(3:7),mu2f(3:7)
      double precision fpheps(-2:2,-6:6,0:nint_max)
      double precision fleptoneps(-2:2,-3:3,0:nint_max)
      double precision tiny
      double precision eps(-2:2)
      parameter(tiny=1d-10)
      double precision fqpre(ngrid_max,-6:6,0:nint_max)
      double precision flpre(ngrid_max,-3:3,0:nint_max)
      common / pretabAPFEL / fqpre,flpre
      double precision Msg(-2:2,2,2,0:nint_max,0:nint_max)
      double precision Mnsp(-2:2,0:nint_max,0:nint_max)
      double precision Mnsm(-2:2,0:nint_max,0:nint_max)
      double precision Mnsv(-2:2,0:nint_max,0:nint_max)
      character*100 pdfsetbkp
      logical EvolOpbkp
*
*     For the LO solution use the standard 'EvolveAPFEL'
*
      if(ipt.eq.0)then
         call ExponentiatedEvolveAPFEL(Q20,Q2)
         return
      endif
*
*     Find initial and final number of flavours
*
      if(Evs.eq."FF")then
         mfi = Nf_FF
         mff = Nf_FF
*
         sign = 1
         mu2i(mfi) = Q20
         mu2f(mff) = Q2
      elseif(Evs.eq."VF")then
*
*     Define maximun number of flavours
*
         mfmax = max(nfMaxPDFs,nfMaxAlpha)
*
         if(Q2.gt.m2th(6))then
            mff = 6
         elseif(Q2.gt.m2th(5))then
            mff = 5
         elseif(Q2.gt.m2th(4))then
            mff = 4
         else
            mff = 3
         endif
         if(mff.gt.mfmax) mff = mfmax
*
         if(Q20.gt.m2th(6))then
            mfi = 6
         elseif(Q20.gt.m2th(5))then
            mfi = 5
         elseif(Q20.gt.m2th(4))then
            mfi = 4
         else
            mfi = 3
         endif
         if(mfi.gt.mfmax) mfi = mfmax
*
         sign = 1
         if(Q2.lt.Q20) sign = - 1
*
         EvolOpbkp = EvolOp
         if(EvolOp) EvolOp = .false.
*
*     Define threshold scales
*
         mu2i(mfi) = Q20
         if(sign.eq.1)then
            do inf=mfi+1,mff
               mu2i(inf) = m2th(inf) + tiny
            enddo
            do inf=mfi,mff-1
               mu2f(inf) = m2th(inf+1) - tiny
            enddo
         elseif(sign.eq.-1)then
            do inf=mfi-1,mff,sign
               mu2i(inf) = m2th(inf+1) + tiny
            enddo
            do inf=mfi,mff+1,sign
               mu2f(inf) = m2th(inf) - tiny
            enddo
         endif
         mu2f(mff) = Q2
      endif
*
*     put epsilon in an array
*
      eps(-2) = - 2d0 * EpsTrunc
      eps(-1) = - EpsTrunc
      eps(0)  = 0d0
      eps(1)  = EpsTrunc
      eps(2)  = 2d0 * EpsTrunc
*
*     Compute Evolution
*
      do igrid=1,ngrid
*     Tabulate PDFs at the initial scale on the grid
         if(pdfset.eq."apfel")then
            do alpha=0,nin(igrid)
               do ipdf=-6,6
                  fqpre(igrid,ipdf,alpha) = fph(igrid,ipdf,alpha)
               enddo
               do ipdf=-3,3
                  flpre(igrid,ipdf,alpha) = flepton(igrid,ipdf,alpha)
               enddo
            enddo
         else
            call initPDFs(Q20)
            do alpha=0,nin(igrid)
               do ipdf=-6,6
                  fqpre(igrid,ipdf,alpha) = f0ph(ipdf,alpha)
               enddo
               do ipdf=-3,3
                  flpre(igrid,ipdf,alpha) = f0lep(ipdf,alpha)
               enddo
            enddo
         endif
*
*     Back up PDF name
*
         pdfsetbkp = pdfset
*
*     Set pretabulated PDFs as input
*
         call SetPDFSet("pretabulated")
*     Loop over the number of flavours
         do inf=mfi,mff,sign
*     Loop over the values of epsilon
            do ieps=-2,2
               EpsEff = eps(ieps)
               if(FastEvol)then
*     Evolve directly PDFs on the grid
                  if(Th.eq."QCD")then
                     call EvolutionQCD(mu2i(inf),mu2f(inf))
                  elseif(Th.eq."QED")then
                     call EvolutionQED(mu2i(inf),mu2f(inf))
                  elseif(Th.eq."QUniD")then
                     call EvolutionUnified(mu2i(inf),mu2f(inf))
                  else
                     write(6,*) "The fast evolution is currently",
     1                          " available"
                     write(6,*) "only for the 'QCD','QED','QUniD'",
     1                          " evolutions."
                     write(6,*) "  "
                     call exit(-10)
                  endif
               else
*     Evaluate evolution operators on the grid
                  if(Th.eq."QCD")then
                     call EvolutionOperatorsQCD(mu2i(inf),mu2f(inf))
                  elseif(Th.eq."QED")then
                     call EvolutionOperatorsQED(mu2i(inf),mu2f(inf))
                  elseif(Th.eq."QCEDP".or.Th.eq."QCEDS".or.
     1                    Th.eq."QECDP".or.Th.eq."QECDS".or.
     2                    Th.eq."QavDP".or.Th.eq."QavDS")then
                     call EvolutionOperatorsQCD(mu2i(inf),mu2f(inf))
                     call EvolutionOperatorsQED(mu2i(inf),mu2f(inf))
                  elseif(Th.eq."QUniD")then
                     call EvolutionOperatorsUnified(mu2i(inf),mu2f(inf))
                  endif
*
*     If the evolution operator computation has been enabled,
*     save evolution operators for each "eps(ieps)".
*
                  if(EvolOpbkp)then
                     do alpha=0,nin(igrid)
                        do beta=0,nin(igrid)
                           do ipdf=1,2
                              do jpdf=1,2
                                 Msg(ieps,ipdf,jpdf,alpha,beta) = 
     1                                MQCDsg(inf,ipdf,jpdf,alpha,beta)
                              enddo
                           enddo
                           Mnsp(ieps,alpha,beta) = 
     1                          MQCDnsp(inf,alpha,beta)
                           Mnsm(ieps,alpha,beta) =
     1                          MQCDnsm(inf,alpha,beta)
                           Mnsv(ieps,alpha,beta) =
     1                          MQCDnsv(inf,alpha,beta)
                        enddo
                     enddo
                  endif
*     Initialize PDFs at the initial scale on the grid
                  call initPDFs(mu2i(inf))
*     Convolute intial PDFs with the evolution operators
                  call EvolvePDFs(igrid)
               endif
*
*     Save PDFs for each "eps(ieps)"
*
               do alpha=0,nin(igrid)
                  do ipdf=-6,6
                     fpheps(ieps,ipdf,alpha) = 
     1                    fph(igrid,ipdf,alpha)
                  enddo
                  do ipdf=-3,3
                     fleptoneps(ieps,ipdf,alpha) = 
     1                    flepton(igrid,ipdf,alpha)
                  enddo
               enddo
            enddo
*
*     Combine PDFs
*
*     NLO
*
            if(ipt.ge.1)then
               do alpha=0,nin(igrid)
                  do ipdf=-6,6
                     fqpre(igrid,ipdf,alpha) = fpheps(0,ipdf,alpha)
     1                    + ( - fpheps(2,ipdf,alpha)
     2                    + 8d0 * fpheps(1,ipdf,alpha)
     3                    - 8d0 * fpheps(-1,ipdf,alpha)
     4                    + fpheps(-2,ipdf,alpha) )
     5                    / 12d0 / EpsTrunc
                  enddo
                  do ipdf=-3,3
                     flpre(igrid,ipdf,alpha) = 
     1                    fleptoneps(0,ipdf,alpha)
     2                    + ( - fleptoneps(2,ipdf,alpha)
     3                    + 8d0 * fleptoneps(1,ipdf,alpha)
     4                    - 8d0 * fleptoneps(-1,ipdf,alpha)
     5                    + fleptoneps(-2,ipdf,alpha))
                  enddo
               enddo
            endif
*
*     NNLO
*
            if(ipt.ge.2)then
               do alpha=0,nin(igrid)
                  do ipdf=-6,6
                     fqpre(igrid,ipdf,alpha) = fqpre(igrid,ipdf,alpha)
     1                    + ( - fpheps(2,ipdf,alpha)
     2                    + 16d0 * fpheps(1,ipdf,alpha)
     3                    - 30d0 * fpheps(0,ipdf,alpha)
     4                    + 16d0 * fpheps(-1,ipdf,alpha)
     5                    - fpheps(-2,ipdf,alpha))
     6                    / 24d0 / EpsTrunc / EpsTrunc
                  enddo
                  do ipdf=-3,3
                     flpre(igrid,ipdf,alpha) = 
     1                    flpre(igrid,ipdf,alpha)
     2                    + ( - fleptoneps(2,ipdf,alpha)
     3                    + 16d0 * fleptoneps(1,ipdf,alpha)
     4                    - 30d0 * fleptoneps(0,ipdf,alpha)
     5                    + 16d0 * fleptoneps(-1,ipdf,alpha) 
     6                    - fleptoneps(-2,ipdf,alpha) )
     7                    / 24d0 / EpsTrunc / EpsTrunc
                  enddo
               enddo
            endif
*
*     Match PDFs at the thresholds
*
            if(sign.eq.1.and.inf.lt.mff)then
               EpsEff = 1d0
               if(FastEvol)then
                  if(Th.eq."QCD")then
                     call EvolutionQCD(mu2f(inf),mu2i(inf+1))
                  elseif(Th.eq."QED")then
                     call EvolutionQED(mu2f(inf),mu2i(inf+1))
                  elseif(Th.eq."QUniD")then
                     call EvolutionUnified(mu2f(inf),mu2i(inf+1))
                  endif
               else
                  if(Th.eq."QCD")then
                     call EvolutionOperatorsQCD(mu2f(inf),mu2i(inf+1))
                  elseif(Th.eq."QED")then
                     call EvolutionOperatorsQED(mu2f(inf),mu2i(inf+1))
                  elseif(Th.eq."QCEDP".or.Th.eq."QCEDS".or.
     1                    Th.eq."QECDP".or.Th.eq."QECDS".or.
     2                    Th.eq."QavDP".or.Th.eq."QavDS")then
                     call EvolutionOperatorsQCD(mu2f(inf),mu2i(inf+1))
                     call EvolutionOperatorsQED(mu2f(inf),mu2i(inf+1))
                  elseif(Th.eq."QUniD")then
                     call EvolutionOperatorsUnified(mu2f(inf),
     1                                              mu2i(inf+1))
                  endif
                  call initPDFs(mu2f(inf))
                  call EvolvePDFs(igrid)
               endif
*
*     Copy PDFs into the pretabulated arrays
*
               do alpha=0,nin(igrid)
                  do ipdf=-6,6
                     fqpre(igrid,ipdf,alpha) = fph(igrid,ipdf,alpha)
                  enddo
                  do ipdf=-3,3
                     flpre(igrid,ipdf,alpha) = flepton(igrid,ipdf,alpha)
                  enddo
               enddo
            endif
*
*     Combine evolution operators
*
*     NLO
*
            if(EvolOpbkp)then
               if(ipt.ge.1)then
                  do alpha=0,nin(igrid)
                     do beta=0,nin(igrid)
                        do ipdf=1,2
                           do jpdf=1,2
                              MQCDsg(inf,ipdf,jpdf,alpha,beta) =
     1                             Msg(0,ipdf,jpdf,alpha,beta)
     2                             + ( - Msg(2,ipdf,jpdf,alpha,beta)
     3                             + 8d0 * Msg(1,ipdf,jpdf,alpha,beta)
     4                             - 8d0 * Msg(-1,ipdf,jpdf,alpha,beta)
     5                             + Msg(-2,ipdf,jpdf,alpha,beta) )
     6                             / 12d0 / EpsTrunc 
                           enddo
                        enddo
                        MQCDnsp(inf,alpha,beta) =
     1                       Mnsp(0,alpha,beta)
     2                       + ( - Mnsp(2,alpha,beta)
     3                       + 8d0 * Mnsp(1,alpha,beta)
     4                       - 8d0 * Mnsp(-1,alpha,beta)
     5                       + Mnsp(-2,alpha,beta) )
     6                       / 12d0 / EpsTrunc 
                        MQCDnsm(inf,alpha,beta) =
     1                       Mnsm(0,alpha,beta)
     2                       + ( - Mnsm(2,alpha,beta)
     3                       + 8d0 * Mnsm(1,alpha,beta)
     4                       - 8d0 * Mnsm(-1,alpha,beta)
     5                       + Mnsm(-2,alpha,beta) )
     6                       / 12d0 / EpsTrunc 
                        MQCDnsv(inf,alpha,beta) =
     1                       Mnsv(0,alpha,beta)
     2                       + ( - Mnsv(2,alpha,beta)
     3                       + 8d0 * Mnsv(1,alpha,beta)
     4                       - 8d0 * Mnsv(-1,alpha,beta)
     5                       + Mnsv(-2,alpha,beta) )
     6                       / 12d0 / EpsTrunc 
                     enddo
                  enddo
               endif
*
*     NNLO
*
               if(ipt.ge.2)then
                  do alpha=0,nin(igrid)
                     do beta=0,nin(igrid)
                        do ipdf=1,2
                           do jpdf=1,2
                              MQCDsg(inf,ipdf,jpdf,alpha,beta) =
     1                             MQCDsg(inf,ipdf,jpdf,alpha,beta)
     2                             + ( - Msg(2,ipdf,jpdf,alpha,beta)
     3                             + 16d0 * Msg(1,ipdf,jpdf,alpha,beta)
     4                             - 30d0 * Msg(0,ipdf,jpdf,alpha,beta)
     5                             + 16d0 * Msg(-1,ipdf,jpdf,alpha,beta)
     6                             - Msg(-2,ipdf,jpdf,alpha,beta) )
     7                             / 24d0 / EpsTrunc / EpsTrunc 
                           enddo
                        enddo
                        MQCDnsp(inf,alpha,beta) =
     1                       MQCDnsp(inf,alpha,beta)
     2                       + ( - Mnsp(2,alpha,beta)
     3                       + 16d0 * Mnsp(1,alpha,beta)
     4                       - 30d0 * Mnsp(0,alpha,beta)
     5                       + 16d0 * Mnsp(-1,alpha,beta) 
     6                       - Mnsp(-2,alpha,beta))
     7                       / 24d0 / EpsTrunc / EpsTrunc 
                        MQCDnsm(inf,alpha,beta) =
     1                       MQCDnsm(inf,alpha,beta)
     2                       + ( - Mnsm(2,alpha,beta)
     3                       + 16d0 * Mnsm(1,alpha,beta)
     4                       - 30d0 * Mnsm(0,alpha,beta)
     5                       + 16d0 * Mnsm(-1,alpha,beta) 
     6                       - Mnsm(-2,alpha,beta))
     7                       / 24d0 / EpsTrunc / EpsTrunc 
                        MQCDnsv(inf,alpha,beta) =
     1                       MQCDnsv(inf,alpha,beta)
     2                       + ( - Mnsv(2,alpha,beta)
     3                       + 16d0 * Mnsv(1,alpha,beta)
     4                       - 30d0 * Mnsv(0,alpha,beta)
     5                       + 16d0 * Mnsv(-1,alpha,beta) 
     6                       - Mnsv(-2,alpha,beta))
     7                       / 24d0 / EpsTrunc / EpsTrunc 
                     enddo
                  enddo
               endif
            endif
         enddo
*
*     Copy combined PDFs into the standard evolved PDFs
*
         do alpha=0,nin(igrid)
            do ipdf=-6,6
                fph(igrid,ipdf,alpha) = fqpre(igrid,ipdf,alpha)
            enddo
            fgamma(igrid,alpha) = flpre(igrid,0,alpha)
            do ipdf=-3,3
                flepton(igrid,ipdf,alpha) = flpre(igrid,ipdf,alpha)
            enddo
         enddo
*
*     Restore PDF name
*
         call SetPDFSet(pdfsetbkp)
*
*     Join evolution operators
*
         if(EvolOpbkp)then
            nfi = mfi
            nff = mff
            call JoinOperatorsQCD(igrid)
         endif
      enddo
*     Join all the subgrids.
*     Join also the operators if the production of the evolution operators
*     has been enabled (For the moment available only for QCD).
      EvolOp = EvolOpbkp
      call JoinGrids
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
      double precision fqpre(ngrid_max,-6:6,0:nint_max)
      double precision flpre(ngrid_max,-3:3,0:nint_max)
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
