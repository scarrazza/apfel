************************************************************************
*
*     EvolutionUnified.f:
*
*     This routine returns the singlet and the non-singlet evolution
*     operators on the x-space grid between the scales m20 and mu2 for
*     nf active flavours given the initial evolution operators in the
*     unified evolution basis.
*
************************************************************************
      subroutine EvolutionUnified(muF20,muF2)
*
      implicit none
*
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/m2th.h"
      include "../commons/grid.h"
      include "../commons/wrap.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/f0ph.h"
      include "../commons/fph.h"
**
*     Input Variables
*
      double precision muF20,muF2
**
*     Internal Variables
*
      integer inf
      integer i,alpha
      integer sgn
      integer nfi,nff
      double precision mu2i(3:6),mu2f(4:7)
      double precision f0ev(0:13,0:nint_max)
      double precision fev(0:13,0:nint_max)
      double precision fpy(-7:6,0:nint_max)
      double precision fevQCD(0:13,0:nint_max)
      double precision Sg1(4,0:nint_max),Sg10(4,0:nint_max)
      double precision Sg2(2,0:nint_max),Sg20(2,0:nint_max)
      double precision Tu1(0:nint_max),Tu10(0:nint_max)
      double precision Tu2(0:nint_max),Tu20(0:nint_max)
      double precision Td1(0:nint_max),Td10(0:nint_max)
      double precision Td2(0:nint_max),Td20(0:nint_max)
      double precision Vu1(0:nint_max),Vu10(0:nint_max)
      double precision Vu2(0:nint_max),Vu20(0:nint_max)
      double precision Vd1(0:nint_max),Vd10(0:nint_max)
      double precision Vd2(0:nint_max),Vd20(0:nint_max)
      double precision tiny
      parameter(tiny=1d-10)
*
*     Initial Conditions
*
      call initPDFs(muF20)
*
*     Rotate initial PDFs into the evolution basis
*
      call PDFphys2evUni(f0bos,f0ph,f0ev)
*
*     Mass scheme
*
*     Fixed Flavour Number Scheme
*
      if(Evs.eq."FF")then
         wnf = Nf_FF
*     If initial and final energies are equal, return immediately the intial conditions
         if(muF2.eq.muF20)then
            do alpha=0,nin(igrid)
               do i=0,13
                  fev(i,alpha) = f0ev(i,alpha)
               enddo
            enddo
            goto 104
         endif
*
         do alpha=0,nin(igrid)
            Sg10(1,alpha) = f0ev(0,alpha)
            Sg10(2,alpha) = f0ev(1,alpha)
            Sg10(3,alpha) = f0ev(2,alpha)
            Sg10(4,alpha) = f0ev(3,alpha)
            Tu10(alpha)   = f0ev(4,alpha)
            Tu20(alpha)   = f0ev(5,alpha)
            Td10(alpha)   = f0ev(6,alpha)
            Td20(alpha)   = f0ev(7,alpha)
            Sg20(1,alpha) = f0ev(8,alpha)
            Sg20(2,alpha) = f0ev(9,alpha)
            Vu10(alpha)   = f0ev(10,alpha)
            Vu20(alpha)   = f0ev(11,alpha)
            Vd10(alpha)   = f0ev(12,alpha)
            Vd20(alpha)   = f0ev(13,alpha)
         enddo
*     Singlet 1
         call odeintsgUnifiedfS1(muF20,muF2,Sg10,Sg1)
*     Singlet 2
         call odeintsgUnifiedfS2(muF20,muF2,Sg20,Sg2)
*     Non-Singlet
         call odeintnsUnifiedf(2,muF20,muF2,Td1,Td10)
         call odeintnsUnifiedf(4,muF20,muF2,Vd1,Vd10)
         if(Nf_FF.eq.3)then
            do alpha=0,nin(igrid)
               Tu1(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
               Tu2(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
               Td2(alpha) = ( Sg1(3,alpha) - Sg1(4,alpha) ) / 2d0
               Vu1(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
               Vu2(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
               Vd2(alpha) = ( Sg2(1,alpha) - Sg2(2,alpha) ) / 2d0
            enddo
         endif
         if(Nf_FF.eq.4)then
            call odeintnsUnifiedf(1,muF20,muF2,Tu1,Tu10)
            call odeintnsUnifiedf(3,muF20,muF2,Vu1,Vu10)
            do alpha=0,nin(igrid)
               Tu2(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
               Td2(alpha) = ( Sg1(3,alpha) - Sg1(4,alpha) ) / 2d0
               Vu2(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
               Vd2(alpha) = ( Sg2(1,alpha) - Sg2(2,alpha) ) / 2d0
            enddo
         endif
         if(Nf_FF.eq.5)then
            call odeintnsUnifiedf(1,muF20,muF2,Tu1,Tu10)
            call odeintnsUnifiedf(2,muF20,muF2,Td2,Td20)
            call odeintnsUnifiedf(3,muF20,muF2,Vu1,Vu10)
            call odeintnsUnifiedf(4,muF20,muF2,Vd2,Vd20)
            do alpha=0,nin(igrid)
               Tu2(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
               Vu2(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
            enddo
         endif
         if(Nf_FF.eq.6)then
            call odeintnsUnifiedf(1,muF20,muF2,Tu1,Tu10)
            call odeintnsUnifiedf(1,muF20,muF2,Tu2,Tu20)
            call odeintnsUnifiedf(2,muF20,muF2,Td2,Td20)
            call odeintnsUnifiedf(3,muF20,muF2,Vu1,Vu10)
            call odeintnsUnifiedf(3,muF20,muF2,Vu2,Vu20)
            call odeintnsUnifiedf(4,muF20,muF2,Vd2,Vd20)
         endif
*
*     Variable Flavour Number Scheme
*
      elseif(Evs.eq."VF")then
*
*     Find initial and final number of flavours
*
         if(muF2.gt.m2th(6))then
            nff = 6
         elseif(muF2.gt.m2th(5))then
            nff = 5
         elseif(muF2.gt.m2th(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfMaxPDFs) nff = nfMaxPDFs
*
         if(muF20.gt.m2th(6))then
            nfi = 6
         elseif(muF20.gt.m2th(5))then
            nfi = 5
         elseif(muF20.gt.m2th(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfMaxPDFs) nfi = nfMaxPDFs
*     If initial and final energies are equal, return immediately the intial conditions
         if(muF2.eq.muF20)then
            do alpha=0,nin(igrid)
               do i=0,13
                  fev(i,alpha) = f0ev(i,alpha)
               enddo
            enddo
            goto 104
         elseif(muF2.gt.muF20)then
            sgn = 1
         elseif(muF2.lt.muF20)then
            sgn = - 1
         endif
*
         mu2i(nfi) = muF20
         if(sgn.eq.1)then
            do inf=nfi+1,nff
               mu2i(inf) = m2th(inf)
            enddo
            do inf=nfi,nff-1
               mu2f(inf) = m2th(inf+1) - tiny
            enddo
         elseif(sgn.eq.-1)then
            do inf=nfi-1,nff,sgn
               mu2i(inf) = m2th(inf+1) + tiny
            enddo
            do inf=nfi,nff+1,sgn
               mu2f(inf) = m2th(inf)
            enddo
         endif
         mu2f(nff) = muF2
*
         do inf=nfi,nff,sgn
            wnf = inf
*
*     Apply matching conditions for the backward evolution
*
            if(sgn.eq.-1)then
               if(inf.lt.nfi)then
*     Rotate to the QCD evolution basis
                  call PDFevUni2evQCD(f0ev,fevQCD)
*     Apply matching conditions
                  call MatchPDFs(inf,fevQCD)
*     Rotate back to the unified evolution basis
                  call PDFevQCD2evUni(fevQCD,f0ev)
               endif
            endif
*
            do alpha=0,nin(igrid)
               Sg10(1,alpha) = f0ev(0,alpha)
               Sg10(2,alpha) = f0ev(1,alpha)
               Sg10(3,alpha) = f0ev(2,alpha)
               Sg10(4,alpha) = f0ev(3,alpha)
               Tu10(alpha)   = f0ev(4,alpha)
               Tu20(alpha)   = f0ev(5,alpha)
               Td10(alpha)   = f0ev(6,alpha)
               Td20(alpha)   = f0ev(7,alpha)
               Sg20(1,alpha) = f0ev(8,alpha)
               Sg20(2,alpha) = f0ev(9,alpha)
               Vu10(alpha)   = f0ev(10,alpha)
               Vu20(alpha)   = f0ev(11,alpha)
               Vd10(alpha)   = f0ev(12,alpha)
               Vd20(alpha)   = f0ev(13,alpha)
            enddo
*     Singlet 1
            call odeintsgUnifiedfS1(mu2i(inf),mu2f(inf),Sg10,Sg1)
*     Singlet 2
            call odeintsgUnifiedfS2(mu2i(inf),mu2f(inf),Sg20,Sg2)
*     Non-Singlet
            call odeintnsUnifiedf(2,mu2i(inf),mu2f(inf),Td1,Td10)
            call odeintnsUnifiedf(4,mu2i(inf),mu2f(inf),Vd1,Vd10)
            if(inf.eq.3)then
               do alpha=0,nin(igrid)
                  Tu1(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
                  Tu2(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
                  Td2(alpha) = ( Sg1(3,alpha) - Sg1(4,alpha) ) / 2d0
                  Vu1(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
                  Vu2(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
                  Vd2(alpha) = ( Sg2(1,alpha) - Sg2(2,alpha) ) / 2d0
               enddo
            endif
            if(inf.eq.4)then
               call odeintnsUnifiedf(1,mu2i(inf),mu2f(inf),Tu1,Tu10)
               call odeintnsUnifiedf(3,mu2i(inf),mu2f(inf),Vu1,Vu10)
               do alpha=0,nin(igrid)
                  Tu2(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
                  Td2(alpha) = ( Sg1(3,alpha) - Sg1(4,alpha) ) / 2d0
                  Vu2(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
                  Vd2(alpha) = ( Sg2(1,alpha) - Sg2(2,alpha) ) / 2d0
               enddo
            endif
            if(inf.eq.5)then
               call odeintnsUnifiedf(1,mu2i(inf),mu2f(inf),Tu1,Tu10)
               call odeintnsUnifiedf(2,mu2i(inf),mu2f(inf),Td2,Td20)
               call odeintnsUnifiedf(3,mu2i(inf),mu2f(inf),Vu1,Vu10)
               call odeintnsUnifiedf(4,mu2i(inf),mu2f(inf),Vd2,Vd20)
               do alpha=0,nin(igrid)
                  Tu2(alpha) = ( Sg1(3,alpha) + Sg1(4,alpha) ) / 2d0
                  Vu2(alpha) = ( Sg2(1,alpha) + Sg2(2,alpha) ) / 2d0
               enddo
            endif
            if(inf.eq.6)then
               call odeintnsUnifiedf(1,mu2i(inf),mu2f(inf),Tu1,Tu10)
               call odeintnsUnifiedf(1,mu2i(inf),mu2f(inf),Tu2,Tu20)
               call odeintnsUnifiedf(2,mu2i(inf),mu2f(inf),Td2,Td20)
               call odeintnsUnifiedf(3,mu2i(inf),mu2f(inf),Vu1,Vu10)
               call odeintnsUnifiedf(3,mu2i(inf),mu2f(inf),Vu2,Vu20)
               call odeintnsUnifiedf(4,mu2i(inf),mu2f(inf),Vd2,Vd20)
            endif
*
*     Update initial coditions for the next step
*
            do alpha=0,nin(igrid)
               f0ev(0,alpha)  = Sg1(1,alpha)
               f0ev(1,alpha)  = Sg1(2,alpha)
               f0ev(2,alpha)  = Sg1(3,alpha)
               f0ev(3,alpha)  = Sg1(4,alpha)
               f0ev(4,alpha)  = Tu1(alpha)
               f0ev(5,alpha)  = Tu2(alpha)
               f0ev(6,alpha)  = Td1(alpha)
               f0ev(7,alpha)  = Td2(alpha)
               f0ev(8,alpha)  = Sg2(1,alpha)
               f0ev(9,alpha)  = Sg2(2,alpha)
               f0ev(10,alpha) = Vu1(alpha)
               f0ev(11,alpha) = Vu2(alpha)
               f0ev(12,alpha) = Vd1(alpha)
               f0ev(13,alpha) = Vd2(alpha)
            enddo
*
*     Apply matching conditions for the forward evolution
*
            if(sgn.eq.1)then
               if(inf.lt.nff)then
*     Rotate to the QCD evolution basis
                  call PDFevUni2evQCD(f0ev,fevQCD)
*     Apply matching conditions
                  call MatchPDFs(inf+1,fevQCD)
*     Rotate back to the unified evolution basis
                  call PDFevQCD2evUni(fevQCD,f0ev)
               endif
            endif
         enddo
      endif
*
*     Put evolved PDFs in a single vector
*
      do alpha=0,nin(igrid)
         fev(0,alpha)  = Sg1(1,alpha)
         fev(1,alpha)  = Sg1(2,alpha)
         fev(2,alpha)  = Sg1(3,alpha)
         fev(3,alpha)  = Sg1(4,alpha)
         fev(4,alpha)  = Tu1(alpha)
         fev(5,alpha)  = Tu2(alpha)
         fev(6,alpha)  = Td1(alpha)
         fev(7,alpha)  = Td2(alpha)
         fev(8,alpha)  = Sg2(1,alpha)
         fev(9,alpha)  = Sg2(2,alpha)
         fev(10,alpha) = Vu1(alpha)
         fev(11,alpha) = Vu2(alpha)
         fev(12,alpha) = Vd1(alpha)
         fev(13,alpha) = Vd2(alpha)
      enddo
*
*     Rotate PDFs back into the physical basis
*
 104  call PDFevUni2phys(fev,fpy)
*
*     Put evolved PDFs in the physical basis in the common used for the interpolation
*
      do alpha=0,nin(igrid)
         do i=-6,6
            fph(igrid,i,alpha)  = fpy(i,alpha)
         enddo
         fgamma(igrid,alpha) = fpy(-7,alpha)
      enddo
*
      return
      end
