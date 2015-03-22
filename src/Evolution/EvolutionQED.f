************************************************************************
*
*     EvolutionQED.f:
*
*     This routine returns the singlet and the non-singlet evolution
*     operators on the x-space grid between the scales m20 and mu2 for
*     nf active flavours given the initial evolution operators in QED.
*
************************************************************************
      subroutine EvolutionQED(muF20,muF2)
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
      double precision mu2i(3:7),mu2f(3:7)
      double precision f0ev(0:13,0:nint_max)
      double precision fev(0:13,0:nint_max)
      double precision fpy(-6:6,0:nint_max)
      double precision Sg(3,0:nint_max),Sg0(3,0:nint_max)
      double precision Duc(0:nint_max),Duc0(0:nint_max)
      double precision Dds(0:nint_max),Dds0(0:nint_max)
      double precision Dsb(0:nint_max),Dsb0(0:nint_max)
      double precision Dct(0:nint_max),Dct0(0:nint_max)
      double precision um(0:nint_max),um0(0:nint_max)
      double precision dm(0:nint_max),dm0(0:nint_max)
      double precision sm(0:nint_max),sm0(0:nint_max)
      double precision cm(0:nint_max),cm0(0:nint_max)
      double precision bm(0:nint_max),bm0(0:nint_max)
      double precision tm(0:nint_max),tm0(0:nint_max)
*
*     Initial Conditions
*
      call initPDFs(muF20)
*
*     Rotate initial PDFs into the evolution basis
*
      call switchGluonPhoton
      call PDFphys2evQED(f0ph,f0ev)
*
*     Mass scheme
*
*     Fixed Flavour Number Scheme
*
      wnl = 3
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
            Sg0(1,alpha) = f0ev(1,alpha)
            Sg0(2,alpha) = f0ev(2,alpha)
            Sg0(3,alpha) = f0ev(3,alpha)
            Duc0(alpha)  = f0ev(4,alpha)
            Dds0(alpha)  = f0ev(5,alpha)
            Dsb0(alpha)  = f0ev(6,alpha)
            Dct0(alpha)  = f0ev(7,alpha)
            um0(alpha)   = f0ev(8,alpha)
            dm0(alpha)   = f0ev(9,alpha)
            sm0(alpha)   = f0ev(10,alpha)
            cm0(alpha)   = f0ev(11,alpha)
            bm0(alpha)   = f0ev(12,alpha)
            tm0(alpha)   = f0ev(13,alpha)
         enddo
*     Singlet
         call odeintsgQEDf(muF20,muF2,Sg0,Sg)
*     Non-Singlet
         call odeintnsQEDf(2,muF20,muF2,Dds0,Dds)
         call odeintnsQEDf(1,muF20,muF2,um0,um)
         call odeintnsQEDf(2,muF20,muF2,dm0,dm)
         call odeintnsQEDf(2,muF20,muF2,sm0,sm)
         if(Nf_FF.eq.3)then
            do alpha=0,nin(igrid)
               Duc(alpha) = ( Sg(2,alpha) + Sg(3,alpha) ) / 2d0
               Dsb(alpha) = ( Sg(2,alpha) - Sg(3,alpha) 
     1                    - 2d0 * Dds(alpha) ) / 4d0
               Dct(alpha) = 0d0
               cm(alpha)  = 0d0
               bm(alpha)  = 0d0
               tm(alpha)  = 0d0
            enddo
         endif
         if(Nf_FF.eq.4)then
            call odeintnsQEDf(1,muF20,muF2,Duc0,Duc)
            call odeintnsQEDf(1,muF20,muF2,cm0,cm)
            do alpha=0,nin(igrid)
               Dsb(alpha) = ( Sg(2,alpha) - Sg(3,alpha) 
     1                    - 2d0 * Dds(alpha) ) / 4d0
               Dct(alpha) = ( Sg(2,alpha) + Sg(3,alpha) 
     1                    - 2d0 * Duc(alpha) ) / 4d0
               bm(alpha)  = 0d0
               tm(alpha)  = 0d0
            enddo
         endif
         if(Nf_FF.eq.5)then
            call odeintnsQEDf(1,muF20,muF2,Duc0,Duc)
            call odeintnsQEDf(2,muF20,muF2,Dsb0,Dsb)
            call odeintnsQEDf(1,muF20,muF2,cm0,cm)
            call odeintnsQEDf(2,muF20,muF2,bm0,bm)
            do alpha=0,nin(igrid)
               Dct(alpha) = ( Sg(2,alpha) + Sg(3,alpha) 
     1                    - 2d0 * Duc(alpha) ) / 4d0
               tm(alpha)  = 0d0
            enddo
         endif
         if(Nf_FF.eq.6)then
            call odeintnsQEDf(1,muF20,muF2,Duc0,Duc)
            call odeintnsQEDf(2,muF20,muF2,Dsb0,Dsb)
            call odeintnsQEDf(1,muF20,muF2,Dct0,Dct)
            call odeintnsQEDf(1,muF20,muF2,cm0,cm)
            call odeintnsQEDf(2,muF20,muF2,bm0,bm)
            call odeintnsQEDf(1,muF20,muF2,tm0,tm)
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
               mu2f(inf) = m2th(inf+1)
            enddo
         elseif(sgn.eq.-1)then
            do inf=nfi-1,nff,sgn
               mu2i(inf) = m2th(inf+1)
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
            do alpha=0,nin(igrid)
               Sg0(1,alpha) = f0ev(1,alpha)
               Sg0(2,alpha) = f0ev(2,alpha)
               Sg0(3,alpha) = f0ev(3,alpha)
               Duc0(alpha)  = f0ev(4,alpha)
               Dds0(alpha)  = f0ev(5,alpha)
               Dsb0(alpha)  = f0ev(6,alpha)
               Dct0(alpha)  = f0ev(7,alpha)
               um0(alpha)   = f0ev(8,alpha)
               dm0(alpha)   = f0ev(9,alpha)
               sm0(alpha)   = f0ev(10,alpha)
               cm0(alpha)   = f0ev(11,alpha)
               bm0(alpha)   = f0ev(12,alpha)
               tm0(alpha)   = f0ev(13,alpha)
            enddo
*     Singlet
            call odeintsgQEDf(mu2i(inf),mu2f(inf),Sg0,Sg)
*     Non-Singlet
            call odeintnsQEDf(2,mu2i(inf),mu2f(inf),Dds0,Dds)
            call odeintnsQEDf(1,mu2i(inf),mu2f(inf),um0,um)
            call odeintnsQEDf(2,mu2i(inf),mu2f(inf),dm0,dm)
            call odeintnsQEDf(2,mu2i(inf),mu2f(inf),sm0,sm)
            if(inf.eq.3)then
               do alpha=0,nin(igrid)
                  Duc(alpha) = ( Sg(2,alpha) + Sg(3,alpha) ) / 2d0
                  Dsb(alpha) = ( Sg(2,alpha) - Sg(3,alpha) 
     1                       - 2d0 * Dds(alpha) ) / 4d0
                  Dct(alpha) = 0d0
                  cm(alpha)  = 0d0
                  bm(alpha)  = 0d0
                  tm(alpha)  = 0d0
               enddo
            endif
            if(inf.eq.4)then
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),Duc0,Duc)
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),cm0,cm)
               do alpha=0,nin(igrid)
                  Dsb(alpha) = ( Sg(2,alpha) - Sg(3,alpha) 
     1                       - 2d0 * Dds(alpha) ) / 4d0
                  Dct(alpha) = ( Sg(2,alpha) + Sg(3,alpha) 
     1                       - 2d0 * Duc(alpha) ) / 4d0
                  bm(alpha)  = 0d0
                  tm(alpha)  = 0d0
               enddo
            endif
            if(inf.eq.5)then
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),Duc0,Duc)
               call odeintnsQEDf(2,mu2i(inf),mu2f(inf),Dsb0,Dsb)
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),cm0,cm)
               call odeintnsQEDf(2,mu2i(inf),mu2f(inf),bm0,bm)
               do alpha=0,nin(igrid)
                  Dct(alpha) = ( Sg(2,alpha) + Sg(3,alpha) 
     1                       - 2d0 * Duc(alpha) ) / 4d0
                  tm(alpha)  = 0d0
               enddo
            endif
            if(inf.eq.6)then
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),Duc0,Duc)
               call odeintnsQEDf(2,mu2i(inf),mu2f(inf),Dsb0,Dsb)
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),Dct0,Dct)
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),cm0,cm)
               call odeintnsQEDf(2,mu2i(inf),mu2f(inf),bm0,bm)
               call odeintnsQEDf(1,mu2i(inf),mu2f(inf),tm0,tm)
            endif
*
*     Update initial coditions for the next step
*
            do alpha=0,nin(igrid)
               f0ev(1,alpha)  = Sg(1,alpha)
               f0ev(2,alpha)  = Sg(2,alpha)
               f0ev(3,alpha)  = Sg(3,alpha)
               f0ev(4,alpha)  = Duc(alpha)
               f0ev(5,alpha)  = Dds(alpha)
               f0ev(6,alpha)  = Dsb(alpha)
               f0ev(7,alpha)  = Dct(alpha)
               f0ev(8,alpha)  = um(alpha)
               f0ev(9,alpha)  = dm(alpha)
               f0ev(10,alpha) = sm(alpha)
               f0ev(11,alpha) = cm(alpha)
               f0ev(12,alpha) = bm(alpha)
               f0ev(13,alpha) = tm(alpha)
            enddo
         enddo
      endif
*
*     Put evolved PDFs in a single vector
*
      do alpha=0,nin(igrid)
         fev(1,alpha)  = Sg(1,alpha)
         fev(2,alpha)  = Sg(2,alpha)
         fev(3,alpha)  = Sg(3,alpha)
         fev(4,alpha)  = Duc(alpha)
         fev(5,alpha)  = Dds(alpha)
         fev(6,alpha)  = Dsb(alpha)
         fev(7,alpha)  = Dct(alpha)
         fev(8,alpha)  = um(alpha)
         fev(9,alpha)  = dm(alpha)
         fev(10,alpha) = sm(alpha)
         fev(11,alpha) = cm(alpha)
         fev(12,alpha) = bm(alpha)
         fev(13,alpha) = tm(alpha)
      enddo
*
*     Rotate PDFs back into the physical basis
*
 104  call PDFevQED2phys(fev,fpy)
*
*     Put evolved PDFs in the physical basis in the common used for the interpolation
*
      do alpha=0,nin(igrid)
         fph(igrid,0,alpha)     = f0lep(0,alpha) 
         do i=1,6
            fph(igrid,i,alpha)  = fpy(i,alpha)
            fph(igrid,-i,alpha) = fpy(-i,alpha)
         enddo
         fgamma(igrid,alpha)    = fpy(0,alpha)
         flepton(igrid,0,alpha) = fpy(0,alpha)
         do i=1,3
            flepton(igrid,i,alpha)  = 0d0
            flepton(igrid,-i,alpha) = 0d0
         enddo
      enddo
*
      return
      end
