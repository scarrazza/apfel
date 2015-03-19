************************************************************************
*
*     initPDFs.f:
*
*     This routine initializes the PDFs at the initial scale Q20
*     on the interpolation grid.
*
************************************************************************
      subroutine initPDFs(Q20)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/pdfset.h"
      include "../commons/f0ph.h"
      include "../commons/fph.h"
      include "../commons/Replica.h"
**
*     Input Variables
*
      double precision Q20
**
*     Internal Variables
*
      integer alpha
      integer ifl,ilept
      integer Nrep
      double precision f0(-6:6),fp0,fext0(-6:7)
      logical has_photon
      external ExternalSetAPFEL
*     Fragmentation functions variables
      integer iset,icharge
      integer ih,ic,io
      double precision ff(-5:5),grad(-5:5,17)
*     Common needed by the DSS fragmentations routine
      double precision fini
      common / fragini / fini
*
*     User defined PDFs
*
      if(pdfset(1:7).eq."private")then
         do alpha=0,nin(igrid)
            call private(xg(igrid,alpha),f0)
            do ifl=-6,6
               f0ph(ifl,alpha) = f0(ifl)
            enddo
            do ilept=-3,3
               f0lep(ilept,alpha) = 0d0
            enddo
         enddo
*
*     In case one wants to use a PDF set previously evolved by APFEL
*     as an input.
*
      elseif(pdfset(1:5).eq."apfel")then
         do alpha=0,nin(igrid)
            do ifl=-6,6
               f0ph(ifl,alpha) = fph(igrid,ifl,alpha)
               if(dabs(f0ph(ifl,alpha)).lt.1d-14) f0ph(ifl,alpha) = 0d0
            enddo
            f0lep(0,alpha) = fgamma(igrid,alpha)
            if(dabs(f0lep(0,alpha)).lt.1d-14) f0lep(0,alpha) = 0d0
            do ilept=1,3
               f0lep(ilept,alpha) = 0d0
               f0lep(-ilept,alpha) = 0d0
            enddo
         enddo
*
*     LHA Toy PDFs
*
      elseif(pdfset(1:5).eq."ToyLH")then
         do alpha=0,nin(igrid)
            call toyLHPDFs(xg(igrid,alpha),f0)
            do ifl=-6,6
               f0ph(ifl,alpha) = f0(ifl)
            enddo
            do ilept=-3,3
               f0lep(ilept,alpha) = 0d0
            enddo
         enddo
*
*     External Set
*
      elseif(pdfset(1:8).eq."external")then
         do alpha=0,nin(igrid)
            call ExternalSetAPFEL(xg(igrid,alpha),dsqrt(Q20),fext0)
            do ifl=-6,6
               f0ph(ifl,alpha) = fext0(ifl)
            enddo
            f0lep(0,alpha) = fext0(7)
            do ilept=1,3
               f0lep(ilept,alpha) = 0d0
               f0lep(-ilept,alpha) = 0d0
            enddo
         enddo
*
*     HKNS07 fragmentation functions
*
      elseif(pdfset(1:8).eq."hknsff07")then
*     pi^+ LO
         if(pdfset(10:15).eq."pip_lo")then
            iset    = 1
            icharge = 1
*     pi^+ NLO
         elseif(pdfset(10:16).eq."pip_nlo")then
            iset    = 2
            icharge = 1
*     pi^- LO
         elseif(pdfset(10:15).eq."pim_lo")then
            iset    = 1
            icharge = 2
*     pi^- NLO
         elseif(pdfset(10:16).eq."pim_nlo")then
            iset    = 2
            icharge = 2
*     K^+ LO
         elseif(pdfset(10:14).eq."Kp_lo")then
            iset    = 3
            icharge = 1
*     K^+ NLO
         elseif(pdfset(10:15).eq."Kp_nlo")then
            iset    = 4
            icharge = 1
*     K^- LO
         elseif(pdfset(10:14).eq."Km_lo")then
            iset    = 3
            icharge = 2
*     K^- NLO
         elseif(pdfset(10:15).eq."Km_nlo")then
            iset    = 4
            icharge = 2
*     proton LO
         elseif(pdfset(10:13).eq."p_lo")then
            iset    = 5
            icharge = 1
*     proton NLO
         elseif(pdfset(10:14).eq."p_nlo")then
            iset    = 6
            icharge = 1
*     antiproton LO
         elseif(pdfset(10:14).eq."pb_lo")then
            iset    = 5
            icharge = 2
*     antiproton NLO
         elseif(pdfset(10:15).eq."pb_nlo")then
            iset    = 6
            icharge = 2
         endif
*
         do alpha=0,nin(igrid)
            call HKNSFF(Q20,xg(igrid,alpha),iset,icharge,ff,grad)
            do ifl=-5,5
               f0ph(ifl,alpha) = xg(igrid,alpha) * ff(ifl)
            enddo
            f0ph(6,alpha)  = 0d0
            f0ph(-6,alpha) = 0d0
            do ilept=-3,3
               f0lep(ilept,alpha) = 0d0
            enddo
         enddo
*
*     DSS fragmentation functions
*
      elseif(pdfset(1:3).eq."dss")then
         fini = 0d0
*     pi^+ LO
         if(pdfset(5:10).eq."pip_lo")then
            ih = 1
            ic = 1
            io = 0
*     pi^+ NLO
         elseif(pdfset(5:11).eq."pip_nlo")then
            ih = 1
            ic = 1
            io = 1
*     pi^- LO
         elseif(pdfset(5:10).eq."pim_lo")then
            ih = 1
            ic = -1
            io = 0
*     pi^- NLO
         elseif(pdfset(5:11).eq."pim_nlo")then
            ih = 1
            ic = -1
            io = 1
*     K^+ LO
         elseif(pdfset(5:9).eq."Kp_lo")then
            ih = 2
            ic = 1
            io = 0
*     K^+ NLO
         elseif(pdfset(5:10).eq."Kp_nlo")then
            ih = 2
            ic = 1
            io = 1
*     K^- LO
         elseif(pdfset(5:9).eq."Km_lo")then
            ih = 2
            ic = -1
            io = 0
*     K^- NLO
         elseif(pdfset(5:10).eq."Km_nlo")then
            ih = 2
            ic = -1
            io = 1
*     proton LO
         elseif(pdfset(5:8).eq."p_lo")then
            ih = 3
            ic = 1
            io = 0
*     proton NLO
         elseif(pdfset(5:9).eq."p_nlo")then
            ih = 3
            ic = 1
            io = 1
*     antiproton LO
         elseif(pdfset(5:9).eq."pb_lo")then
            ih = 3
            ic = -1
            io = 0
*     antiproton NLO
         elseif(pdfset(5:10).eq."pb_nlo")then
            ih = 3
            ic = -1
            io = 1
*     hadron LO
         elseif(pdfset(5:8).eq."h_lo")then
            ih = 4
            ic = 1
            io = 0
*     hadron NLO
         elseif(pdfset(5:9).eq."h_nlo")then
            ih = 4
            ic = 1
            io = 1
*     antihadron LO
         elseif(pdfset(5:9).eq."hb_lo")then
            ih = 4
            ic = -1
            io = 0
*     antihadron NLO
         elseif(pdfset(5:10).eq."hb_nlo")then
            ih = 4
            ic = -1
            io = 1
         endif
*
         do alpha=0,nin(igrid)
            call fDSS (ih,ic,io,xg(igrid,alpha),Q20,ff(2),ff(-2),ff(1),
     1                 ff(-1),ff(3),ff(-3),ff(4),ff(5),ff(0))
            ff(-4) = ff(4)
            ff(-5) = ff(5)
            do ifl=-5,5
               f0ph(ifl,alpha) = ff(ifl)
            enddo
            f0ph(6,alpha)  = 0d0
            f0ph(-6,alpha) = 0d0
            do ilept=-3,3
               f0lep(ilept,alpha) = 0d0
            enddo
         enddo
*
*     Kretzer's parametrization at Q2 = 0.4 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0003177.
*
      elseif(pdfset(1:7).eq."kretzer")then
         do alpha=0,nin(igrid)
            call KretzerFFs(xg(igrid,alpha),f0)
            do ifl=-6,6
               f0ph(ifl,alpha) = f0(ifl)
            enddo
            do ilept=-3,3
               f0lep(ilept,alpha) = 0d0
            enddo
         enddo
*
*     HKNS parametrization at Q2 = 1 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0702250 (used for the benchmark against MELA).
*
      elseif(pdfset(1:4).eq."MELA")then
         do alpha=0,nin(igrid)
            call HKNSFFs(xg(igrid,alpha),f0)
            do ifl=-6,6
               f0ph(ifl,alpha) = f0(ifl)
            enddo
            do ilept=-3,3
               f0lep(ilept,alpha) = 0d0
            enddo
         enddo
*
*     LHAPDF set
*
      else
         if(igrid.eq.1)then
            call InitPDFsetbyname(pdfset)
            call numberPDF(Nrep)
            if(irep.lt.0.or.irep.gt.Nrep)then
               write(6,*) "Replica requested out of range:"
               write(6,*) "- irep=",irep
               write(6,*) "- Nrep=",Nrep
               call exit(-10)
            endif
         endif
         call InitPDF(irep)
         if(has_photon())then
            do alpha=0,nin(igrid)
               call evolvePDFphoton(xg(igrid,alpha),dsqrt(Q20),f0,fp0)
               do ifl=-6,6
                  f0ph(ifl,alpha) = f0(ifl)
               enddo
               do ilept=1,3
                  f0lep(ilept,alpha) = 0d0
                  f0lep(-ilept,alpha) = 0d0
               enddo
               f0lep(0,alpha) = fp0
            enddo
         else
            do alpha=0,nin(igrid)
               call evolvePDF(xg(igrid,alpha),dsqrt(Q20),f0)
               fp0 = 0d0
               do ifl=-6,6
                  f0ph(ifl,alpha) = f0(ifl)
               enddo
               do ilept=-3,3
                  f0lep(ilept,alpha) = 0d0
               enddo
            enddo
         endif
      endif
*
      return
      end
