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
      double precision f0(-6:6),fp0,fext0(-6:7),flext0(-3:3)
      logical has_photon
      external ExternalSetAPFEL
      external ExternalSetAPFEL1
      external ExternalSetAPFELLept
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
            do ilept=-3,3
               f0lep(ilept,alpha) = flepton(igrid,ilept,alpha)
               if(dabs(f0lep(ilept,alpha)).lt.1d-14)
     1              f0lep(ilept,alpha) = 0d0
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
         if(pdfset(9:9).eq."1")then
            do alpha=0,nin(igrid)
               call ExternalSetAPFEL1(xg(igrid,alpha),dsqrt(Q20),fext0)
               do ifl=-6,6
                  f0ph(ifl,alpha) = fext0(ifl)
               enddo
               f0lep(0,alpha) = fext0(7)
               do ilept=1,3
                  f0lep(ilept,alpha) = 0d0
                  f0lep(-ilept,alpha) = 0d0
               enddo
            enddo
         else
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
         endif
*
*     External Set with leptons
*
      elseif(pdfset(1:12).eq."leptexternal")then
         do alpha=0,nin(igrid)
            call ExternalSetAPFELLept(xg(igrid,alpha),dsqrt(Q20),
     1                                flext0,fext0)
            do ifl=-6,6
               f0ph(ifl,alpha) = fext0(ifl)
            enddo
            do ilept=-3,3
               f0lep(ilept,alpha) = flext0(ilept)
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
               f0lep(0,alpha) = fp0
               do ilept=1,3
                  f0lep(ilept,alpha) = 0d0
                  f0lep(-ilept,alpha) = 0d0
               enddo
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
