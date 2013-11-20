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
      include "../commons/Replica.h"
**
*     Input Variables
*
      double precision Q20
**
*     Internal Variables
*
      integer alpha
      integer ifl
      integer Nrep
      double precision f0(-6:6),fp0
      logical has_photon
*
*     User defined PDFs
*
      if(pdfset(1:7).eq."private")then
         do alpha=0,nin(igrid)
            call private(xg(igrid,alpha),f0)
            do ifl=-6,6
               f0ph(ifl,alpha)  = f0(ifl)
            enddo
            f0bos(alpha) = 0d0
         enddo
*
*     LHA Toy PDFs
*
      elseif(pdfset(1:5).eq."ToyLH")then
         do alpha=0,nin(igrid)
            call toyLHPDFs(xg(igrid,alpha),f0)
            do ifl=-6,6
               f0ph(ifl,alpha)  = f0(ifl)
            enddo
            f0bos(alpha) = 0d0
         enddo
*
*     LHAPDF set
*
      else
         if(igrid.eq.1)then
            call InitPDFsetbyName(pdfset)
            call numberPDF(Nrep)
            if(irep.lt.0.or.irep.gt.Nrep)then
               write(6,*) "Replica requested ot of range:"
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
                  f0ph(ifl,alpha)  = f0(ifl)
               enddo
               f0bos(alpha) = fp0
            enddo
         else
            do alpha=0,nin(igrid)
               call evolvePDF(xg(igrid,alpha),dsqrt(Q20),f0)
               fp0 = 0d0
               do ifl=-6,6
                  f0ph(ifl,alpha)  = f0(ifl)
               enddo
               f0bos(alpha) = fp0
            enddo
         endif
      endif
*
      return
      end
