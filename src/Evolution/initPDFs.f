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
               f0ph(ifl,alpha) = f0(ifl)
            enddo
            f0bos(alpha) = 0d0
         enddo
*
*     "delta" PDFs for the PDFs independent evolution tables.
*
*     In this case the the index irep, which in the other cases
*     is associated to the PDF replica, is instead associated
*     to the second grid index.
*
      elseif(pdfset(1:5).eq."delta")then
c         if(irep.lt.0.or.irep.gt.nin(igrid))then
c            write(6,*) "Replica requested out of range:"
c            write(6,*) "- irep=",irep
c            write(6,*) "- Nrep=",nin(igrid)
c            call exit(-10)
c         endif
         do alpha=0,nin(igrid)
            do ifl=-6,6
               f0ph(ifl,alpha) = 0d0
            enddo
            f0bos(alpha) = 0d0
         enddo
         do ifl=-6,6
            f0ph(ifl,irep) = 1d0
         enddo
         f0bos(irep) = 1d0
*
*     In case one wants to use a set PDFs previously evolved by APFEL
*
      elseif(pdfset(1:7).eq."apfel")then
         do alpha=0,nin(igrid)
            do ifl=-6,6
               f0ph(ifl,alpha) = fph(igrid,ifl,alpha)
               if(dabs(f0ph(ifl,alpha)).lt.1d-14) f0ph(ifl,alpha) = 0d0
            enddo
            f0bos(alpha) = fgamma(igrid,alpha)
            if(dabs(f0bos(alpha)).lt.1d-14) f0bos(alpha) = 0d0
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
               f0bos(alpha) = fp0
            enddo
         else
            do alpha=0,nin(igrid)
               call evolvePDF(xg(igrid,alpha),dsqrt(Q20),f0)
               fp0 = 0d0
               do ifl=-6,6
                  f0ph(ifl,alpha) = f0(ifl)
               enddo
               f0bos(alpha) = fp0
            enddo
         endif
      endif
*
      return
      end
