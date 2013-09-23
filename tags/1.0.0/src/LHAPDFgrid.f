************************************************************************
*
*     LHAPDFgrid.f:
*
*     This subrotine produces the file *.LHgrid that can be used with 
*     LHAPDF.
*
************************************************************************
      subroutine LHAPDFgrid(Nrep,Qin,fname)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/m2th.h"
      include "../commons/kren.h"
      include "../commons/Th.h"
      include "../commons/LHAgrid.h"
**
*     Input Variables
*
      integer Nrep
      double precision Qin
      character*30 fname
**
*     Intenal Variables
*
      integer ln
      integer i,ix,iq2,ipdf,krep
      double precision qref,mth(4:6)
      double precision lref,lambda3,lambda4,lambda5,lambdaNF
      double precision xbLHA(nxLHA),q2LHA(nq2LHA)
      double precision xpdfLHA(-6:6,nxLHA,nq2LHA)
      double precision xgammaLHA(nxLHA,nq2LHA)
      double precision xPDF,xgamma
      parameter(lref=0.239d0)
*
*     Specify initialization for the LHgrid
*
      call SetQLimits(dsqrt(q2minLHA),dsqrt(q2maxLHA))
      call SetNumberOfGrids(3)
      call SetGridParameters(1,100,3,xminLHA)
      call SetGridParameters(2,50,5,xmLHA)
      call SetGridParameters(3,20,5,8d-1)
*
      call initializeAPFEL
*
*     Compute the values of lambdaQCD
*
c      lambda5 = lref
c      call lambda(4,lref,lambda4)
c      call lambda(3,lref,lambda3)
      lambda5 = lambdaNF(5,alpha_ref_qcd,q2_ref_qcd)
      lambda4 = lambdaNF(4,alpha_ref_qcd,q2_ref_qcd)
      lambda3 = lambdaNF(3,alpha_ref_qcd,q2_ref_qcd)
*
*     Open the output file
*
      ln = index(fname,char(0)) - 1
      open(unit=13,status="unknown",file=fname(1:ln)//".LHgrid")
*
*     Write header
*
      write(13,*) "'Version'"," '5.8'"
      write(13,*) "'Description:'"
      write(13,*) "'",fname(1:ln),".LHgrid'"
      write(13,*) "'Reference: APFEL PDFs'"
      write(13,*) "'APFEL Collaboration: '"
      write(13,*) "'V. Bertone, S. Carrazza'"
      write(13,*)"' '"
      write(13,*)"' '"
      write(13,*)"'arXiv:yymm.xxxx'"
      write(13,*)"'This set has ",Nrep+1," member PDF'"
      write(13,*)"'mem=0 --> average on replicas. '"
      write(13,*)"'Evolution: pol. interpolation on the LH grid.'"
*
      write(13,*) "'Alphas:'"
      if(ipt.eq.0)then
         write(13,*) "'Variable',","'lo',","'EvolCode'"
      elseif(ipt.eq.1)then
         write(13,*) "'Variable',","'nlo',","'EvolCode'"
      elseif(ipt.eq.2)then
         write(13,*) "'Variable',","'nnlo',","'EvolCode'"
      endif
*
      qref = dsqrt(q2_ref_qcd)
      do i=3,6
         mth(i) = dsqrt(m2th(i))
      enddo
*
      write(13,*) 1,",",qref,",",mth(4),",",mth(5),",",mth(6)
*
      write(13,*) "'MinMax:'"
      write(13,*) Nrep,",",1
      write(13,*) xminLHA,",",xmaxLHA,",",q2minLHA,",",q2maxLHA
*
      write(13,*) "'QCDparams:'"
      write(13,*) Nrep,",",1
      write(13,*) lambda4,",",lambda5
*
      write(13,*) "'Parameterlist:'"
      write(13,*) "'list',",Nrep,",",1
*
      do i=0,Nrep
         write(13,*) alpha_ref_qcd
      enddo
*
      write(13,* ) "'Evolution:'"
      if(ipt.eq.0)then
         write(13,*) "'lo',",q2minLHA,",",kren
      elseif(ipt.eq.1)then
         write(13,*) "'nlo',",q2minLHA,",",kren
      elseif(ipt.eq.2)then
         write(13,*) "'nnlo',",q2minLHA,",",kren
      endif
      if(Th.eq."QCD")then
         write(13,*) "'NNPDF20int'"
      else
         write(13,*) "'NNPDF20intqed'"
      endif
      write(13,*) Nrep,",",1
*
*     Compute and write the grids in x and Q2
*
      write(13,*) nxLHA
      do ix=1,nxLHA
         if(ix.le.nxmLHA)then
            xbLHA(ix) = xminLHA * ( xmLHA / xminLHA )
     1                **( 2d0 * dble( ix-1 ) / dble( nxLHA - 1 ) )
         else
            xbLHA(ix) = xmLHA + ( xmaxLHA - xmLHA )
     1                * ( dble( ix - nxmLHA - 1 )
     2                / dble( nxLHA - nxmLHA - 1 ) )
         endif
         write(13,*) xbLHA(ix)
      enddo
*
      write(13,*) nq2LHA
      write(13,*) q2minLHA
      do iq2=1,nq2LHA
         q2LHA(iq2) = q2minLHA * ( q2maxLHA / q2minLHA )
     1              **( dble( iq2 - 1 ) / dble( nq2LHA - 1 ) )
         write(13,*) q2LHA(iq2)
      enddo
*
*     EvolvePPDFs starting from Qin and write them on file
*
      write(13,*) Nrep
      do krep=0,Nrep
         write(6,*) "Evaluating replica",krep
         call SetReplica(krep)
         do iq2=1,nq2LHA
            call EvolveAPFEL(Qin,dsqrt(q2LHA(iq2)))
            do ix=1,nxLHA
               do ipdf=-6,6
                  xpdfLHA(ipdf,ix,iq2) = xPDF(ipdf,xbLHA(ix))
               enddo
               xgammaLHA(ix,iq2) = xgamma(xbLHA(ix))
            enddo
         enddo
*
         if(Th.eq."QCD")then
            do ix=1,nxLHA
               do iq2=1,nq2LHA
                  write(13,44) (xpdfLHA(ipdf,ix,iq2),ipdf=-6,6)
               enddo
            enddo
         else
            do ix=1,nxLHA
               do iq2=1,nq2LHA
                  write(13,45) (xpdfLHA(ipdf,ix,iq2),ipdf=-6,6),
     1                          xgammaLHA(ix,iq2)
               enddo
            enddo
         endif
      enddo
* 
      write(13,* ) "'End:'"
      close(13)
*
      write(6,*) "File ",fname(1:ln),".LHgrid produced!"
      write(6,*) "  "
*
 44   format(1x,13((es14.7),1x))
 45   format(1x,14((es14.7),1x))
*
      return
      end
