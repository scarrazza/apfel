************************************************************************
*
*     LHAPDFgridDerivative.f:
*
*     If LHAPDF5 is found, this subrotine produces the *.LHgrid file 
*     that can be used with LHAPDF v5.9. If instead LHAPDF6 is found
*     a folder contained the PDF set in a suitable format is produced.
*
************************************************************************
      subroutine LHAPDFgridDerivative(Nrep,fname)
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
      logical islhapdf6
      character*30 fname
**
*     Intenal Variables
*
      integer ln
      integer i,ix,iq2,ipdf,krep
      double precision qref,mth(4:6),alphasPDF
      double precision lref,lambda3,lambda4,lambda5,lambdaNF
      double precision xbLHA(nxLHA),q2LHA(nq2LHA)
      double precision dxpdfLHA(-6:6,nxLHA,nq2LHA)
      double precision dxgammaLHA(nxLHA,nq2LHA)
      double precision dxPDF,dxgamma
      double precision AlphaQCD
      character*30 str
      parameter(lref=0.239d0)
*
*     For the moment only if the fast evolution is off
*
      call SetFastEvolution(.false.)
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
*     LHAPDF5 output
*
      if (islhapdf6().eqv..false.) then
*
*     Open the output file
*
         ln = index(fname,char(0)) - 1
         open(unit=13,status="unknown",file=fname(1:ln)//".LHgrid")
*
*     Write header
*
         write(13,*) "'Version'"," '5.9'"
         write(13,*) "'Description:'"
         write(13,*) "'",fname(1:ln),".LHgrid'"
         write(13,*) "'Reference: arXiv:1310.1394'"
         write(13,*) "'APFEL: A PDF Evolution Library'"
         write(13,*) "'V. Bertone, S. Carrazza, J. Rojo'"
         write(13,*)"' ***********************'"
         write(13,*)"' *** PDF DERIVATIVES ***'"
         write(13,*)"' ***********************'"
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
c     qref = dsqrt(q2_ref_qcd)
         qref = 91.2d0          ! Z mass
         do i=4,6
            call GetQmass(i,mth(i))
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
c     write(13,*) alpha_ref_qcd
            write(13,*) alphasPDF(qref)
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
     1              **( 2d0 * dble( ix-1 ) / dble( nxLHA - 1 ) )
            else
               xbLHA(ix) = xmLHA + ( xmaxLHA - xmLHA )
     1              * ( dble( ix - nxmLHA - 1 )
     2              / dble( nxLHA - nxmLHA - 1 ) )
            endif
            write(13,*) xbLHA(ix)
         enddo
*     
         write(13,*) nq2LHA
         write(13,*) q2minLHA
         do iq2=1,nq2LHA
            q2LHA(iq2) = q2minLHA * ( q2maxLHA / q2minLHA )
     1           **( dble( iq2 - 1 ) / dble( nq2LHA - 1 ) )
            write(13,*) q2LHA(iq2)
         enddo
*     
*     EvolvePPDFs starting from Qin and write them on file
*     
         write(13,*) Nrep
         do krep=0,Nrep
            write(6,*) "Evaluating replica",krep," ..."
            call SetReplica(krep)
            do iq2=1,nq2LHA
               call DeriveAPFEL(dsqrt(q2LHA(iq2)))
               do ix=1,nxLHA
                  do ipdf=-6,6
                     dxpdfLHA(ipdf,ix,iq2) = dxPDF(ipdf,xbLHA(ix))
                  enddo
                  dxgammaLHA(ix,iq2) = dxgamma(xbLHA(ix))
               enddo
            enddo
*     
            if(Th.eq."QCD")then
               do ix=1,nxLHA
                  do iq2=1,nq2LHA
                     write(13,44) (dxpdfLHA(ipdf,ix,iq2),ipdf=-6,6)
                  enddo
               enddo
            else
               do ix=1,nxLHA
                  do iq2=1,nq2LHA
                     write(13,45) (dxpdfLHA(ipdf,ix,iq2),ipdf=-6,6),
     1                    dxgammaLHA(ix,iq2)
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
      else
*
*     LHAPDF6 output
*
         ln = index(fname,char(0)) - 1
*     creating main folder
         call mkdir(fname(1:ln))
*     creating info file
         open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1        //fname(1:ln)//".info")
*     
*     Write header
*     
         write(13,*) "SetDesc: ",fname(1:ln)
         write(13,*) "Authors: APFEL a PDF Evolution Library"
         write(13,*) "Reference: 1310.1394"
         write(13,*) "Format: lhagrid1"
         write(13,*) "DataVersion: 1"
         write(13,*) "NumMembers:", Nrep+1
         write(13,*) "Particle: 2212"
         if(Th.eq."QCD") then
            write(13,*) "Flavors: [-6, -5, -4, -3, -1, -2,",
     1           "21, 2, 1, 3, 4, 5, 6]"
         else
            write(13,*) "Flavors: [-6, -5, -4, -3, -1, -2, 21, 2, 1,", 
     1           "3, 4, 5, 6, 22]"
         endif
         write(13,*) "OrderQCD:",ipt
         write(13,*) "FlavorScheme: variable"
         write(13,*) "NumFlavors: 6"
         write(13,*) "ErrorType: replicas"
         write(13,*) "XMin:", xminLHA
         write(13,*) "XMax:", xmaxLHA
         write(13,*) "QMin:", dsqrt(q2minLHA)
         write(13,*) "QMax:", dsqrt(q2maxLHA)
         write(13,*) "MZ: 91.1876"
         write(13,*) "MUp: 0"
         write(13,*) "MDown: 0"
         write(13,*) "MStrange: 0"
         write(13,*) "MCharm:",dsqrt(m2th(4))
         write(13,*) "MBottom:",dsqrt(m2th(5))
         write(13,*) "MTop:",dsqrt(m2th(6))
         write(13,*) "AlphaS_MZ:",alpha_ref_qcd
         write(13,*) "AlphaS_OrderQCD:",ipt
         write(13,*) "AlphaS_Type: ipol"
         
         do ix=1,nxLHA
            if(ix.le.nxmLHA)then
               xbLHA(ix) = xminLHA * ( xmLHA / xminLHA )
     1              **( 2d0 * dble( ix-1 ) / dble( nxLHA - 1 ) )
            else
               xbLHA(ix) = xmLHA + ( xmaxLHA - xmLHA )
     1              * ( dble( ix - nxmLHA - 1 )
     2              / dble( nxLHA - nxmLHA - 1 ) )
            endif
         enddo

         do iq2=1,nq2LHA
            q2LHA(iq2) = q2minLHA * ( q2maxLHA / q2minLHA )
     1           **( dble( iq2 - 1 ) / dble( nq2LHA - 1 ) )
         enddo

         write(13,*)"AlphaS_Qs: [",(dsqrt(q2LHA(iq2)),",",iq2=1,nq2LHA),
     1        "]" 
         write(13,*) "AlphaS_Vals: [", 
     1        (AlphaQCD(dsqrt(q2LHA(iq2))),",",iq2=1,nq2LHA), "]"
         write(13,*) "AlphaS_Lambda4:", lambda4
         write(13,*) "AlphaS_Lambda5:", lambda5
         
         close(13)
*
*     Now Loop over all replicas and print to file
*
         do krep=0,Nrep
            write(6,*) "Evaluating replica",krep," ..."            

            if (krep.lt.10) then
               write(str,'(i1)' ) krep
               open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1              //fname(1:ln)//"_000"//str(1:1)//".dat")
            elseif (krep.lt.100) then
               write(str,'(i2)' ) krep
               open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1              //fname(1:ln)//"_00"//str(1:2)//".dat")
            elseif (krep.lt.1000) then               
               write(str,'(i3)' ) krep
               open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1              //fname(1:ln)//"_0"//str(1:3)//".dat")
            else
               write(str,'(i4)' ) krep
               open(unit=13,status="unknown",file=fname(1:ln)//"/"
     1              //fname(1:ln)//"_"//str(1:4)//".dat")
            endif
            if (krep.eq.0) then
               write(13,*) "PdfType: central"
            else
               write(13,*) "PdfType: replica"
            endif
            write(13,*) "Format: lhagrid1"
            write(13,*) "---"
            write(13,*) (xbLHA(ix), ix=1,nxLHA)
            write(13,*) (dsqrt(q2LHA(iq2)), iq2=1,nq2LHA)
*
            if(Th.eq."QCD") then
               write(13,"(a)") "-6 -5 -4 -3 -1 -2 21 2 1 3 4 5 6"
            else
               write(13,"(a)") "-6 -5 -4 -3 -1 -2 21 2 1 3 4 5 6 22"
            endif
*
            call SetReplica(krep)
            do iq2=1,nq2LHA
               call DeriveAPFEL(dsqrt(q2LHA(iq2)))
               do ix=1,nxLHA
                  do ipdf=-6,6
                     dxpdfLHA(ipdf,ix,iq2) = dxPDF(ipdf,xbLHA(ix))
                  enddo
                  dxgammaLHA(ix,iq2) = dxgamma(xbLHA(ix))
               enddo
            enddo
*     
            if(Th.eq."QCD")then
               do ix=1,nxLHA
                  do iq2=1,nq2LHA
                     write(13,44) (dxpdfLHA(ipdf,ix,iq2),ipdf=-6,6)
                  enddo
               enddo
            else
               do ix=1,nxLHA
                  do iq2=1,nq2LHA
                     write(13,45) (dxpdfLHA(ipdf,ix,iq2),ipdf=-6,6),
     1                    dxgammaLHA(ix,iq2)
                  enddo
               enddo
            endif
            write(13,*) "---"
            close(13)
         enddo
*     
         write(6,*) "File ",fname(1:ln)," grid produced!"
         write(6,*) "  "     
      endif
*
 44   format(1x,13((es14.7),1x))
 45   format(1x,14((es14.7),1x))
*
      return
      end
