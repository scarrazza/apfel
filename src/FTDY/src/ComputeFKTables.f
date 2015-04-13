************************************************************************
*
*     ComputeFKTables.f:
*
*     Put in the FK form and split into sets the evolution tables
*     produced by the project EvolutionKernels for Drell-Yan.
*
************************************************************************
      subroutine ComputeFKTables(inputfile,Q0,flmap)
*
      implicit none
*
      include "../commons/mxdata.h"
      include "../commons/kinematics.h"
      include "../commons/set.h"
**
*     Input variables
*
      double precision Q0
      character*(*) inputfile
**
*     Internal Variables
*
      integer idat
      integer ounit
      integer ls
      integer ipdf,jpdf
      double precision t1,t2
**
*     Output Variables
*
      integer flmap(0:13,0:13)
*
*     Initialize APFEL
*
      call EnableWelcomeMessage(.false.)
      call SetFastEvolution(.false.)
      call EnableEvolutionOperator(.true.)
      call LockGrids(.true.)
      call initParametersDIS()
      call InitializeAPFEL
*
*     Initialize Drell-Yan couplings
*
      call initDYcouplings
*
*     Read partonic cross sections
*
      call readcDY(inputfile)
*
      if(set.eq."DYE886R")then
         set = "DYE886R_P"
      elseif(set.eq."CDFWASYM")then
         set = "CDFWASYM_WP"
      endif
*
*     Write preamble on unit ounit (I still need to settle the k-factors)
*
 103  call cpu_time(t1)
*
      ounit = 16
      ls = index(set," ") - 1
      open(unit=ounit,file="FK_"//set(1:ls)//".dat",status="unknown")
c      open(unit=ounit,file="FK_"//set(1:ls)//".dat",status="old",
c     1     access="append")
      write(6,*) "Processing ",set(1:ls)," set ..."
      do idat=1,ndata
         write(6,"(a,i4,a,i4)") " Convoluting data number =",idat,
     1                          " /",ndata
         if(set.eq."DYE886R_P")then
            obs(idat) = "DYP_E886P"
         elseif(set.eq."DYE886R_D")then
            obs(idat) = "DYP_E886D"
         elseif(set.eq."CDFWASYM_WP")then
            obs(idat) = "EWK_WASYM_WP"
         elseif(set.eq."CDFWASYM_WM")then
            obs(idat) = "EWK_WASYM_WM"
         endif
*
*     Evaluate sigma tables
*
         call sigmafk_dy(idat,Q0)
*
*     Write FK tables and compute predictions
*
         call writeFK(idat,ounit,flmap)
      enddo
      close(ounit)
*
      write(6,*) "Flavour map:"
      do ipdf=0,13
         write(6,"(14(i2,x))") (flmap(ipdf,jpdf),jpdf=0,13)
      enddo
      write(6,*) " "
*
      if(set.eq."DYE886R_P")then
         set = "DYE886R_D"
         goto 103
      elseif(set.eq."CDFWASYM_WP")then
         set = "CDFWASYM_WM"
         goto 103
      endif
*
      call cpu_time(t2)
*
      write(6,*) "Time taken for the FK table generation =",t2-t1," s"
      write(6,*) "   "
*
      return
      end
