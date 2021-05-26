************************************************************************
*
*     ComputeFKTables.f:
*
*     Put in the FK form and split into sets the evolution tables
*     produced by the project EvolutionKernels for Drell-Yan.
*
************************************************************************
      subroutine ComputeFKTables(inputfile,outputpath,Q0,flmap)
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
      character*50 outputpath
**
*     Internal Variables
*
      integer idat
      integer ounit
      integer ls,lp
      integer ipdf,jpdf
      double precision t1,t2
**
*     Output Variables
*
      integer flmap(0:13,0:13)
*
*     Initialize APFEL
*
c      call EnableWelcomeMessage(.false.)
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
      elseif(set.eq."DYE906R_BIN01")then
         set = "DYE906R_P_BIN01"
      elseif(set.eq."DYE906R_BIN02")then
         set = "DYE906R_P_BIN02"
      elseif(set.eq."DYE906R_BIN03")then
         set = "DYE906R_P_BIN03"
      elseif(set.eq."DYE906R_BIN04")then
         set = "DYE906R_P_BIN04"
      elseif(set.eq."DYE906R_BIN05")then
         set = "DYE906R_P_BIN05"
      elseif(set.eq."DYE906R_BIN06")then
         set = "DYE906R_P_BIN06"
      elseif(set.eq."DYE906R_BIN07")then
         set = "DYE906R_P_BIN07"
      elseif(set.eq."DYE906R_BIN08")then
         set = "DYE906R_P_BIN08"
      elseif(set.eq."DYE906R_BIN09")then
         set = "DYE906R_P_BIN09"
      elseif(set.eq."DYE906R_BIN10")then
         set = "DYE906R_P_BIN10"
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
*     Find the length of the outputpath.
*     (char(0) is the null char (for c++) while char(32) is space (for fortran))
      lp = index(outputpath,char(0)) - 1
      if(lp.eq.-1) lp = index(outputpath,char(32)) - 1
      open(unit=ounit,file=outputpath(1:lp)//"/FK_"//set(1:ls)//".dat",
     1     status="unknown")
      write(6,*) "Processing ",set(1:ls)," set ..."
      do idat=1,ndata
         write(6,"(a,i4,a,i4)") " Convoluting data number =",idat,
     1                          " /",ndata
         if(set.eq."DYE886R_P")then
            obs(idat) = "DYP_E886P"
         elseif(set.eq."DYE886R_D")then
            obs(idat) = "DYP_E886D"
         elseif(set.eq."DYE906R_P_BIN01"
     1           .or.set.eq."DYE906R_P_BIN02"
     1           .or.set.eq."DYE906R_P_BIN03"
     1           .or.set.eq."DYE906R_P_BIN03"
     1           .or.set.eq."DYE906R_P_BIN04"
     1           .or.set.eq."DYE906R_P_BIN05"
     1           .or.set.eq."DYE906R_P_BIN06"
     1           .or.set.eq."DYE906R_P_BIN07"
     1           .or.set.eq."DYE906R_P_BIN08"
     1           .or.set.eq."DYE906R_P_BIN09"
     1           .or.set.eq."DYE906R_P_BIN10")then
            obs(idat) = "DYP_E906P"
         elseif(set.eq."DYE906R_D_BIN01"
     1           .or.set.eq."DYE906R_D_BIN02"
     1           .or.set.eq."DYE906R_D_BIN03"
     1           .or.set.eq."DYE906R_D_BIN04"
     1           .or.set.eq."DYE906R_D_BIN05"
     1           .or.set.eq."DYE906R_D_BIN06"
     1           .or.set.eq."DYE906R_D_BIN07"
     1           .or.set.eq."DYE906R_D_BIN08"
     1           .or.set.eq."DYE906R_D_BIN09"
     1           .or.set.eq."DYE906R_D_BIN10")then
            obs(idat) = "DYP_E906D"   
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
      elseif(set.eq."DYE906R_P_BIN01")then
         set = "DYE906R_D_BIN01"
         goto 103
      elseif(set.eq."DYE906R_P_BIN02")then
         set = "DYE906R_D_BIN02"
         goto 103
      elseif(set.eq."DYE906R_P_BIN03")then
         set = "DYE906R_D_BIN03"
         goto 103
      elseif(set.eq."DYE906R_P_BIN04")then
         set = "DYE906R_D_BIN04"
         goto 103
      elseif(set.eq."DYE906R_P_BIN05")then
         set = "DYE906R_D_BIN05"
         goto 103
      elseif(set.eq."DYE906R_P_BIN06")then
         set = "DYE906R_D_BIN06"
         goto 103
      elseif(set.eq."DYE906R_P_BIN07")then
         set = "DYE906R_D_BIN07"
         goto 103
      elseif(set.eq."DYE906R_P_BIN08")then
         set = "DYE906R_D_BIN08"
         goto 103
      elseif(set.eq."DYE906R_P_BIN09")then
         set = "DYE906R_D_BIN09"
         goto 103
      elseif(set.eq."DYE906R_P_BIN10")then
         set = "DYE906R_D_BIN10"
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
