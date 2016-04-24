************************************************************************
*
*     CheckAPFEL.f:
*
*     Code that tests APFEL against a pretabulated set of results for:
*     - alpha_s,
*     - evolution,
*     - structure functions.
*
************************************************************************
      function CheckAPFEL()
*
      implicit none
*
      include "../commons/CheckAPFEL.h"
**
*     Internal Variables
*
      integer ilha,iQ,iref
      double precision Q0,Q
      double precision eps
      double precision toll
      double precision AlphaQCD
      double precision xPDFj
      double precision F2light,F2charm,F2bottom,F2total
      double precision FLlight,FLcharm,FLbottom,FLtotal
      double precision F3light,F3charm,F3bottom,F3total
      double precision Act(1048)
      double precision dev
      character*6 apfelversion
      logical succ
      parameter(eps=1d-10)
      parameter(toll=1d-6)
**
*     Output Variables
*
      logical CheckAPFEL
*
      call getapfelversion(apfelversion)
      write(6,*) achar(27)//"[34m"
      write(6,*) "Checking APFEL v",apfelversion," ...",achar(27)//"[0m"
*
*     Settings
*
c      call EnableWelcomeMessage(.false.)
      call SetMassScheme("FONLL-C")
      Q0  = dsqrt(Q20) - eps
*
*     Compute predictions
*
      call InitializeAPFEL_DIS
*
*     Compare the output with the reference
*
      iref = 0
      do iQ=1,4
         Q = dsqrt(Q2(iQ))
*     alpha_s
         iref = iref + 1
         Act(iref) = AlphaQCD(Q)
*     PDF evolution
         call EvolveAPFEL(Q0,Q)
         do ilha=1,9
            iref = iref + 1
            Act(iref) = xPDFj(2,xlha(ilha))-xPDFj(-2,xlha(ilha))
            iref = iref + 1
            Act(iref) = xPDFj(1,xlha(ilha))-xPDFj(-1,xlha(ilha))
            iref = iref + 1
            Act(iref) = 2d0*(xPDFj(-1,xlha(ilha))+xPDFj(-2,xlha(ilha)))
            iref = iref + 1
            Act(iref) = xPDFj(4,xlha(ilha))+xPDFj(-4,xlha(ilha)) 
            iref = iref + 1
            Act(iref) = xPDFj(0,xlha(ilha))
         enddo
*      NC Structure Functions
         call SetProcessDIS("NC")
         call ComputeStructureFunctionsAPFEL(Q0,Q)
         do ilha=1,9
*     F_2
            iref = iref + 1
            Act(iref) = F2light(xlha(ilha))
            iref = iref + 1
            Act(iref) = F2charm(xlha(ilha))
            iref = iref + 1
            Act(iref) = F2bottom(xlha(ilha))
            iref = iref + 1
            Act(iref) = F2total(xlha(ilha))
*     F_L
            iref = iref + 1
            Act(iref) = FLlight(xlha(ilha))
            iref = iref + 1
            Act(iref) = FLcharm(xlha(ilha))
            iref = iref + 1
            Act(iref) = FLbottom(xlha(ilha))
            iref = iref + 1
            Act(iref) = FLtotal(xlha(ilha))
*     F_3
            iref = iref + 1
            Act(iref) = F3light(xlha(ilha))
            iref = iref + 1
            Act(iref) = F3charm(xlha(ilha))
            iref = iref + 1
            Act(iref) = F3bottom(xlha(ilha))
            iref = iref + 1
            Act(iref) = F3total(xlha(ilha))
         enddo
*      CC Structure Functions
         call SetProcessDIS("CC")
         call ComputeStructureFunctionsAPFEL(Q0,Q)
         do ilha=1,9
*     F_2
            iref = iref + 1
            Act(iref) = F2light(xlha(ilha))
            iref = iref + 1
            Act(iref) = F2charm(xlha(ilha))
            iref = iref + 1
            Act(iref) = F2bottom(xlha(ilha))
            iref = iref + 1
            Act(iref) = F2total(xlha(ilha))
*     F_L
            iref = iref + 1
            Act(iref) = FLlight(xlha(ilha))
            iref = iref + 1
            Act(iref) = FLcharm(xlha(ilha))
            iref = iref + 1
            Act(iref) = FLbottom(xlha(ilha))
            iref = iref + 1
            Act(iref) = FLtotal(xlha(ilha))
*     F_3
            iref = iref + 1
            Act(iref) = F3light(xlha(ilha))
            iref = iref + 1
            Act(iref) = F3charm(xlha(ilha))
            iref = iref + 1
            Act(iref) = F3bottom(xlha(ilha))
            iref = iref + 1
            Act(iref) = F3total(xlha(ilha))
         enddo
      enddo
*
*     Compare present predictions with the reference
*
      succ = .true.
      do iref=1,1048
         dev = dabs((Ref(iref)-Act(iref))/Ref(iref))
         if(dev.gt.toll) succ = .false.
c         if(.not.succ) write(6,*) iref,dev,Ref(iref),Act(iref)
c         succ = .true.
      enddo
*
      if(succ)then
         write(6,*) "Check ... ",
     1        achar(27)//"[1;32m","succeded",achar(27)//"[0m"
      else
         write(6,*) "Check ... ",
     1        achar(27)//"[1;31m","failed",achar(27)//"[0m"
      endif
      write(6,*)
      CheckAPFEL = succ
*
*     Code to write the reference table
*
c$$$      open(unit=9,file="CheckAPFEL.h")
c$$$      write(9,'(a)') "*     -*-fortran-*-"
c$$$      write(9,'(a)') "*"
c$$$      write(9,'(a)') "*     Tabutalet values to check APFEL"
c$$$      write(9,'(a)') "*"
c$$$      write(9,'(a)') "      double precision Q20,Q2(4)"
c$$$      write(9,'(a)') "      double precision xlha(9)"
c$$$      write(9,'(a)') "      double precision Ref(1049)"
c$$$      write(9,'(a,es10.3,a)') "      data Q20  /",Q20,"/"
c$$$      write(9,'(a,4(es10.3,a))') "      data Q2   /",Q2(1),",",
c$$$     1                                               Q2(2),",",
c$$$     2                                               Q2(3),",",
c$$$     3                                               Q2(4),"/"
c$$$      write(9,'(a,5(es10.3,a))') "      data xlha /",xlha(1),",",
c$$$     1                                               xlha(2),",",
c$$$     2                                               xlha(3),",",
c$$$     3                                               xlha(4),",",
c$$$     4                                               xlha(5),","
c$$$      write(9,'(a,4(es10.3,a))') "     1           ",xlha(6),",",
c$$$     1                                               xlha(7),",",
c$$$     2                                               xlha(8),",",
c$$$     3                                               xlha(9),"/"
c$$$      write(9,'(a)') "      data Ref /"
c$$$      do iQ=1,4
c$$$         Q = dsqrt(Q2(iQ))
c$$$*     alpha_s
c$$$         write(9,'(a,es14.7,a)') "     1   ",AlphaQCD(Q),","
c$$$*     PDF evolution
c$$$         call EvolveAPFEL(Q0,Q)
c$$$         do ilha=1,9
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           xPDFj(2,xlha(ilha)) - xPDFj(-2,xlha(ilha)),",",
c$$$     2           xPDFj(1,xlha(ilha)) - xPDFj(-1,xlha(ilha)),",",
c$$$     3           2d0*(xPDFj(-1,xlha(ilha)) + xPDFj(-2,xlha(ilha))),",",
c$$$     4           xPDFj(4,xlha(ilha)) + xPDFj(-4,xlha(ilha)),","
c$$$            write(9,'(a,es14.7,a)') "     1   ",
c$$$     1           xPDFj(0,xlha(ilha)),","
c$$$         enddo
c$$$*      NC Structure Functions
c$$$         call SetProcessDIS("NC")
c$$$         call ComputeStructureFunctionsAPFEL(Q0,Q)
c$$$         do ilha=1,9
c$$$*     F_2
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           F2light(xlha(ilha)),",",
c$$$     2           F2charm(xlha(ilha)),",",
c$$$     3           F2bottom(xlha(ilha)),",",
c$$$     4           F2total(xlha(ilha)),","
c$$$*     F_L
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           FLlight(xlha(ilha)),",",
c$$$     2           FLcharm(xlha(ilha)),",",
c$$$     3           FLbottom(xlha(ilha)),",",
c$$$     4           FLtotal(xlha(ilha)),","
c$$$*     F_3
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           F3light(xlha(ilha)),",",
c$$$     2           F3charm(xlha(ilha)),",",
c$$$     3           F3bottom(xlha(ilha)),",",
c$$$     4           F3total(xlha(ilha)),","
c$$$         enddo
c$$$*      CC Structure Functions
c$$$         call SetProcessDIS("CC")
c$$$         call ComputeStructureFunctionsAPFEL(Q0,Q)
c$$$         do ilha=1,9
c$$$*     F_2
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           F2light(xlha(ilha)),",",
c$$$     2           F2charm(xlha(ilha)),",",
c$$$     3           F2bottom(xlha(ilha)),",",
c$$$     4           F2total(xlha(ilha)),","
c$$$*     F_L
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           FLlight(xlha(ilha)),",",
c$$$     2           FLcharm(xlha(ilha)),",",
c$$$     3           FLbottom(xlha(ilha)),",",
c$$$     4           FLtotal(xlha(ilha)),","
c$$$*     F_3
c$$$            write(9,'(a,4(es14.7,a))') "     1   ",
c$$$     1           F3light(xlha(ilha)),",",
c$$$     2           F3charm(xlha(ilha)),",",
c$$$     3           F3bottom(xlha(ilha)),",",
c$$$     4           F3total(xlha(ilha)),","
c$$$         enddo
c$$$      enddo
c$$$      write(9,'(a)') "     1    0.0000000E+00/"
c$$$      close(9)
*
      return
      end

