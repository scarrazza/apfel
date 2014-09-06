************************************************************************
*
*     DISObservables.f:
*
*     Test program for the DIS observables.
*
************************************************************************
      program DISObservables
*
      implicit none
**
*     Variables
*
      integer pto
      integer ip,np,irep
      double precision x,q2,q,q20,q0,y,pol
      double precision xin,xfi
      double precision F2(3:7),F3(3:7),FL(3:7),sigma(3:7)
      double precision t2,t1
      character*2  proc
      character*5  scheme
      character*53 pdfset
      character*9  target
      character*12 proj
*
      call cpu_time(t1)
*
*     Define kinematics
*
      q20 = 1d0
      q2  = 10d0
      y   = 0.5d0
      pol = 0d0
      xin = 1d-7
      xfi = 0.9999d0
      np  = 50
      q0  = dsqrt(q20)
      q   = dsqrt(q2)
*
*     Define parameters of the computation
*
      proj   = "ELECTRON"
      target = "PROTON"
      pto    = 2
      scheme = "FONLL"
      irep   = 0
      if(pto.eq.0)then
         pdfset = "toyLH_LO.LHgrid"
      elseif(pto.eq.1)then
c         pdfset = "toyLH_NLO.LHgrid"
         pdfset = "NNPDF30_nlo_as_0118.LHgrid"
      elseif(pto.eq.2)then
c         pdfset = "toyLH_NNLO.LHgrid"
         pdfset = "NNPDF30_nnlo_as_0118.LHgrid"
      endif
c      pdfset = "APFEL"
c      call SetPoleMasses(dsqrt(2d0),1d5,1d5)
*
      do irep=0,100
*     
*     Electro-magnetic Observables
*     
         proc   = "NC"
*     
         write(6,*) irep
c         write(6,*) "   x    ",
c     1        "   Q2   ",
c     2        "    F2    ",
c     3        "    FL    ",
c     4        "    F3    "
c         write(6,*) "   "
*     
         x = xin
         do ip=1,np
            call DIS_xsec(x,q0,q,y,pol,proc,scheme,pto,pdfset,irep,
     1           target,proj,F2,F3,FL,sigma)
*     
c            write(*,44) x,q2,F2(7),FL(7),F3(7)
            write(*,*) x,q2,F2(7),FL(7),F3(7)
            x = x * dexp( 1d0 / ( np - 1d0 ) * dlog( xfi / xin ) )
         enddo
         write(6,*) "   "
         write(6,*) "   "
c$$$*     
c$$$*     Charge-current Obserables
c$$$*     
c$$$         proc   = "CC"
c$$$*     
c$$$         write(6,*) "   x    ",
c$$$     1        "   Q2   ",
c$$$     2        "    F2    ",
c$$$     3        "    FL    ",
c$$$     4        "    F3    "
c$$$         write(6,*) "   "
c$$$*     
c$$$         x = xin
c$$$         do ip=1,np
c$$$            call DIS_xsec(x,q0,q,y,pol,proc,scheme,pto,pdfset,irep,
c$$$     1           target,proj,F2,F3,FL,sigma)
c$$$*     
c$$$            write(*,44) x,q2,F2(7),FL(7),F3(7)
c$$$            x = x * dexp( 1d0 / ( np - 1d0 ) * dlog( xfi / xin ) )
c$$$         enddo
c$$$         write(6,*) "   "
      enddo
*
      call cpu_time(t2)
      write(6,*) "time elapsed:",t2-t1," s"
      write(6,*) "   "
*
 44   format(f9.7,1x,f8.4,1x,3(f11.8,1x))
*
      end
