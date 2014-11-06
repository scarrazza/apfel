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
      q20 = 2d0
      q2  = 10d0
      y   = 0.5d0
      pol = 0d0
      xin = 1d-5
      xfi = 1d-1
      np  = 5
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
         pdfset = "toyLH_NLO.LHgrid"
      elseif(pto.eq.2)then
         pdfset = "toyLH_NNLO.LHgrid"
      endif
c      pdfset = "APFEL"
c      call SetPoleMasses(dsqrt(2d0),1d5,1d5)
*
*     Electro-magnetic Observables
*
      proc   = "EM"
*
      write(6,*) "   x    ",
     1           "   Q2   ",
     2           "   y    ",
     3           "    F2l   ",
     4           "    F2c   "
      write(6,*) "   "
*
      x = xin
      do ip=1,np
         call DIS_xsec(x,q0,q,y,pol,proc,scheme,pto,pdfset,irep,target,
     1                 proj,F2,F3,FL,sigma)
*
         write(*,44) x,q2,y,F2(3),F2(4)
         x = x * dexp( 1d0 / ( np - 1d0 ) * dlog( xfi / xin ) )
      enddo
      write(6,*) "   "
*
*     Charge-current Obserables
*
      proc   = "CC"
*
      write(6,*) "   x    ",
     1           "   Q2    ",
     2           "   y    ",
     3           "      F2c     ",
     4           "      FLc     ",
     5           "      F3c     ",     
     6           "    sigmac    "
      write(6,*) "   "
*
      x = xin
      do ip=1,np
         call DIS_xsec(x,q0,q,y,pol,proc,scheme,pto,pdfset,irep,target,
     1                 proj,F2,F3,FL,sigma)
*
         write(*,45) x,q2,y,F2(4),FL(4),F3(4),sigma(4)
         x = x * dexp( 1d0 / ( np - 1d0 ) * dlog( xfi / xin ) )
      enddo
      write(6,*) "   "
*
      call cpu_time(t2)
      write(6,*) "time elapsed:",t2-t1," s"
      write(6,*) "   "
*
 44   format(f7.5,1x,f8.4,1x,f7.5,1x,7(f9.6,1x))
 45   format(f7.5,1x,f8.4,1x,f7.5,1x,4(f13.8,1x))
*
      end
