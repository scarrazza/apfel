************************************************************************
*
*     toyLHPDFs.f:
*
*     This routine returns the toyLH PDFs at the intitial scale
*     which is supposed to be Q = sqrt(2) GeV.
*
*     Appended at the end of this file we also provide a "private" set
*     of PDF that can be used for private tests.
*
************************************************************************
      subroutine toyLHPDFs(x,xpdf)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer ipdf
      double precision xuv,xdv,xg,xdbar,xubar,xs,xsbar
      double precision N_uv,auv,buv,N_dv,adv,bdv,N_g,ag
      double precision bg,N_db,adb,bdb,fs
**
*     Output Variables
*
      double precision xpdf(-6:6)
*
*     Parameters of the User defined PDFs
*
      N_uv = 5.107200d0
      auv  = 0.8d0
      buv  = 3d0
      N_dv = 3.064320d0
      adv  = 0.8d0
      bdv  = 4d0
      N_g  = 1.7d0
      ag   = -0.1d0
      bg   = 5d0
      N_db = 0.1939875d0
      adb  = -0.1d0
      bdb  = 6d0
      fs   = 0.2d0
*
*     User defined PDFs
*
      xuv   = N_uv * x**auv * ( 1d0 - x )**buv
      xdv   = N_dv * x**adv * ( 1d0 - x )**bdv
      xg    = N_g  * x**ag  * ( 1d0 - x )**bg 
      xdbar = N_db * x**adb * ( 1d0 - x )**bdb
      xubar = xdbar * ( 1d0 - x )
      xs    = fs * ( xdbar + xubar )
      xsbar = xs
*
*     Initialize PDFs to zero
*
      do ipdf=-6,6
         xpdf(ipdf) = 0d0
      enddo
*
      if(x.gt.1d0) return
*
      xpdf(3)  = xs
      xpdf(2)  = xuv + xubar
      xpdf(1)  = xdv + xdbar
      xpdf(0)  = xg
      xpdf(-1) = xdbar
      xpdf(-2) = xubar
      xpdf(-3) = xsbar
*
      return
      end
*
************************************************************************
      subroutine private(x,xpdf)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer ipdf
      double precision xuv,xdv,xg,xSS,xdbar,xubar,xs,xsbar
      double precision N_uv,auv,buv,N_dv,adv,bdv,N_g,ag,bg
      double precision N_S,aS,bS
**
*     Output Variables
*
      double precision xpdf(-6:6)
*
*     Parameters of the User defined PDFs
*
      N_uv = 35d0/16d0
      auv  = 0.5d0
      buv  = 3d0
      N_dv = 315d0/256d0
      adv  = 0.5d0
      bdv  = 4d0
      N_g  = 1.90836d0
      ag   = -0.2d0
      bg   = 5d0
      N_S  = 0.673345d0
      aS   = -0.2d0
      bS   = 7d0
c      N_uv = 2d0 * 3.69141d0 / 3d0
c      auv  = 0.5d0
c      buv  = 4d0
c      N_dv = 3.69141d0 / 3d0
c      adv  = auv
c      bdv  = buv
c      N_g  = 2.04996d0
c      ag   = - 0.18d0
c      bg   = 5d0
c      N_S  = N_g / 3d0
c      aS   = ag
c      bS   = bg
*
*     User defined PDFs
*
      xuv   = N_uv * x**auv * ( 1d0 - x )**buv
      xdv   = N_dv * x**adv * ( 1d0 - x )**bdv
      xg    = N_g  * x**ag  * ( 1d0 - x )**bg 
      xSS   = N_S  * x**aS  * ( 1d0 - x )**bS
      xubar = xSS / 6d0
      xdbar = xubar
      xsbar = xubar
      xs    = xsbar
*
*     Initialize PDFs to zero
*
      do ipdf=-6,6
         xpdf(ipdf) = 0d0
      enddo
*
      if(x.gt.1d0) return
*
      xpdf(3)  = xs
      xpdf(2)  = xuv + xubar
      xpdf(1)  = xdv + xdbar
      xpdf(0)  = xg
      xpdf(-1) = xdbar
      xpdf(-2) = xubar
      xpdf(-3) = xsbar
*
      return
      end
*
************************************************************************
*
*     Kretzer's parametrization at Q2 = 0.4 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0003177.
*
************************************************************************
      subroutine KretzerFFs(x,xff)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer iff
      double precision ag,bg,N_g
      double precision as,bs,N_s
      double precision al,bl,N_l
      double precision beta
**
*     Output Variables
*
      double precision xff(-6:6)
*
*     Security cutoff
*
      if(x.gt.1d0) x = 1d0
*
*     Parameters of the User defined PDFs
*
      al  = - 0.829d0
      bl  = 0.949d0
      N_l = 0.264d0 / beta(al+2d0,bl+1d0)
      as  = al
      bs  = bl + 1d0
      N_s = 0.165d0 / beta(as+2d0,bs+1d0)
      ag  = 4.374d0
      bg  = 9.778d0
      N_g = 0.215d0 / beta(ag+2d0,bg+1d0)
*
*     Initialize PDFs to zero
*
      do iff=-6,6
         xff(iff) = 0d0
      enddo
*
      xff(3)  = x * N_s * x**as * ( 1d0 - x )**bs
      xff(2)  = x * N_l * x**al * ( 1d0 - x )**bl
      xff(1)  = xff(3)
      xff(0)  = x * N_g * x**ag * ( 1d0 - x )**bg
      xff(-1) = xff(2)
      xff(-2) = xff(1)
      xff(-3) = xff(3)
*
      return
      end
*
************************************************************************
*
*     HKNS parametrization at Q2 = 1 GeV^2 of the light partons
*     for pi+ taken at NLO from hep-ph/0702250.
*
************************************************************************
      subroutine HKNSFFs(x,xff)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer iff
      double precision ag,bg,N_g
      double precision as,bs,N_s
      double precision al,bl,N_l
      double precision beta
**
*     Output Variables
*
      double precision xff(-6:6)
*
*     Security cutoff
*
      if(x.gt.1d0) x = 1d0
*
*     Parameters of the User defined PDFs
*
      al  = - 0.963d0
      bl  = 1.370d0
      N_l = 0.401d0 / beta(al+2d0,bl+1d0)
      as  = 0.718d0
      bs  = 6.266d0
      N_s = 0.094d0 / beta(as+2d0,bs+1d0)
      ag  = 1.943d0
      bg  = 8.000d0
      N_g = 0.238d0 / beta(ag+2d0,bg+1d0)
*
*     Initialize PDFs to zero
*
      do iff=-6,6
         xff(iff) = 0d0
      enddo
*
      xff(3)  = x * N_s * x**as * ( 1d0 - x )**bs
      xff(2)  = x * N_l * x**al * ( 1d0 - x )**bl
      xff(1)  = xff(3)
      xff(0)  = x * N_g * x**ag * ( 1d0 - x )**bg
      xff(-1) = xff(2)
      xff(-2) = xff(1)
      xff(-3) = xff(3)
*
      return
      end
