************************************************************************
*
*     gauleg.f:
*
*     Routine fo the Gauss integration.
*
************************************************************************
      subroutine gauleg(x1,x2,x,w,n)
*
      implicit none
*
      integer n
      double precision x1,x2,x(n),w(n)
      double precision EPS
      parameter (EPS=3.d-14)
      integer i,j,m
      double precision p1,p2,p3,pp,xl,xm,z,z1
*
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
c      print*,'m xm xl ',m,xm,xl
      do i=1,m
         z=cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
 1       continue
         p1=1.d0
         p2=0.d0
         do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         end do
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp
         if (abs(z-z1).gt.EPS) goto 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
      end do
      return
      end
