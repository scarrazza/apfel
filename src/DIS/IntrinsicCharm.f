************************************************************************
*
*     O(as) coefficient functions for massive initial state charm.
*     The following expressions hold both for NC and CC processes.
*     The z-independent kinematics is bassed by common.
*
************************************************************************
*     Regular functions that multiplies 1 / ( 1 - z )
************************************************************************
      function c11ICR(Q2,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision s1h
      double precision Lxi,Ixi
      double precision DeltaFun
      double precision f1hat,N1
**
*     Output Variables
*
      double precision c11ICR
*
      Del2  = Del * Del
      s1h   = ( 1d0 - z ) * ( ( Del - Spm ) * z + Del + Spm ) / 2d0 / z
      s1h2  = s1h * s1h
      Delp  = DeltaFun(m12,s1h+m12,-Q2)
      Delp2 = Delp * Delp
      Lxi   = dlog( ( Spp + s1h - Delp ) / ( Spp + s1h + Delp ) )
      Ixi   = ( ( s1h + 2d0 * m22 ) / s1h2 + ( s1h + m22 )
     1      / Delp / s1h2 * Spp * Lxi )
*
      f1hat = 8d0 / Delp2 * ( - Del2 * ( Splus * Spp
     1     - 2d0 * m1 * m2 * Sminus ) * Ixi(z,Q)
     2     + 2d0 * m1 * m2 * Sminus * ( 1d0 / s1h * ( Delp2
     3     + 4d0 * m22 * Spm )
     4     + 2d0 * Spm - Smp + (Spp + s1h) / 2d0
     5     + ( s1h + m22 ) / Delp / s1h * ( Delp2
     6     + 2d0 * Spm * Spp + ( m22 + Q2 ) * s1h ) * Lxi )
     7     + Splus * ( ( - m22 * Spp ) / ( ( s1h + m22 ) * s1h )
     8     * ( Del2 + 4d0 * m22 * Spm)
     9     - 1d0 / 4d0 / ( s1h + m22 )
     1     * ( 3d0 * Spp**2d0 * Smp
     2     + 4d0 * m22 * (10d0 * Spp * Spm - Spm * Smp
     3     - m12 * Spp)
     4     + s1h * ( - 7d0 * Spp * Smp + 18d0 * Del2
     5     - 4d0 * m12 * ( 7d0 * Q2 - 4 * m22
     6     + 7d0 * m12 ) )
     7     + 3d0 * s1h2 * ( Spm - 2d0 * m12 ) - s1h**3d0 )
     8     + ( s1h + m22 ) / 2d0 / Delp
     9     * ( - 2d0 / s1h * Spp * ( Del2 + 2d0 * Spm * Spp )
     1     + ( 4d0 * m12 * m22 - 7d0 * Spm * Spp )
     2     - 4d0 * Spm * s1h - s1h*2d0c ) * Lxi ) )
*
      N1 = ( Splus * Spp - 2d0 * m1 * m2 * Sminus ) / 2d0 /  Del
*
      c11ICR = fact1 * 2d0 * CF * ( 1 - z ) * s1h * f1hat / N1 / 8d0 
     1       / ( s1h + m22 )
*
      return
      end
*
************************************************************************
      function c21ICR(Q2,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision s1h
      double precision Lxi,Ixi
      double precision DeltaFun
      double precision f2hat,N2
**
*     Output Variables
*
      double precision c21ICR
*
      Del2  = Del * Del
      s1h   = ( 1d0 - z ) * ( ( Del - Spm ) * z + Del + Spm ) / 2d0 / z
      s1h2  = sih * s1h
      Delp  = DeltaFun(m12,s1h+m12,-Q2)
      Delp2 = Delp * Delp
      Lxi   = dlog( ( Spp + s1h - Delp ) / ( Spp + s1h + Delp ) )
      Ixi   = ( ( s1h + 2d0 * m22 ) / s1h2 + ( s1h + m22 )
     1      / Delp / s1h2 * Spp * Lxi )
*
      f2hat = 16d0 / Delp**4d0 * ( - 2d0 * Del2**4d0 * Splus * Ixi
     1     + 2d0 * m1 * m2 * Sminus * ( ( ( s1h + m22 ) / Delp )
     2     * ( Delp2 - 6d0 * m12 * Q2 ) * Lxi
     3     - Delp2 * ( s1h + Spp ) / 2d0 / ( s1h + m22 )
     4     + ( 2d0 * Delp2 - 3d0 * Q2 * ( s1h + Spp) ) )
     5     + Splus * ( - 2d0 * ( Del2 - 6d0 * m12 * Q2 )
     6     * ( s1h + m22 ) - 2d0 * ( m12 + m22 ) * s1h2
     7     - 9d0 * m22 * Spm**2d0 + Del2 * ( 2d0 * Spp - m22 )
     8     + 2d0 * s1h * ( 2d0 * Del2 + ( m12 - 5d0 * m22 ) * Spm )
     9     + ( Delp2 - 6d0 * Q2 * ( m22 + s1h ) )
     1     * Spp * ( s1h + Spp ) / 2d0 / ( s1h + m22 )
     2     - 2d0 * Del2 / s1h * ( Del2
     3     + 2d0 * ( 2d0 * m22 + s1h ) * Spm )
     4     + ( s1h + m22 ) / Delp * ( - 2d0 / s1h * Del2
     5     * ( Del2 + 2d0 * Spm * Spp )
     6     - 2d0 * s1h * ( Del2 - 6d0 * m12 * Q2 )
     7     - ( Delp2 - 18d0 * m12 * Q2 ) * Spp
     8     - 2d0 * Del2 * ( Spp + 2d0 * Spm) ) * Lxi ) )
*
      N2 = 2d0 * Splus * Del / Delp2
*
      c21ICR = fact2 * 2d0 * CF * ( 1 - z ) * s1h * f2hat / N2 / 8d0 
     1       / ( s1h + m22 ) 
      return
      end
*
************************************************************************
      function c31ICR(Q2,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision s1h
      double precision Lxi,Ixi
      double precision DeltaFun
      double precision f3hat,N3
**
*     Output Variables
*
      double precision c31ICR
*
      Del2  = Del * Del
      s1h   = ( 1d0 - z ) * ( ( Del - Spm ) * z + Del + Spm ) / 2d0 / z
      s1h2  = sih * s1h
      Delp  = DeltaFun(m12,s1h+m12,-Q2)
      Delp2 = Delp * Delp
      Lxi   = dlog( ( Spp + s1h - Delp ) / ( Spp + s1h + Delp ) )
      Ixi   = ( ( s1h + 2d0 * m22 ) / s1h2 + ( s1h + m22 )
     1      / Delp / s1h2 * Spp * Lxi )
*
      f3hat = 16d0 / Delp2 * ( - 2d0 * Del2 * Rplus * Ixi
     1     + 2d0 * m1 * m2 * Rminus * ( 1d0 - Smp / s1h + ( s1h + m22 )
     2     * ( s1h + Spm ) / Delp / s1h * Lxi )
     3     + Rplus * ( Smp - 3d0 * Spm - 2d0 / s1h * ( Del2
     4     + 2d0 * m22 * Spm ) - ( s1h - Smp ) * ( s1h + Spp ) / 2d0
     5     / ( s1h + m22 ) + ( s1h + m22 ) / Delp / s1h * ( - s1h2
     6     + 4d0 * ( m12 * Smp - Del2 ) - 3d0 * s1h * Spm ) * Lxi )
*
      N3 = 2d0 * Rplus / Delp
*
      c31ICR = fact3 * 2d0 * CF * ( 1 - z ) * s1h * f3hat / N3 / 8d0 
     1       / ( s1h + m22 ) 
      return
      end
*
************************************************************************
      function cL1ICR(Q2,z)
*
      implicit none
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision c11ICR,c21ICR
**
*     Output Variables
*
      double precision cL1ICR
*
      cL1ICR = factL * ( c21ICR(Q2,z) - c11ICR(Q2,z) )
*
      return
      end
*
************************************************************************
*     Local terms proportional to delta(1 - z)
*     + the residual term coming frm the plus prescripted term.
************************************************************************
      function c11ICL(Q2,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision c11ICR
      double precision one
      parameter(one=0.99999999d0)
**
*     Output Variables
*
      double precision c11ICL
*
      c11ICL = fact1 * ( S1 + V1 ) + c11ICR(Q2,one) * dlog(1d0-z)
*
      return
      end
*
************************************************************************
      function c21ICL(Q2,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision Q2,z
      double precision one
      parameter(one=0.99999999d0)
**
*     Internal Variables
*
      double precision c21ICR
**
*     Output Variables
*
      double precision c21ICL
*
      c21ICL = fact2 * ( S2 + V2 ) + c21ICR(Q2,one) * dlog(1d0-z)
*
      return
      end
*
************************************************************************
      function c31ICL(Q2,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrapIC.h"
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision c31ICR
      double precision one
      parameter(one=0.99999999d0)
**
*     Output Variables
*
      double precision c31ICL
*
      c31ICL = fact3 * ( S3 + V3 ) + c31ICR(Q2,one) * dlog(1d0-z)
*
      return
      end
*
************************************************************************
      function cL1ICL(Q2,z)
*
      implicit none
**
*     Input Variables
*
      double precision Q2,z
**
*     Internal Variables
*
      double precision c11ICL,c21ICL
**
*     Output Variables
*
      double precision cL1ICL
*
      cL1ICL = factL * ( c21ICL(Q2,z) - c11ICL(Q2,z) )
*
      return
      end
*
************************************************************************
      function DeltaFun(a,b,c)
*
      implicit none
**
*     Input Variables
*
      double precision a,b,c
**
*     Output Variables
*
      double precision DeltaFun
*
      DeltaFun = dsqrt( a * a + b * b + c * c
     1         - 2d0 * ( a * b + b * c + c * a ) )
*
      return
      end
