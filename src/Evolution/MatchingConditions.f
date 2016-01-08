************************************************************************
*
*     MatchingConditions.f:
*
*     This is a collection of funtions that are used to match PDFs and
*     FFs.
*
************************************************************************
*
*     Space-like matching conditions (PDFs)
*     Reference: Appendix B of hep-ph/9612398
*
************************************************************************
*     Eq. (B.1)
************************************************************************
      function APS2Hq(z)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision lnz,z2,Li21mz,S121mz
      double precision wgplg,ddilog
      double precision A0
**
*     Output Variables
*
      double precision APS2Hq
*
*     some definitions
*
      lnz    = dlog(z)
      z2     = z * z
      Li21mz = ddilog(1d0-z)
      S121mz = wgplg(1,2,1d0-z)
*
      A0 = ( 1d0 + z ) * ( 32d0 * S121mz + 16d0 * lnz * Li21mz 
     1   - 16d0 * zeta2 * lnz - 4d0 * lnz**3 / 3d0 ) 
     2   + ( 32d0 / 3d0 / z + 8d0 - 8d0 * z - 32d0 * z2 / 3d0 ) 
     3   * ( Li21mz - zeta2 ) + ( 2d0 + 10d0 * z + 16d0 * z2 / 3d0 ) 
     4   * lnz**2 - ( 56d0 / 3d0 + 88d0 * z / 3d0 
     5   + 448d0 * z2 / 9d0 ) * lnz - 448d0 / 27d0 / z - 4d0 / 3d0
     6   - 124d0 * z / 3d0 + 1600d0 * z2 / 27d0
*
      APS2Hq = CF * TR * A0
*
      return
      end
*
************************************************************************
*     Eq. (B.2) ...without log (needed for MSbar masses)
************************************************************************
      function AS1Hg(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision B0
**
*     Output Variables
*
      double precision AS1Hg
*
      B0 = 4d0 * ( z**2 + ( 1d0 - z )**2 )
*
      AS1Hg = TR * B0
*
      return
      end
*
************************************************************************
*     Eq. (B.3)
************************************************************************
      function AS2Hg(z)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
      include "../commons/mass_scheme.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision z2
      double precision S121mz,S12mz,S211mz,S21mz,S111mz,S11mz
      double precision lnz,lnz2,lnz3
      double precision ln1mz,ln1mz2,ln1mz3
      double precision ln1pz,ln1pz2
      double precision wgplg,ddilog
      double precision A01,A02,B01,B02
      double precision h1,AS1Hg
**
*     Output Variables
*
      double precision AS2Hg
*
*     some definitions
*
      S121mz = wgplg(1,2,1d0-z) 
      S12mz  = wgplg(1,2,-z) 
      S211mz = wgplg(2,1,1d0-z) 
      S21mz  = wgplg(2,1,-z) 
      S111mz = ddilog(1d0-z)
      S11mz  = ddilog(-z)
*
      z2     = z * z
      lnz    = dlog(z)
      lnz2   = lnz * lnz 
      lnz3   = lnz2 * lnz
      ln1mz  = dlog(1d0 - z)
      ln1mz2 = ln1mz * ln1mz 
      ln1mz3 = ln1mz2 * ln1mz 
      ln1pz  = dlog(1d0 + z)
      ln1pz2 = ln1pz * ln1pz 
*     CF * TR  part
      A01 = ( 1d0 - 2d0 * z + 2d0 * z2 ) * ( 8d0 * zeta3 
     1    + 4d0 * ln1mz3 / 3d0 - 8d0 * ln1mz * s111mz 
     2    + 8d0 * zeta2 * lnz - 4d0 * lnz * ln1mz2 + 2d0 * lnz3 / 3d0 
     3    - 8d0 * lnz * s111mz + 8d0 * s211mz - 24d0 * s121mz )               
      B01 = - ( 4d0 + 96d0 * z - 64d0 * z2 ) * s111mz 
     1    - ( 4d0 - 48d0 * z + 40d0 * z2 ) * zeta2 
     2    - ( 8d0 + 48d0 * z - 24d0 * z2 ) * lnz * ln1mz 
     3    + ( 4d0 + 8d0 * z - 12d0 * z2 ) * ln1mz2 
     4    - ( 1d0 + 12d0 * z - 20d0 * z2 ) * lnz2 
     5    - ( 52d0 * z - 48d0 * z2 ) * ln1mz 
     6    - ( 16d0 + 18d0 * z + 48d0 * z2 ) * lnz 
     7    + 26d0 - 82d0 * z + 80d0 * z2 + z2 * ( - 16d0 * zeta2 * lnz 
     8    + 4d0 * lnz3 / 3d0 +  16d0 * lnz * s111mz +  32d0 * s121mz )
*     CA * TR  part
      A02 = ( 1d0 - 2d0 * z + 2d0 * z2 ) * ( - 4d0 * ln1mz3 / 3d0 
     1    + 8d0 * ln1mz * s111mz - 8d0 * s211mz ) 
     2    + ( 1d0 + 2d0 * z + 2d0 * z2 ) * ( - 8d0 * zeta2 * ln1pz 
     3    - 16d0 * ln1pz * s11mz - 8d0 * lnz * ln1pz2 
     4    + 4d0 * lnz2 * ln1pz + 8d0 * lnz * s11mz - 8d0 * s21mz 
     5    - 16d0 * s12mz ) + ( 16d0 + 64d0 * z ) * ( 2d0 * s121mz 
     6    + lnz * s111mz ) - ( 4d0 + 8d0 * z ) * lnz3 / 3d0 
     7    + ( 8d0 - 32d0 * z + 16d0 * z2 ) * zeta3 
     8    - ( 16d0 + 64d0 * z ) * zeta2 * lnz
      B02 = ( 16d0 * z + 16d0 * z2 ) * ( s11mz + lnz * ln1pz ) ! There is a typo here in the e-Print ( 16 -> 16z )
     1    + ( 32d0 / z / 3d0 + 12d0 + 64d0 * z - 272d0 * z2 / 3d0 ) 
     2    * s111mz - ( 12d0 + 48d0 * z - 260d0 * z2 / 3d0 
     3    + 32d0 / z / 3d0 ) * zeta2 - 4d0 * z2 * lnz * ln1mz 
     4    - ( 2d0 + 8d0 * z - 10d0 * z2 ) * ln1mz2 
     5    + ( 2d0 + 8d0 * z + 46d0 * z2 / 3d0 ) * lnz2 
     6    + ( 4d0 + 16d0 * z - 16d0 * z2 ) * ln1mz
     7    - ( 56d0 / 3d0 + 172d0 * z / 3d0 + 1600d0 * z2 / 9d0 ) * lnz 
     8    - 448d0 / z / 27d0 - 4d0 / 3d0 - 628d0 * z / 3d0 
     9    + 6352d0 * z2 / 27d0
*
      AS2Hg = TR * ( CF * ( A01 + B01 )  +  CA * ( A02 + B02 ) )
*
*     In case the the MSbar for the heavy querk masses is chosen
*
      if(mass_scheme.eq."MSbar")then
         h1 = 4d0 * CF
*
         AS2Hg = AS2Hg - 2d0 * h1 * AS1Hg(z)
      endif
*
      return
      end
*
************************************************************************
*     Eq. (B.4) Regular piece
************************************************************************
      function ANS2qqH_R(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision z2,lnz,lnz2,ln1mz
      double precision A0
**
*     Output Variables
*
      double precision ANS2qqH_R
*
*     some definitions
*
      z2    = z * z
      lnz   = dlog(z)
      lnz2  = lnz * lnz
      ln1mz = dlog(1d0 - z)
*
      A0 = ( 1d0 + z2 ) * ( 2d0 * lnz2 / 3d0 
     1   + 20d0 * lnz / 9d0 ) / ( 1d0 - z ) 
     2   + 8d0 * ( 1d0 - z ) * lnz / 3d0 
     3   + 44d0 / 27d0 - 268d0 * z / 27d0 
*
      ANS2qqH_R = CF * TR * A0
*
      return
      end
*
************************************************************************
*     Eq. (B.4) Singular piece
************************************************************************
      function ANS2qqH_S(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision A0
**
*     Output Variables
*
      double precision ANS2qqH_S
*
      A0 = 224d0 / 27d0
*
      ANS2qqH_S = CF * TR * A0 / ( 1d0 - z )
*
      return
      end
*
************************************************************************
*     Eq. (B.4) Local piece
************************************************************************
      function ANS2qqH_L(z)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision ln1mz
      double precision A0
**
*     Output Variables
*
      double precision ANS2qqH_L
*
*     some definitions
*
      ln1mz = dlog(1d0 - z)
*
      A0 = - 8d0 * zeta3 / 3d0 + 40d0 * zeta2 / 9d0 + 73d0 / 18d0
     1   + 224d0 * ln1mz / 27d0
*
      ANS2qqH_L = CF * TR * A0
*
      return
      end
*
************************************************************************
*     Eq. (B.5)
************************************************************************
      function AS2gqH(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision ln1mz
      double precision A0
**
*     Output Variables
*
      double precision AS2gqH
*
*     some definitions
*
      ln1mz = dlog(1d0 - z)
*
      A0 = 4d0 * ( 2d0 / z - 2d0 + z ) * ln1mz**2 / 3d0 
     1   + 8d0 * ( 10d0 / z - 10d0 + 8d0 * z ) * ln1mz / 9d0 
     2   + ( 448d0 / z - 448d0 + 344d0 * z ) / 27d0
*
      AS2gqH = CF * TR * A0
*
      return
      end
*
************************************************************************
*     Eq. (B.6) ...without log (needed for MSbar masses)
************************************************************************
      function AS1ggH_L()
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Internal Variables
*
      double precision B0
**
*     Output Variables
*
      double precision AS1ggH_L
*
      B0 = - 4d0 / 3d0
*
      AS1ggH_L = TR * B0
*
      return
      end
*
************************************************************************
*     Eq. (B.7) Regular piece
************************************************************************
      function AS2ggH_R(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision z2,lnz,lnz2,lnz3,ln1mz
      double precision A01,A02
**
*     Output Variables
*
      double precision AS2ggH_R
*
*     some definitions
*
      z2    = z * z
      lnz   = dlog(z)
      lnz2  = lnz * lnz; 
      lnz3  = lnz2 * lnz
      ln1mz = log(1d0 - z)
*     CF * TR part
      A01 = 4d0 * ( 1d0 + z ) * lnz3 / 3d0 + ( 6d0 + 10d0 * z) * lnz2 
     1    + ( 32d0 + 48d0 * z ) * lnz - 8d0 / z + 80d0 - 48d0 * z 
     2    - 24 * z2 
*     CA * TR part
      A02 = 4d0 * ( 1d0 + z ) * lnz2 / 3d0 
     1    + ( 52d0 + 88d0 * z ) * lnz / 9d0
     2    - 4d0 * z * ln1mz / 3d0 
     3    + ( 556d0 / z - 628d0 + 548d0 * z - 700d0 * z2 ) / 27d0 
*
      AS2ggH_R = TR * ( CF * A01 + CA * A02 )
*
      return
      end
*
************************************************************************
*     Eq. (B.7) Singular piece
************************************************************************
      function AS2ggH_S(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision A0
**
*     Output Variables
*
      double precision AS2ggH_S
*
      A0 = 224d0 / 27d0
*
      AS2ggH_S = CA * TR * A0 / ( 1d0 - z )
*
      return
      end
*
************************************************************************
*     Eq. (B.7) Local piece
************************************************************************
      function AS2ggH_L(z)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/ColorFactors.h"
      include "../commons/mass_scheme.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision ln1mz
      double precision A01,A02
      double precision h1,AS1ggH_L
**
*     Output Variables
*
      double precision AS2ggH_L
*
*     some definitions
*
      ln1mz = dlog(1d0 - z)
*     CF * TR part
      A01 = - 15d0
*     CA * TR part
      A02 = 10d0 / 9d0 + 224d0 * ln1mz / 27d0
*
      AS2ggH_L = TR * ( CF * A01 + CA * A02 )
*
*     In case the the MSbar for the heavy querk masses is chosen
*
      if(mass_scheme.eq."MSbar")then
         h1 = 4d0 * CF
*
         AS2ggH_L = AS2ggH_L - 2d0 * h1 * AS1ggH_L()
      endif
*
      return
      end
*
************************************************************************
*
*     Time-like matching conditions (FFs)
*
*     Reference: hep-ph/0504192
*
************************************************************************
*     Eq. (15)
************************************************************************
      function AS1HgT(z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      double precision z
**
*     Internal Variables
*
      double precision B0
**
*     Output Variables
*
      double precision AS1HgT
*
      B0 = ( 1d0 + ( 1d0 - z )**2 ) * ( - 1d0 - 2d0 * dlog(z) ) / z
*
      AS1HgT = 2d0 * CF * B0
*
      return
      end
*
************************************************************************
*                                                                      *
*     Additional functions required for matching nf to nf+1 flavours   *
*     at mass scales *not* equal to the mass of the heavy quark.       *
*     They are all proportional to at least power of ln(m2ph/m2th).    * 
*     (Thanks to Andrew Papanastasiou and Marco Bonvini)               *
*                                                                      *
************************************************************************
*     Eq. (B.2)                                                        *
************************************************************************
      function AS1Hg_mass(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision B0
      double precision lnk
**
*     Output Variables
*
      double precision AS1Hg_mass
*
      lnk = - dlog(k2th(nf))
*
      B0 = - 4d0 * ( z**2 + ( 1d0 - z )**2 ) * lnk
*
      AS1Hg_mass = TR * B0
*
      return
      end
*
************************************************************************
*     Eq. (B.6)                                                        *
************************************************************************
      function AS1ggH_mass_L(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
*
*     Internal Variables
*
      double precision B0
      double precision lnk
**
*     Output Variables
*
      double precision AS1ggH_mass_L
*
      lnk  = - dlog(k2th(nf))
*
      B0 = 4d0 * lnk / 3d0
*
      AS1ggH_mass_L = TR * B0
*
      return
      end
*
************************************************************************
*     Eq. (B.1)                                                        *
************************************************************************
      function APS2Hq_mass(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision lnz,lnz2,z2
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision APS2Hq_mass
*
      lnz   = dlog(z)
      lnz2  = lnz * lnz 
      z2    = z * z
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = 8d0 * ( 1d0 + z ) * lnz2
     1      - ( 8d0 + 40d0 * z + 64d0 * z2 / 3d0 ) * lnz
     2      - 160d0 / 9d0 / z + 16d0 - 48d0 * z + 448d0 * z2 / 9d0
      omeL1 = omeL1 * lnk
*
      omeL2 = - 8d0 * ( 1d0 + z ) * lnz - 16d0 / 3d0 / z - 4d0 + 4d0 * z
     1      + 16d0 * z2 / 3d0
      omeL2 = omeL2 * lnk2
*
      APS2Hq_mass = CF * TR * ( omeL2 + omeL1 )
*
      return
      end
*
************************************************************************
*     Eq. (B.3)                                                        *
************************************************************************
      function AS2Hg_mass(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/consts.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision z2
      double precision S11mz
      double precision lnz,lnz2
      double precision ln1mz,ln1mz2,ln1pz
      double precision ddilog
      double precision lnk,lnk2
      double precision omeL1,omeL1p1,omeL1p2
      double precision omeL2,omeL2p1,omeL2p2,omeL2p3
**
*     Output Variables
*
      double precision AS2Hg_mass
*
*     some definitions
*
      z2     = z * z
      lnz    = dlog(z)
      lnz2   = lnz * lnz 
      ln1mz  = dlog(1d0 - z)
      ln1mz2 = ln1mz * ln1mz 
      ln1pz  = dlog(1d0 + z)
      S11mz  = ddilog(-z)
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1p1 = CF * TR * ( ( 8d0 - 16d0 * z + 16d0 * z2 )
     1        * ( 2d0 * lnz * ln1mz - ln1mz2 + 2d0 * zeta2 ) 
     2        - ( 4d0 - 8d0 * z + 16d0 * z2 ) * lnz2
     3        - 32d0 * z * ( 1d0 - z ) * ln1mz
     4        - ( 12d0 - 16d0 * z + 32d0 * z2 ) * lnz - 56d0 + 116d0 * z
     5        - 80d0 * z2 )
*
      omeL1p2 = CA * TR * ( ( 16d0 + 32d0 * z + 32d0 * z2 )
     1        * ( S11mz + lnz * ln1pz ) + ( 8d0 - 16d0 * z + 16d0 * z2)
     2        * ln1mz2 + ( 8d0 + 16d0 * z ) * lnz2
     3        + 32d0 * z * zeta2 + 32d0 * z * ( 1d0 - z ) * ln1mz 
     4        - ( 8d0 + 64d0 * z + 352d0 * z2 / 3d0 ) * lnz
     5        - 160d0 / 9d0 / z + 16d0 - 200d0 * z + 1744d0 * z2 / 9d0 )
*
      omeL1 = ( omeL1p1 + omeL1p2 ) * lnk
*
      omeL2p1 = CF * TR * ( ( 8d0 - 16d0 * z + 16d0 * z2 ) * ln1mz
     1        - ( 4d0 - 8d0 * z + 16d0 * z2 ) * lnz
     2        - ( 2d0 - 8d0 * z ) ) 
      omeL2p2 = CA * TR * ( - ( 8d0 - 16d0 * z + 16d0 * z2 ) * ln1mz
     1        - ( 8d0 + 32d0 * z ) * lnz - 16d0 / 3d0 / z - 4d0
     2        - 32d0 * z + 124d0 * z2 / 3d0 )
      omeL2p3 = TR * TR * ( - 16d0 * ( z2 + ( 1d0 - z )**2 ) / 3d0 )
*
      omeL2 = ( omeL2p1 + omeL2p2 + omeL2p3 ) * lnk2
*
      AS2Hg_mass = omeL2 + omeL1
*
      return
      end
*
************************************************************************
*     Eq. (B.4) Regular                                                *
************************************************************************
      function ANS2qqH_mass_R(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision z2,lnz,lnz2,ln1mz
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision ANS2qqH_mass_R
*
*     some definitions
*
      z2    = z * z
      lnz   = dlog(z)
      lnz2  = lnz * lnz
      ln1mz = dlog(1d0 - z)
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = 8d0 * ( 1d0 + z2 ) * lnz / 3d0 / ( 1d0 - z ) + 8d0 / 9d0 
     1      - 88d0 * z / 9d0
      omeL1 = omeL1 * lnk 
*
      omeL2 = - 4d0 / 3d0 - 4d0 * z / 3d0 
      omeL2 = omeL2 * lnk2
*
      ANS2qqH_mass_R = CF * TR * ( omeL2 + omeL1 )
*
      return
      end
*
************************************************************************
*     Eq. (B.4) Singular                                               *
************************************************************************
      function ANS2qqH_mass_S(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision ANS2qqH_mass_S
*
*     some definitions
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = 80d0 / 9d0 / ( 1d0 - z )
      omeL1 = omeL1 * lnk 
*
      omeL2 = 8d0 / 3d0 / ( 1d0 - z ) 
      omeL2 = omeL2 * lnk2
*
      ANS2qqH_mass_S = CF * TR * ( omeL2 + omeL1 )
*
      return
      end
*
************************************************************************
*     Eq. (B.4) Local                                                  *
************************************************************************
      function ANS2qqH_mass_L(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/consts.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision ln1mz
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision ANS2qqH_mass_L
*
*     some definitions
*
      ln1mz = dlog(1d0 - z)
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = 80d0 * ln1mz / 9d0 + 16d0 * zeta2 / 3d0 + 2d0 / 3d0
      omeL1 = omeL1 * lnk 
*
      omeL2 = 8d0 * ln1mz / 3d0 + 2d0  
      omeL2 = omeL2 * lnk2
*
      ANS2qqH_mass_L = CF * TR * ( omeL2 + omeL1 )
*
      return
      end
*
************************************************************************
*     Eq. (B.5)                                                        *
************************************************************************
      function AS2gqH_mass(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision ln1mz
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision AS2gqH_mass
*
*     some definitions
*
      ln1mz = dlog(1d0 - z)
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = 160d0 / 9d0 / z - 160d0 / 9d0 + 128d0 * z / 9d0 
     1      + ( 32d0 / 3d0 / z - 32d0 / 3d0 + 16d0 * z / 3d0 ) * ln1mz
      omeL1 = omeL1 * lnk
*
      omeL2 = 16d0 / 3d0 / z - 16d0 / 3d0 + 8d0 * z / 3d0 
      omeL2 = omeL2 * lnk2
*
      AS2gqH_mass = CF * TR * ( omeL2 + omeL1 )
*
      return
      end
*
************************************************************************
*     Eq. (B.7) Regular                                                *
************************************************************************
      function AS2ggH_mass_R(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision z2,lnz,lnz2,ln1mz
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision AS2ggH_mass_R
*
*     some definitions
*
      z2    = z * z
      lnz   = dlog(z)
      lnz2  = lnz * lnz; 
      ln1mz = log(1d0 - z)
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = CF * TR * ( 8d0 * ( 1d0 + z ) * lnz2 + ( 24d0 + 40d0 * z )
     1      * lnz - 16d0 / 3d0 / z + 64d0 - 32d0 * z - 80d0 * z2 / 3d0 )
     2      + CA * TR * ( 16d0 * ( 1d0 + z ) * lnz / 3d0
     3      + 184d0 / 9d0 / z - 232d0 / 9d0 + 152d0 * z / 9d0
     4      - 184d0 * z2 / 9d0 )
      omeL1 = omeL1 * lnk 
*
      omeL2 = CF * TR * ( 8d0 * ( 1d0 + z ) * lnz + 16d0 / 3d0 / z + 4d0
     1      - 4d0 * z - 16d0 * z2 / 3d0 )
     2      + CA * TR * ( 8d0 / 3d0 / z - 16d0 / 3d0 + 8d0 * z / 3d0
     3      - 8d0 * z2 / 3d0 )
      omeL2 = omeL2 * lnk2
*
      AS2ggH_mass_R = omeL2 + omeL1
*
      return
      end
*
************************************************************************
*     Eq. (B.7) Singular                                               *
************************************************************************
      function AS2ggH_mass_S(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision AS2ggH_mass_S
*
*     some definitions
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = 80d0 / 9d0 / ( 1d0 - z )
      omeL1 = omeL1 * lnk 
*
      omeL2 = 8d0 / 3d0 / ( 1d0 - z )
      omeL2 = omeL2 * lnk2
*
      AS2ggH_mass_S = CA * TR * ( omeL2 + omeL1 )
*
      return
      end
*
************************************************************************
*     Eq. (B.7) Local                                                  *
************************************************************************
      function AS2ggH_mass_L(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision ln1mz
      double precision lnk,lnk2
      double precision omeL1,omeL2
**
*     Output Variables
*
      double precision AS2ggH_mass_L
*
*     some definitions
*
      ln1mz = log(1d0 - z)
*
      lnk  = - dlog(k2th(nf))
      lnk2 = lnk * lnk
*
      omeL1 = CF * TR * ( 4d0 ) 
     1      + CA * TR * ( 16d0 / 3d0 + 80d0 * ln1mz / 9d0 )
      omeL1 = omeL1 * lnk 
*
      omeL2 = TR * TR * ( 16d0 / 9d0 ) 
     1      + CA * TR * ( 8d0 * ln1mz / 3d0 )
      omeL2 = omeL2 * lnk2
*
      AS2ggH_mass_L = omeL2 + omeL1
*
      return
      end
*
************************************************************************
*
*     Time-like matching conditions (FFs)
*
*     Reference: hep-ph/0504192
*
************************************************************************
*     Eq. (15)
************************************************************************
      function AS1HgT_mass(nf,z)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision z
**
*     Internal Variables
*
      double precision lnk
      double precision B0
**
*     Output Variables
*
      double precision AS1HgT_mass
*
      lnk = - dlog(k2th(nf))
*
      B0 = ( 1d0 + ( 1d0 - z )**2 ) * ( - lnk ) / z
*
      AS1HgT_mass = 2d0 * CF * B0
*
      return
      end
