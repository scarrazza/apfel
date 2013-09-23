****************************************************************
*
*     a_QED.f:
*
*     It returns a_QED = alpha_QED / 4 / pi .
*
*     Note that the only argument of this routine is the 
*     factorization scale mu2F which is convereted into the 
*     renormalization scale inside the routine itself.
*
****************************************************************
      function a_QED(mu2f)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha_ref_QED.h"
      include "../commons/kren.h"
      include "../commons/m2th.h"
      include "../commons/mass_scheme.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Evs.h"
      include "../commons/MaxFlavourAlpha.h"
**
*     Input Variables
*
      double precision mu2F
**
*     Internal Variables
*
      integer i
      integer nfi,nff
      integer dnf,snf
      double precision mur2th(4:6)
      double precision mur2,mur20
      double precision aqedi,aqedr0
      double precision alphaqedev
      external alphaqedev
**
*     Output Variables
*
      double precision a_QED
*
      aqedr0 = alpha_ref_QED / 4d0 / pi
*     Uncomment to switch off the alpha running
c      a_QED  = aqedr0
c      return
*
      mur20 = q2_ref_QED
      mur2  = kren * mu2F
      do i=4,6
         mur2th(i) = kren * m2th(i)
      enddo
*
      if(Evs.eq."FF")then
         nfi = Nf_FF
         nff = Nf_FF
      elseif(Evs.eq."VF")then
         if(mur2.ge.mur2th(6))then
            nff = 6
         elseif(mur2.ge.mur2th(5))then
            nff = 5
         elseif(mur2.ge.mur2th(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfMaxAlpha) nff = nfMaxAlpha
*
         if(mur20.gt.mur2th(6))then
            nfi = 6
         elseif(mur20.gt.mur2th(5))then
            nfi = 5
         elseif(mur20.gt.mur2th(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfMaxAlpha) nfi = nfMaxAlpha
      endif
*
 10   if(nff.eq.nfi) then
         a_QED = alphaqedev(nfi,mur2,mur20,aqedr0)
         return
      else
         if(nff.gt.nfi)then
            dnf = 1
            snf = 1
         else
            dnf = -1
            snf = 0
         endif
*
         aqedi = alphaqedev(nfi,mur2th(nfi+snf),mur20,aqedr0)
*
         aqedr0  = aqedi
         mur20 = mur2th(nfi+snf)
         nfi  = nfi + dnf
         goto 10
      endif
*
      end
*
****************************************************************
      FUNCTION ALPHAQEDEV(NF,Q2F,Q2I,ALPHAREF)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION Q2F,Q2I
      DOUBLE PRECISION ALPHAREF,BETA0QED
**
*     Internal Variables
*
      DOUBLE PRECISION BETA0
      DOUBLE PRECISION L
**
*     Output Variables
*
      DOUBLE PRECISION ALPHAQEDEV
*
      BETA0 = BETA0QED(NF)
      L = DLOG( Q2F / Q2I )
*
      ALPHAQEDEV = ALPHAREF / ( 1D0 - ALPHAREF * BETA0 * L )
*
      RETURN
      END
*
****************************************************************
      function beta0qed(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer nc
      double precision sumch2(3:6)
**
*     Output Variables
*
      double precision beta0qed
*
*     Number of colours
*
      nc = 3
*
*     Sum of the first 3, 4, 5 and 6 squared electric charges
*
      sumch2(3) = 2d0 / 3d0
      sumch2(4) = 10d0 / 9d0
      sumch2(5) = 11d0 / 9d0
      sumch2(6) = 5d0 / 3d0
*
      beta0qed = 8d0 / 3d0 * ( nc * sumch2(nf) )!+ 3d0 )
*
      return
      end
