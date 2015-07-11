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
      include "../commons/TauMass.h"
      include "../commons/mass_scheme.h"
      include "../commons/Nf_FF.h"
      include "../commons/Evs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/LeptEvol.h"
**
*     Input Variables
*
      double precision mu2F
**
*     Internal Variables
*
      integer i
      integer nfi,nff
      integer Nl_FF,nli,nlf
      integer dnf,snf
      double precision mur2th(4:6),MTau2ren
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
*     Uncomment to switch off the running
c      a_QED  = aqedr0
c      return
*
      mur20 = q2_ref_QED
      mur2  = kren * mu2F
*
      if(Evs.eq."FF")then
*     In the FFNS assume that there is no tau contribution
         Nl_FF = 0
         if(LeptEvol) Nl_FF = 2
*
         a_QED = alphaqedev(Nf_FF,Nl_FF,mur2,mur20,aqedr0)
         return
      elseif(Evs.eq."VF")then
         do i=4,6
            mur2th(i) = kren * m2th(i)
         enddo
*
         nli = 0
         nlf = 0
         Mtau2ren = kren * MTau * MTau
         if(LeptEvol)then
            nli = 2
            if(mur20.gt.Mtau2ren) nli = 3
            nlf = 2
            if(mur2.ge.Mtau2ren) nlf = 3
         endif
*
*     If the final and initial and numbers are equal ...
*
         if(nli.eq.nlf)then
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
*
 10         if(nff.eq.nfi)then
               a_QED = alphaqedev(nfi,nli,mur2,mur20,aqedr0)
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
               aqedi = alphaqedev(nfi,nli,mur2th(nfi+snf),mur20,aqedr0)
*
               aqedr0  = aqedi
               mur20 = mur2th(nfi+snf)
               nfi  = nfi + dnf
               goto 10
            endif
*
*     ... if else the final and initial and numbers are not equal ...
*
         else
*
*     ...first evolve from the initial scale to the tau scale with "nli"
*     active leptons and the from the tau scale to the final scale with
*     "nlf" active leptons.
*
            mur20 = q2_ref_QED
            mur2  = Mtau2ren
 12         if(mur2.ge.mur2th(6))then
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
*
 11         if(nff.eq.nfi)then
               a_QED = alphaqedev(nfi,nli,mur2,mur20,aqedr0)
               if(nli.eq.nlf)then
                  return
               else
                  nli    = nlf
                  mur20  = Mtau2ren
                  mur2   = kren * mu2F
                  aqedr0 = a_QED
                  goto 12
               endif
            else
               if(nff.gt.nfi)then
                  dnf = 1
                  snf = 1
               else
                  dnf = -1
                  snf = 0
               endif
*
               aqedi = alphaqedev(nfi,nli,mur2th(nfi+snf),mur20,aqedr0)
*
               aqedr0 = aqedi
               mur20  = mur2th(nfi+snf)
               nfi    = nfi + dnf
               goto 11
            endif
         endif
      endif
*
      end
*
****************************************************************
      FUNCTION ALPHAQEDEV(NF,NL,Q2F,Q2I,ALPHAREF)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER NF,NL
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
      BETA0 = BETA0QED(NF,NL)
      L = DLOG( Q2F / Q2I )
*
      ALPHAQEDEV = ALPHAREF / ( 1D0 - ALPHAREF * BETA0 * L )
*
      RETURN
      END
*
****************************************************************
      function beta0qed(nf,nl)
*
      implicit none
**
*     Input Variables
*
      integer nf,nl
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
      beta0qed = 8d0 / 3d0 * ( nc * sumch2(nf) + nl )
*
      return
      end
*
****************************************************************************
      function a_as_exact(nf,nl,as0,a0,as,ipt)
*
      implicit none
**
*     input variables
*
      integer nf,nl,ipt
      double precision a0,as0
      double precision as
**
*     inernal variables
*
      integer nstep,n
      double precision a
      double precision h
      double precision k1,k2,k3,k4
      double precision fbetaQED,fbeta

      parameter(nstep=10)
**
*     output variables
*
      double precision a_as_exact
*
      a = a0
      h = ( as - as0 ) / nstep
*
*     Fourth-order runge-kutta beyond the leading order
*
      do n=1,nstep
         k1 = h * fbetaQED(a,nf,nl) / fbeta(as,nf,ipt)
         k2 = h * fbetaQED(a+k1/2d0,nf,nl) / fbeta(as+h/2d0,nf,ipt)
         k3 = h * fbetaQED(a+k2/2d0,nf,nl) / fbeta(as+h/2d0,nf,ipt)
         k4 = h * fbetaQED(a+k3,nf,nl) / fbeta(as+h,nf,ipt)
*
         a  = a + ( k1 + 2d0 * k2 + 2d0 * k3 + k4 ) / 6d0
         as = as + h
      enddo
*
      a_as_exact = a
*
      return
      end
*
****************************************************************************
*
*     QED beta function.
*
****************************************************************************
      function fbetaQED(a,nf,nl)
*
      implicit none
**
*     Input Variables
*
      double precision a
      integer nf,nl
**
*     Internal Variables
*
      double precision beta0qed
**
*     Output Variables
*
      double precision fbetaQED
*
      fbetaQED = a**2d0 * beta0qed(nf,nl)
*
      return
      end
