************************************************************************
*
*     Hellx.f:
*
*     This function returns the small-x resummed singlet splitting
*     functions.
*     The inputs are:
*     - k   = 4: qq, 5: qg, 6: gq, 7: gg (splittinf function),
*     - ipt = perturbative order to be matched to
*     - as  = value of alphas / ( 4 * pi ),
*     - y   = Bjorken's variable.
*
************************************************************************
      function Hellx(k,ipt,as,y)
*
      implicit none
**
*     Input Variables
*
      integer k
      integer ipt
      double precision as
      double precision y
**
*     Internal Variables
*
      double precision Hellx
**
*     Input Variables
*
      Hellx = 0d0
*
      return
      end
