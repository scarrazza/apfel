*     -*-fortran-*-
*
*     Switch for the FONLL damping factor
*
      logical DampingFONLL
      character*4 InDampingFONLL
*
      common / DampingFONLLAPFEL / DampingFONLL,InDampingFONLL
*
*     Power of the damping
*
      integer DampPowerFONLL(4:6)
      character*4 InDampPowerFONLL
*
      common / DampPowerFONLLAPFEL / DampPowerFONLL,InDampPowerFONLL
