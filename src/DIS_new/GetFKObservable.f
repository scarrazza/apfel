************************************************************************
*
*     GetFKObservable.f:
*
*     This subroutine returns the observable set using SetFKObservable.
*
************************************************************************
      subroutine GetFKObservable()
*
      implicit none
*
      include "../commons/FKObservable.h"
*
      write(6,*) FKObservable
*
      return
      end
