************************************************************************
*
*     LockGrids.f:
*
*     This subroutine enables or disables the the locking of the internal
*     subgrids (if more than one).
*
************************************************************************
      subroutine LockGrids(lg)
*
      implicit none
*
      include "../commons/lock.h"
*
      logical lg
*
      lock   = lg
      InLock = "done"
*
      return
      end
