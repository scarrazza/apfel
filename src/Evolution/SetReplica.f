************************************************************************
*
*     SetReplica.f:
*
*     This subroutine sets the replica to be used as intial PDFs.
*     Used only if an LHAPDF set hes been chosen.
*
************************************************************************
      subroutine SetReplica(nr)
*
      implicit none
*
      include "../commons/Replica.h"
*
*     Variables
*
      integer nr
*
      irep  = nr
      InRep = "done"
*
      return
      end
