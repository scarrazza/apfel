************************************************************************
*
*     SetScaleVariationProcedure.f:
*
*     This subroutine sets the the procedure to be used to vary
*     factorisation and remormalisation scale:
*     - 0: consistent scale variation in DIS and evolution 
*     - 1: variation only in the DIS structure functions
*
************************************************************************
      subroutine SetScaleVariationProcedure(svproc)
*
      implicit none
*
      include "../commons/ScaleVariationProcedure.h"
*
*     Variables
*
      integer svproc
*
      ScVarProc = svproc
      InScVarProc = "done"
*
      return
      end
