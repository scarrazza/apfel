************************************************************************
*
*     SetLambdaQCDRef.f:
*
*     This subroutine sets the reference value of LambdaQCD for a
*     given number of flavours.
*
************************************************************************
      subroutine SetLambdaQCDRef(lambdaref,nref)
*
      implicit none
*
      include "../commons/lambda_ref_QCD.h"
*
*     Variables
*
      integer i
      integer nref
      double precision lambdaref
*
      lambda_ref_QCD = lambdaref
      n_ref_QCD      = nref
      InLambdaQCD    = "done"
*
*     Initialize values of LambdaQCD for all the number of flavours
*     at the reference value.
*
      do i=3,6
         LambdaQCD(i) = lambda_ref_QCD
      enddo
*
      return
      end
