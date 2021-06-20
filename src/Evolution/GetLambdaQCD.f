************************************************************************
*
*     GetLambdaQCD.f:
*
*     This function returns the value of LambdaQCD with i (i = 3, 4 ,5,
*     6) active flavours.
*
************************************************************************
      function GetLambdaQCD(i)
*
      implicit none
*
      include "../commons/lambda_ref_QCD.h"
**
*     Input Variables
*
      integer i
**
*     Output Variables
*
      double precision GetLambdaQCD
*
*     Check the consistency of the input
*
      if(i.lt.3.and.i.gt.6)then
         write(6,*) "In GetLambdaQCD.f:"
         write(6,*) "Invalid quark index i =",i
         call exit(-10)
      endif
*
      GetLambdaQCD = LambdaQCD(i)
*
      return
      end
