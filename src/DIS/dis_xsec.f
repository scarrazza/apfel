************************************************************************
*
*     dis_xsec.f:
*
*     Driver that calls nc_dis and cc_dis.
*
************************************************************************
      SUBROUTINE DIS_XSEC(X,QI,QF,Y,POL,PROC,SCHEME,PTO,PDFSET,IREP,
     1                    TARGET,PROJ,F2,F3,FL,SIGMA)
*
      IMPLICIT NONE
**
*     Input Varibles
*
      INTEGER PTO,IREP
      DOUBLE PRECISION X,QI,QF,Y,POL
      CHARACTER*2  PROC
      CHARACTER*5  SCHEME
      CHARACTER*53 PDFSET
      CHARACTER*9  TARGET
      CHARACTER*12 PROJ
**
*     Output Variables
*
      DOUBLE PRECISION F2(3:7),F3(3:7),FL(3:7),SIGMA(3:7)
*
      IF(PROC(1:2).EQ."EM".OR.PROC(1:2).EQ."NC")THEN
         CALL NC_DIS(X,QI,QF,Y,POL,PROC,SCHEME,PTO,PDFSET,IREP,TARGET,
     1               PROJ,F2,F3,FL,SIGMA)
      ELSEIF(PROC(1:2).EQ."CC")THEN
         CALL CC_DIS(X,QI,QF,Y,POL,SCHEME,PTO,PDFSET,IREP,TARGET,
     1               PROJ,F2,F3,FL,SIGMA)
      ELSE
         WRITE(6,*) "In dis_xsec.f:"
         WRITE(6,*) "Invalid process, PROC =",PROC
         CALL EXIT(-10)
      ENDIF
*
      CALL CleanUp
*
      RETURN
      END
