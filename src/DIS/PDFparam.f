************************************************************************
*
*     PDFparam.f:
*
*     Initialize all the parameters depending on PDFs, i.e. PDF set,
*     alphas and heavy quark masses.
*
************************************************************************
      SUBROUTINE PDFPARAM(PDFSET,IREP)
*
      IMPLICIT NONE
*
      include "../commons/scale.h"
      include "../commons/heavythr.h"
      include "../commons/asq.h"
      include "../commons/nf.h"
      include "../commons/vfns.h"
      include "../commons/pdf.h"
      include "../commons/jpt.h"
**
*     Input Varibles
*
      INTEGER IREP
      CHARACTER*53 PDFSET
**
*     Internal Variables
*
      INTEGER I,NREP
      DOUBLE PRECISION alphasPDF,AlphaQCD
      DOUBLE PRECISION HeavyQuarkMass
      DOUBLE PRECISION QTH(4:6)
      DOUBLE PRECISION PI
      PARAMETER(PI = 3.14159265358979D0)
*
      IF(PDFSET(1:5).EQ."APFEL")THEN
         PDF = "APF"
*
*     Initialize APFEL
*
         call EnableWelcomeMessage(.false.)
         call SetPerturbativeOrder(ipt)
         IF(VFNS.EQ."FFNS".OR.VFNS.EQ."FFN0")THEN
            call SetFFNS(3)
         ELSEIF(VFNS.EQ."ZMVN".OR.VFNS.EQ."FONLL")THEN
            call SetVFNS
         ENDIF
*
         call InitializeAPFEL
*
*     Evolve PDFs with APFEL
*
         call EvolveAPFEL(Q0,Q)
*
*     Get Value of alphas
*
         call SetReplica(IREP)
         ASQ = AlphaQCD(Q) / 4D0 / PI
*
*     Get Heavy quark masses
*
         Q2TH(4) = HeavyQuarkMass(4,Q)**2D0
         Q2TH(5) = HeavyQuarkMass(5,Q)**2D0
         Q2TH(6) = HeavyQuarkMass(6,Q)**2D0
*
      ELSE
         PDF = "LHA"
*
*     Initialize PDF set with LHAPDF
*
         CALL InitPDFsetbyName(PDFSET)
         CALL numberPDF(NREP)
         CALL InitPDF(IREP)
*
*     Get alphas from LHAPDF
*
         ASQ = alphasPDF(Q) / 4D0 / PI
*
*     Get heavy quark thresholds from LHAPDF
*
         DO I=4,6
            CALL GetQmass(I,QTH(I))
            Q2TH(I) = QTH(I)**2D0
         ENDDO
      ENDIF
*
*     Count number of active flavours
*
      NF = 3
      DO I=4,6
         IF((Q*Q).GT.Q2TH(I)) NF = I
      ENDDO
*
*     In case one the FFNSs is chosen, set by hand the values
*     of the heavy quark masses to be used in the coefficient
*     functions.
*
      IF(VFNS.EQ."FFNS".OR.VFNS.EQ."FFN0")THEN
         Q2TH(4) = 2D0
         Q2TH(5) = 1D10
         Q2TH(6) = 1D10
      ENDIF
*
      RETURN
      END
