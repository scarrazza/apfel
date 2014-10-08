        FUNCTION BETA(P,Q)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p    --- Parameter  ( p > 0 )
C                q    --- Parameter  ( q > 0 )
C       Output:  BETA --- B(p,q)
C       Routine called: GAMMA for computing â(x)
C       ==========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        BETA=GAMMA(P)*GAMMA(Q)/GAMMA(P+Q)
        RETURN
        END
