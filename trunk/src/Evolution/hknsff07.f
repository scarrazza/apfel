C *********************************************************************
C  hknsff07.f  Version 1.2                                2009/DEC/19
C
C  [Package for the HKNS fragmentaion functions]
C  Reference:
C    Determination of fragmentation functions and their uncertainties
C    M. Hirai, S. Kumano, T.-H. Nagai, and K. Sudoh
C    Phys. Rev. D75 (2007) 094009. (hep-ph/0702250) 
C *********************************************************************
C ---------------------------------------------------------------------
C  SUBROUTINE HKNSFF(Q2,X,ISET,ICHARGE,FF,GRAD):
C
C   Subroutine HKNSFF returns the values of fragmentaion functions
C   and their gradient terms at specified Q^2 and x point 
C   by interpolating the grid data.
C   [ Log(Q^2): LINEAR INTERPOLATION, x: CUBIC SPLINE INTERPOLATION ]
C
C   INPUT:
C     Q2, X ... Q^2 and x values at which the functions are calculated. 
C               Available range: 10^-2 <= X <= 1.0,
C                           1.0 GeV^2 <= Q^2 <= 10^8 GeV^2.
C     ISET=1: Pion LO Fragmentation functions and their gradient terms 
C          2: Pion NLO
C          3: Kaon LO
C          4: Kaon NLO
C          5: Proton LO
C          6: Proton NLO
C
C     ICHARGE=1: pi^+, K^+, or proton
C     ICHARGE=2: pi^-, K^-, or neutron
C     ICHARGE=3: pi^0=[pi^+ + pi^-]/2, [K^0+K^0b]/2, or [p+pb]/2
C        If you want to obtain the fragmentaion functions for each
C        K0, K0b, pb, or nb, you may use the relations in Appendix of
C        hep-ph/0702250.
C
C   OUTPUT: Arrays FF(-5:5) & GRAD(I,J)
C
C     FF(I) --> HKNS fragmentation functions (FFs).
C      I = -5 ... b-bar quark (D_b-bar = D_b)
C          -4 ... c-bar quark (D_c-bar = D_c)
C          -3 ... s-bar quark
C          -2 ... d-bar quark 
C          -1 ... u-bar quark 
C           0 ... gluon D_g(x)
C           1 ... u quark 
C           2 ... d quark 
C           3 ... s quark 
C           4 ... c quark 
C           5 ... b quark 
C
C     GRAD(I,J) --> Gradient terms of HKSN FFs
C      I is the same index as the one in FF(I).
C      J indicates the parameter index for a gradient term dFF(I)/da_J
C      (a_J = parameter).
C
C Pion,J= 1..2: g (2ndM, alpha)           2ndM = second moment
C         3..5: u (2ndM, alpha, beta)
C         6..8: d (2ndM, alpha, beta)
C        9..11: c (2ndM, alpha, beta)
C       12..14: b (2ndM, alpha, beta)
C   For example, the above J=14 indicates d D_b^{pi^+}/d beta_b^{pi^+}.
C
C Kaon,J= 1..2: g (2ndM, beta)
C         3..5: u (2ndM, alpha, beta)
C         6..8: d (2ndM, alpha, beta)
C        9..11: sb(2ndM, alpha, beta)
C       12..14: c (2ndM, alpha, beta)
C       15..17: b (2ndM, alpha, beta)
C
CProton,J=1..2: g (2ndM, beta)
C         3..5: u (2nsM, alpha, beta)
C         6..8: qb(2ndM, alpha, beta)
C        9..11: c (2ndM, alpha, beta)
C       11..13: b (2ndM, alpha, beta)
C
C   NOTE: The returned values are not multiplied by x.
C
C      *  Error matrix can be used by declaring a common block:
C         COMMON/ERRM/EM(17,17). This matrix is defined as
C         the inverse matrix of Hessian multiplied by Delta chi^2:
C         EM(i,j)=Delta chi^2*H_ij^-1.
C         The values of Delta chi^2 are as follows:
C         15.9359730(pion), 19.1977555(kaon), 14.8470228(proton).
C *********************************************************************
      SUBROUTINE HKNSFF(Q2,X,ISET,ICHARGE,FF,GRAD)
C ---------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NQ=33, NX=53, ND=146, NFF=8, NSET=3, IFILE=20)
      DIMENSION IREAD(2), QG(NQ), XG(NX),PDFJ1(ND), PDFJ2(ND)
     +         ,FF(-5:5), GRADFF(8,17), GRAD(-5:5,17)
     +         ,BXG(NX,NQ,ND), CXG(NX,NQ,ND), DXG(NX,NQ,ND)
     +         ,PDFG(NX,NQ,ND), EMI(17,17,6)

      COMMON/ERRM/EM(17,17)
      SAVE IREAD, NPAR, BXG, CXG, DXG, PDFG 
      DATA IREAD /1, 1/ 
      DATA CTHRE,BTHRE/2.0449D0,18.490D0/
      DATA EPS/1.D-12/

      INCLUDE '../FragFunc/HKNS/EM.inc'
      INCLUDE '../FragFunc/HKNS/hkns_pilo.inc'
      INCLUDE '../FragFunc/HKNS/hkns_pinlo.inc'
      INCLUDE '../FragFunc/HKNS/hkns_klo.inc'
      INCLUDE '../FragFunc/HKNS/hkns_knlo.inc'
      INCLUDE '../FragFunc/HKNS/hkns_plo.inc'
      INCLUDE '../FragFunc/HKNS/hkns_pnlo.inc'

C Q2 AND X GRID.
      DATA QG /
     +  1.000000D+00, 1.467799D+00, 2.154435D+00,
     +  3.162278D+00, 4.641589D+00, 6.812921D+00,
     +  1.000000D+01, 1.467799D+01, 2.154435D+01,
     +  3.162278D+01, 4.641589D+01, 6.812921D+01,
     +  1.000000D+02, 1.778279D+02, 3.162278D+02, 5.623413D+02,
     +  1.000000D+03, 1.778279D+03, 3.162278D+03, 5.623413D+03,
     +  1.000000D+04, 1.778279D+04, 3.162278D+04, 5.623413D+04,
     +  1.000000D+05, 1.778279D+05, 3.162278D+05, 5.623413D+05,
     +  1.000000D+06, 4.641589D+06, 
     +  1.000000D+07, 4.641589D+07,  
     +  1.000000D+08  /

      DATA XG / 
     +  1.000000D-02, 1.154782D-02, 1.333521D-02, 1.539927D-02,
     +  1.778279D-02, 2.053525D-02, 2.371374D-02, 2.738420D-02,
     +  3.162278D-02, 3.651741D-02, 4.216965D-02, 4.869675D-02,
     +  5.623413D-02, 6.493816D-02, 7.498942D-02, 8.659643D-02,
     +  1.000000D-1, 1.250000D-1, 1.500000D-1, 1.750000D-1,
     +  2.000000D-1, 2.250000D-1, 2.500000D-1, 2.750000D-1,
     +  3.000000D-1, 3.250000D-1, 3.500000D-1, 3.750000D-1,
     +  4.000000D-1, 4.250000D-1, 4.500000D-1, 4.750000D-1, 
     +  5.000000D-1, 5.250000D-1, 5.500000D-1, 5.750000D-1,
     +  6.000000D-1, 6.250000D-1, 6.500000D-1, 6.750000D-1,
     +  7.000000D-1, 7.250000D-1, 7.500000D-1, 7.750000D-1,
     +  8.000000D-1, 8.250000D-1, 8.500000D-1, 8.750000D-1,
     +  9.000000D-1, 9.250000D-1, 9.500000D-1, 9.750000D-1,
     +  1.000000D+0 /

C CALCULATE SPLINE COEFFICIENTS.
      IF((IREAD(1).NE.1).AND.(IREAD(2).EQ.ISET)) GO TO 20
      CTHRE=CTHRE-EPS
      BTHRE=BTHRE-EPS

C READ GRID DATA AND CALCULATE SPLINE COEFFICIENTS.
      IF(ISET.EQ.1) THEN
C        OPEN(UNIT=IFILE,FILE='hkns_pilo.grd',STATUS='OLD')
C        OPEN(UNIT=IFILE+1,FILE='gradpilo.grd',STATUS='OLD')
        NPAR=14
        DO J=1,NQ
          DO K=1,NX-1
            DO I=1,8
              PDFG(K,J,I) = PDFG_PILO(K,J,I)
            ENDDO
            DO NR=1,NPAR
              NI=10+NFF*(NR-1)
              DO I=NI,NI+NFF-1
                PDFG(K,J,I) = PDFG_PILO(K,J,I)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISET.EQ.2) THEN
C        OPEN(UNIT=IFILE,FILE='hkns_pinlo.grd',STATUS='OLD')
C        OPEN(UNIT=IFILE+1,FILE='gradpinlo.grd',STATUS='OLD')
        NPAR=14
        DO J=1,NQ
          DO K=1,NX-1
            DO I=1,8
              PDFG(K,J,I) = PDFG_PINLO(K,J,I)
            ENDDO
            DO NR=1,NPAR
              NI=10+NFF*(NR-1)
              DO I=NI,NI+NFF-1
                PDFG(K,J,I) = PDFG_PINLO(K,J,I)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISET.EQ.3) THEN
C        OPEN(UNIT=IFILE,FILE='hkns_klo.grd',STATUS='OLD')
C        OPEN(UNIT=IFILE+1,FILE='gradklo.grd',STATUS='OLD')
        NPAR=17
        DO J=1,NQ
          DO K=1,NX-1
            DO I=1,8
              PDFG(K,J,I) = PDFG_KLO(K,J,I)
            ENDDO
            DO NR=1,NPAR
              NI=10+NFF*(NR-1)
              DO I=NI,NI+NFF-1
                PDFG(K,J,I) = PDFG_KLO(K,J,I)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISET.EQ.4) THEN
c        OPEN(UNIT=IFILE,FILE='hkns_knlo.grd',STATUS='OLD')
c        OPEN(UNIT=IFILE+1,FILE='gradknlo.grd',STATUS='OLD')
        NPAR=17
        DO J=1,NQ
          DO K=1,NX-1
            DO I=1,8
              PDFG(K,J,I) = PDFG_KNLO(K,J,I)
            ENDDO
            DO NR=1,NPAR
              NI=10+NFF*(NR-1)
              DO I=NI,NI+NFF-1
                PDFG(K,J,I) = PDFG_KNLO(K,J,I)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISET.EQ.5) THEN
C        OPEN(UNIT=IFILE,FILE='hkns_plo.grd',STATUS='OLD')
C        OPEN(UNIT=IFILE+1,FILE='gradplo.grd',STATUS='OLD')
        NPAR=13
        DO J=1,NQ
          DO K=1,NX-1
            DO I=1,8
              PDFG(K,J,I) = PDFG_PLO(K,J,I)
            ENDDO
            DO NR=1,NPAR
              NI=10+NFF*(NR-1)
              DO I=NI,NI+NFF-1
                PDFG(K,J,I) = PDFG_PLO(K,J,I)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISET.EQ.6) THEN
C        OPEN(UNIT=IFILE,FILE='hkns_pnlo.grd',STATUS='OLD')
C        OPEN(UNIT=IFILE+1,FILE='gradpnlo.grd',STATUS='OLD')
        NPAR=13
        DO J=1,NQ
          DO K=1,NX-1
            DO I=1,8
              PDFG(K,J,I) = PDFG_PNLO(K,J,I)
            ENDDO
            DO NR=1,NPAR
              NI=10+NFF*(NR-1)
              DO I=NI,NI+NFF-1
                PDFG(K,J,I) = PDFG_PNLO(K,J,I)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        WRITE(*,1010) ISET
 1010   FORMAT(' ','HKNSFF ERROR: ISET =', I3)
        STOP
      END IF

      DO I=1,17
        DO J=1,17
          EM(I,J)=EMI(I,J,ISET)
        ENDDO
      ENDDO

C      DO J=1,NQ
C        DO K=1,NX-1
C          READ(IFILE,1025) (PDFG(K,J,I), I=1,8)
C          DO NR=1,NPAR
C            NI=10+NFF*(NR-1)
C            READ(IFILE+1,1025) (PDFG(K,J,I),I=NI,NI+NFF-1)
C          ENDDO
C        ENDDO
C      ENDDO
C 1025 FORMAT(1X,8(1PE14.6))
C      CLOSE(IFILE)
C      CLOSE(IFILE+1)

      DO I=1,ND
        DO J=1,NQ
          PDFG(NX,J,I)=0.D0 ! x=1 FF=0.D0
          CALL LSPLINE(NX,XG,PDFG,BXG,CXG,DXG,ISET,I,J)
        ENDDO
      ENDDO

      IREAD(1)=2
      IREAD(2)=ISET
   20 CONTINUE

c$$$      OPEN(UNIT=12,FILE="hkns_pilo.inc",STATUS="UNKNOWN")
c$$$      WRITE(12,*) "     DIMENSION PDFG_PILO(NX,NQ,ND)"
c$$$      WRITE(12,*) "      "
c$$$      DO J=1,NQ
c$$$        DO K=1,NX-1
c$$$          WRITE(12,1111) "      DATA (PDFG_PILO(",K,",",J,
c$$$     1                   ",I),I=1,8)/"
c$$$          WRITE(12,1112) "     &",(PDFG(K,J,I),I=1,4)
c$$$          WRITE(12,1122) "     &",(PDFG(K,J,I),I=5,8)
c$$$
c$$$          DO NR=1,NPAR
c$$$             WRITE(12,1113) "      DATA (PDFG_PILO(",K,",",J,
c$$$     1                      ",I),I=",NI,",",NI+NFF-1,")/"
c$$$             NI=10+NFF*(NR-1)
c$$$             WRITE(12,1112) "     &",(PDFG(K,J,I),I=NI,NI+4-1)
c$$$             WRITE(12,1122) "     &",(PDFG(K,J,I),I=NI+4,NI+NFF-1)
c$$$          ENDDO
c$$$
c$$$        ENDDO
c$$$      ENDDO
c$$$      CLOSE(12)
c$$$ 1111 FORMAT(A,I2,A,I2,A)
c$$$ 1112 FORMAT(A,4(1PE14.6,","))
c$$$ 1122 FORMAT(A,3(1PE14.6,","),1PE14.6,"/")
c$$$ 1113 FORMAT(A,I2,A,I2,A,I3,A,I3,A)
c$$$      STOP

      DO I=1,11
        FF(I-6)=0.D0
      END DO
      DO J=1,17
        DO I=1,8
          GRADFF(I,J)=0.D0
        ENDDO
      ENDDO

C CHECK X AND Q2 VALUES.
      IF(ICHARGE.GT.3) THEN
        WRITE(*,1026) ICHARGE
 1026   FORMAT(' ','HKNSFF ERROR: ICHARGE =', I3)
        STOP
      ENDIF

      IF((X.LT.1.D-2).OR.(X.GT.1.000000001D0)) THEN
        WRITE(*,1030) X
 1030   FORMAT (' ','FF WARNING: OUT OF RANGE --> X =', 1PE12.3)
        STOP
      ENDIF
      IF((Q2.LT.1.D0).OR.(Q2.GT.1.D8)) THEN
        WRITE(*,1040) Q2
 1040   FORMAT (' ','FF WARNING: OUT OF RANGE --> Q2 =', 1PE12.3)
        STOP
      ENDIF

C INTERPOLATION.
C X: CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION.
      J=ISERCH(NQ,QG,Q2)
      IF(J.EQ.NQ) J=NQ-1
      K=ISERCH(NX,XG,X)
      DO I=1,ND
        DX=X-XG(K)
        PDFJ1(I)=PDFG(K,J,I)
     >       +DX*(BXG(K,J,I)+DX*(CXG(K,J,I)+DX*DXG(K,J,I)))
        PDFJ2(I)=PDFG(K,J+1,I)
     >       +DX*(BXG(K,J+1,I)+DX*(CXG(K,J+1,I)+DX*DXG(K,J+1,I)))
      ENDDO

C -- Fragmentation functions --
      T=(DLOG(Q2)-DLOG(QG(J)))/(DLOG(QG(J+1))-DLOG(QG(J)))
      DO I=1,3
        FF(I-1)=(1.D0-T)*PDFJ1(I)+T*PDFJ2(I)     ! g, u, d
        FF(-I)=(1.D0-T)*PDFJ1(I+3)+T*PDFJ2(I+3)  ! ub, [db, s]^pi+, [sb, s]^K; 
        FF(I+2)=(1.D0-T)*PDFJ1(I+5)+T*PDFJ2(I+5) ! s, c, b
      ENDDO
      IF(Q2.LT.CTHRE) FF(4)=0.D0
      IF(Q2.LT.BTHRE) FF(5)=0.D0
      FF(-4)=FF(4) ! cb=c
      FF(-5)=FF(5) ! bb=b

      IF((ISET.EQ.3).or.(ISET.EQ.4)) THEN  ! K+ 
        FF(-3)=FF(-2) ! FF(-2)=sb -> FF(-3)
        FF(-2)=FF(2)  ! db=d
      ENDIF

      IF(ICHARGE.EQ.2) THEN ! Negative charge
        IF(ISET.LT.3) THEN  ! pi:u<->ub, d<->db
          TEMP1=FF(-1);TEMP2=FF(2)
          FF(-1)=FF(1);FF(2)=FF(-2)
          FF(1)=TEMP1;FF(-2)=TEMP2

        ELSE IF((ISET.EQ.3).OR.(ISET.EQ.4)) THEN ! K:u<->ub, s<->sb
          TEMP1=FF(-1);TEMP2=FF(3)
          FF(-1)=FF(1);FF(3)=FF(-3)
          FF(1)=TEMP1;FF(-3)=TEMP2

        ELSE IF(ISET.GT.4) THEN ! neutron, u<->d, ub<->db
          TEMP1=FF(-1);TEMP2=FF(1)
          FF(-1)=FF(-2);FF(1)=FF(2)
          FF(-2)=TEMP1;FF(2)=TEMP2
        ENDIF

      ELSE IF(ICHARGE.EQ.3) THEN ! (pi^+ + pi^-)/2, (K^0 + K^0b)/2, (p+pb)/2
        IF((ISET.EQ.3).OR.(ISET.EQ.4)) THEN ! K^0: u<->d, ub<->db
          TEMP1=FF(-1);TEMP2=FF(1)
          FF(-1)=FF(-2);FF(1)=FF(2)
          FF(-2)=TEMP1;FF(2)=TEMP2
        ENDIF
        DO I=1,5 ! (pi^+ + pi^-)/2, (K^0 + K^0b)/2, (p+pb)/2
          FF(I)=(FF(I)+FF(-I))/2.D0
          FF(-I)=FF(I)
        ENDDO
      ENDIF

C -- Gradient terms of parameters for fragmentation functions --
      DO J=1,NPAR
        DO I=1,NFF
          NI=9+NFF*(J-1)
          GRADFF(I,J)=(1.D0-T)*PDFJ1(NI+I)+T*PDFJ2(NI+I)
        ENDDO
      ENDDO

      IF(Q2.LT.CTHRE) THEN
        DO J=1,NPAR
          GRADFF(NFF-1,J)=0.D0
          IF(J.GT.NPAR-6) THEN
            DO I=1,NFF-2
              GRADFF(I,J)=0.D0
            ENDDO
          ENDIF
        ENDDO
      ENDIF 

      IF(Q2.LT.BTHRE) THEN
        DO J=1,NPAR
          GRADFF(NFF,J)=0.D0
          IF(J.GT.NPAR-3) THEN
            DO I=1,NFF-1
              GRADFF(I,J)=0.D0
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      DO J=1,NPAR
        DO I=1,3
          GRAD(I-1,J)=GRADFF(I,J)   ! g, u, d
          GRAD(-I,J)=GRADFF(I+3,J)  ! ub, [db, s]^pi+, [sb, s]^K+
          GRAD(I+2,J)=GRADFF(I+5,J) ! s, c, b
        ENDDO
        GRAD(-4,J)=GRAD(4,J) ! cb=c      
        GRAD(-5,J)=GRAD(5,J) ! bb=b     

        IF((ISET.EQ.3).OR.(ISET.EQ.4)) THEN ! K+
          GRAD(-3,J)=GRAD(-2,J)
          GRAD(-2,J)=GRAD(2,J)
        ENDIF
      ENDDO

      DO J=1,NPAR
        IF(ICHARGE.EQ.2) THEN
          IF(ISET.LT.3) THEN
            TEMP1=GRAD(-1,J);TEMP2=GRAD(2,J)
            GRAD(-1,J)=GRAD(1,J);GRAD(2,J)=GRAD(-2,J)
            GRAD(1,J)=TEMP1;GRAD(-2,J)=TEMP2

          ELSE IF((ISET.EQ.3).OR.(ISET.EQ.4)) THEN
            TEMP1=GRAD(-1,J);TEMP2=GRAD(3,J)
            GRAD(-1,J)=GRAD(1,J);GRAD(3,J)=GRAD(-3,J)
            GRAD(1,J)=TEMP1;GRAD(-3,J)=TEMP2

          ELSE IF(ISET.GT.4) THEN
            TEMP1=GRAD(-1,J);TEMP2=GRAD(1,J)
            GRAD(-1,J)=GRAD(-2,J);GRAD(1,J)=GRAD(2,J)
            GRAD(-2,J)=TEMP1;GRAD(2,J)=TEMP2
          ENDIF

        ELSE IF(ICHARGE.EQ.3) THEN
          IF((ISET.EQ.3).OR.(ISET.EQ.4)) THEN
            TEMP1=GRAD(-1,J);TEMP2=GRAD(1,J)
            GRAD(-1,J)=GRAD(-2,J);GRAD(1,J)=GRAD(2,J)
            GRAD(-2,J)=TEMP1;GRAD(2,J)=TEMP2
          ENDIF
          DO I=1,5
            GRAD(I,J)=(GRAD(I,J)+GRAD(-I,J))*0.5D0
            GRAD(-I,J)=GRAD(I,J)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE LSPLINE(N,X,Y,B,C,D,ISET,I,J)
C ---------------------------------------------------------------------
C CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
C INTERPOLATION SUBROUTINES ARE TAKEN FROM
C G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
C COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NQ=33, NX=53, ND=146)
      DIMENSION Y(NX,NQ,ND),B(NX,NQ,ND),C(NX,NQ,ND),D(NX,NQ,ND)
     1         ,X(NX) 
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1,J,I)=X(2)-X(1)
      C(2,J,I)=(Y(2,J,I)-Y(1,J,I))/D(1,J,I)
      DO 210 K=2,NM1
        D(K,J,I)=X(K+1)-X(K)
        B(K,J,I)=2.0D0*(D(K-1,J,I)+D(K,J,I))
        C(K+1,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
        C(K,J,I)=C(K+1,J,I)-C(K,J,I)
  210 CONTINUE
      B(1,J,I)=-D(1,J,I)
      B(N,J,I)=-D(N-1,J,I)
      C(1,J,I)=0.0D0
      C(N,J,I)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1,J,I)=C(3,J,I)/(X(4)-X(2))-C(2,J,I)/(X(3)-X(1))
      C(N,J,I)=C(N-1,J,I)/(X(N)-X(N-2))-C(N-2,J,I)/(X(N-1)-X(N-3))
      C(1,J,I)=C(1,J,I)*D(1,J,I)**2.0D0/(X(4)-X(1))
      C(N,J,I)=-C(N,J,I)*D(N-1,J,I)**2.0D0/(X(N)-X(N-3))
  215 CONTINUE
      DO 220 K=2,N
        T=D(K-1,J,I)/B(K-1,J,I)
        B(K,J,I)=B(K,J,I)-T*D(K-1,J,I)
        C(K,J,I)=C(K,J,I)-T*C(K-1,J,I)
  220 CONTINUE
      C(N,J,I)=C(N,J,I)/B(N,J,I)
      DO 230 IB=1,NM1
        K=N-IB
        C(K,J,I)=(C(K,J,I)-D(K,J,I)*C(K+1,J,I))/B(K,J,I)
  230 CONTINUE
      B(N,J,I)=(Y(N,J,I)-Y(NM1,J,I))/D(NM1,J,I)
     1        +D(NM1,J,I)*(C(NM1,J,I)+2.0D0*C(N,J,I))
      DO 240 K=1,NM1
        B(K,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
     1          -D(K,J,I)*(C(K+1,J,I)+2.0D0*C(K,J,I))
        D(K,J,I)=(C(K+1,J,I)-C(K,J,I))/D(K,J,I)
        C(K,J,I)=3.0D0*C(K,J,I)
  240 CONTINUE
      C(N,J,I)=3.0D0*C(N,J,I)
      D(N,J,I)=D(N-1,J,I)
      RETURN
  250 CONTINUE
      B(1,J,I)=(Y(2,J,I)-Y(1,J,I))/(X(2)-X(1))
      C(1,J,I)=0.0D0
      D(1,J,I)=0.0D0
      B(2,J,I)=B(1,J,I)
      C(2,J,I)=0.0D0
      D(2,J,I)=0.0D0
      RETURN
      END
C ---------------------------------------------------------------------
      INTEGER FUNCTION ISERCH(N,X,Y)
C ---------------------------------------------------------------------
C THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
C X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(53)

      MIN=1
      MAX=N+1

   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10

      ISERCH=MIN

      RETURN
      END
C *********************************************************************
C THE END OF THE PROGRAM.
C *********************************************************************
