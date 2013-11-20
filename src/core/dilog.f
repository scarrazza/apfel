c$$$****************************************************
c$$$      REAL*8 FUNCTION DDILOG(X)
c$$$      REAL*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
c$$$      REAL*8 C(0:18),H,ALFA,B0,B1,B2
c$$$      DATA ZERO /0.0D0/, ONE /1.0D0/
c$$$      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
c$$$      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
c$$$      DATA C( 0) / 0.42996 69356 08137 0D0/
c$$$      DATA C( 1) / 0.40975 98753 30771 1D0/
c$$$      DATA C( 2) /-0.01858 84366 50146 0D0/
c$$$      DATA C( 3) / 0.00145 75108 40622 7D0/
c$$$      DATA C( 4) /-0.00014 30418 44423 4D0/
c$$$      DATA C( 5) / 0.00001 58841 55418 8D0/
c$$$      DATA C( 6) /-0.00000 19078 49593 9D0/
c$$$      DATA C( 7) / 0.00000 02419 51808 5D0/
c$$$      DATA C( 8) /-0.00000 00319 33412 7D0/
c$$$      DATA C( 9) / 0.00000 00043 45450 6D0/
c$$$      DATA C(10) /-0.00000 00006 05784 8D0/
c$$$      DATA C(11) / 0.00000 00000 86121 0D0/
c$$$      DATA C(12) /-0.00000 00000 12443 3D0/
c$$$      DATA C(13) / 0.00000 00000 01822 6D0/
c$$$      DATA C(14) /-0.00000 00000 00270 1D0/
c$$$      DATA C(15) / 0.00000 00000 00040 4D0/
c$$$      DATA C(16) /-0.00000 00000 00006 1D0/
c$$$      DATA C(17) / 0.00000 00000 00000 9D0/
c$$$      DATA C(18) /-0.00000 00000 00000 1D0/
c$$$      IF(X .EQ. ONE) THEN
c$$$         DDILOG=PI6
c$$$         RETURN
c$$$      ELSE IF(X .EQ. MONE) THEN
c$$$         DDILOG=MALF*PI6
c$$$         RETURN
c$$$      END IF
c$$$      T=-X
c$$$      IF(T .LE. MTWO) THEN
c$$$         Y=MONE/(ONE+T)
c$$$         S=ONE
c$$$         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
c$$$      ELSE IF(T .LT. MONE) THEN
c$$$         Y=MONE-T
c$$$         S=MONE
c$$$         A=LOG(-T)
c$$$         A=-PI6+A*(A+LOG(ONE+ONE/T))
c$$$      ELSE IF(T .LE. MALF) THEN
c$$$         Y=(MONE-T)/T
c$$$         S=ONE
c$$$         A=LOG(-T)
c$$$         A=-PI6+A*(MALF*A+LOG(ONE+T))
c$$$      ELSE IF(T .LT. ZERO) THEN
c$$$         Y=-T/(ONE+T)
c$$$         S=MONE
c$$$         A=HALF*LOG(ONE+T)**2
c$$$      ELSE IF(T .LE. ONE) THEN
c$$$         Y=T
c$$$         S=ONE
c$$$         A=ZERO
c$$$      ELSE
c$$$         Y=ONE/T
c$$$         S=MONE
c$$$         A=PI6+HALF*LOG(T)**2
c$$$      END IF
c$$$      H=Y+Y-ONE
c$$$      ALFA=H+H
c$$$      B1=ZERO
c$$$      B2=ZERO
c$$$      DO 1 I = 18,0,-1
c$$$         B0=C(I)+ALFA*B1-B2
c$$$         B2=B1
c$$$    1 B1=B0
c$$$      DDILOG=-(S*(B0-H*B2)+A)
c$$$      RETURN
c$$$      END
c$$$

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION DDILOG(X)
      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/
      IF(X .EQ. ONE) THEN
         DDILOG=PI6
         RETURN
      ELSE IF(X .EQ. MONE) THEN
         DDILOG=MALF*PI6
         RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
         Y=MONE/(ONE+T)
         S=ONE
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
         Y=MONE-T
         S=MONE
         A=LOG(-T)
         A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
         Y=(MONE-T)/T
         S=ONE
         A=LOG(-T)
         A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
         Y=-T/(ONE+T)
         S=MONE
         A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
         Y=T
         S=ONE
         A=ZERO
      ELSE
         Y=ONE/T
         S=MONE
         A=PI6+HALF*LOG(T)**2
      END IF
      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 1 I = 18,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
    1 B1=B0
      DDILOG=-(S*(B0-H*B2)+A)
      RETURN
      END
