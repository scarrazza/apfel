*     -*-fortran-*-
*
*     Physical basis:
*      -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6 
*      gm  tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
*
*     Unified evolution basis:
*     0   1   2   3   4   5   6   7   8   9   10  11  12  13
*     g   gm  Sig Dsg Tu1 Tu2 Td1 Td2 V   DV  Vu1 Vu2 Vd1 Vd2
*
*     QCD Evolution basis:
*     0   1   2   3   4   5   6   7   8   9  10  11  12  13
*     gm  Sg   g   V  V3  V8 V15 V24 V35  T3  T8 T15 T24 T35
*
*     Tranformation from physical basis to Unified evolution basis
*
*    {{0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0}, 
*     {1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 
*     {0,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1}, 
*     {0,  1, -1,  1, -1,  1, -1,  0, -1,  1, -1,  1, -1,  1}, 
*     {0,  0,  0, -1,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0}, 
*     {0, -2,  0,  1,  0,  1,  0,  0,  0,  1,  0,  1,  0, -2}, 
*     {0,  0,  0,  0, -1,  0,  1,  0,  1,  0, -1,  0,  0,  0}, 
*     {0,  0, -2,  0,  1,  0,  1,  0,  1,  0,  1,  0, -2,  0}, 
*     {0, -1, -1, -1, -1, -1, -1,  0,  1,  1,  1,  1,  1,  1}, 
*     {0, -1,  1, -1,  1, -1,  1,  0, -1,  1, -1,  1, -1,  1}, 
*     {0,  0,  0,  1,  0, -1,  0,  0,  0,  1,  0, -1,  0,  0}, 
*     {0,  2,  0, -1,  0, -1,  0,  0,  0,  1,  0,  1,  0, -2}, 
*     {0,  0,  0,  0,  1,  0, -1,  0,  1,  0, -1,  0,  0,  0}, 
*     {0,  0,  2,  0, -1,  0, -1,  0,  1,  0,  1,  0, -2,  0}}
*
*     Tranformation from Unified evolution and physical basis (times 12)
*
*    {{0, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
*     {0,  0,  1,  1,  0, -2,  0,  0, -1, -1,  0,  2,  0,  0},
*     {0,  0,  1, -1,  0,  0,  0, -2, -1,  1,  0,  0,  0,  2},
*     {0,  0,  1,  1, -3,  1,  0,  0, -1, -1,  3, -1,  0,  0},
*     {0,  0,  1, -1,  0,  0, -3,  1, -1,  1,  0,  0,  3, -1},
*     {0,  0,  1,  1,  3,  1,  0,  0, -1, -1, -3, -1,  0,  0},
*     {0,  0,  1, -1,  0,  0,  3,  1, -1,  1,  0,  0, -3, -1},
*     {12, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
*     {0,  0,  1, -1,  0,  0,  3,  1,  1, -1,  0,  0,  3,  1},
*     {0,  0,  1,  1,  3,  1,  0,  0,  1,  1,  3,  1,  0,  0},
*     {0,  0,  1, -1,  0,  0, -3,  1,  1, -1,  0,  0, -3,  1},
*     {0,  0,  1,  1, -3,  1,  0,  0,  1,  1, -3,  1,  0,  0},
*     {0,  0,  1, -1,  0,  0,  0, -2,  1, -1,  0,  0,  0, -2},
*     {0,  0,  1,  1,  0, -2,  0,  0,  1,  1,  0, -2,  0,  0}}
*
*     Tranformation from QCD evolution basis to Unified evolution basis
*
*    {{0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,    1,  1/3, -1/3,  1/5, -1/5},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,  1/2,  1/6,  1/3,    0,    0},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,  1/2,  1/6, -1/6, 1/10,  2/5},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0, -1/2,  1/2,    0,    0,    0},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0, -1/2, -1/6,  1/6,  1/2,    0},
*     {0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0,    1,  1/3, -1/3,  1/5, -1/5,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0,  1/2,  1/6,  1/3,    0,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0,  1/2,  1/6, -1/6, 1/10,  2/5,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0, -1/2,  1/2,    0,    0,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0, -1/2, -1/6,  1/6,  1/2,    0,    0,    0,    0,    0,    0}}
*
*     Tranformation from unified evolution basis to QCD evolution basis
*
*    {{0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,  1/3,  1/2,  1/6, -1/2, -1/6},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,  1/3,  1/2,  1/6,  3/2, -1/6},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0, -2/3,    2, -1/3,    0,  1/3},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,  2/3,    0,  1/3,    0,  5/3},
*     {0,    0,    0,    0,    0,    0,    0,    0,    0,   -1,    0,    2,    0,    0},
*     {0,    0,    0,  1/3,  1/2,  1/6, -1/2, -1/6,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,  1/3,  1/2,  1/6,  3/2, -1/6,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0, -2/3,    2, -1/3,    0,  1/3,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,  2/3,    0,  1/3,    0,  5/3,    0,    0,    0,    0,    0,    0},
*     {0,    0,    0,   -1,    0,    2,    0,    0,    0,    0,    0,    0,    0,    0}}
*
      double precision Tev2phUni(0:13,0:13)
      double precision Tph2evUni(0:13,0:13)
*
      double precision TevQCD2evUni(0:13,0:13)
      double precision TevUni2evQCD(0:13,0:13)
*
      data Tph2evUni( 0, 0)   / 0d0 /
      data Tph2evUni( 0, 1)   / 0d0 /
      data Tph2evUni( 0, 2)   / 0d0 /
      data Tph2evUni( 0, 3)   / 0d0 /
      data Tph2evUni( 0, 4)   / 0d0 /
      data Tph2evUni( 0, 5)   / 0d0 /
      data Tph2evUni( 0, 6)   / 0d0 /
      data Tph2evUni( 0, 7)   / 1d0 /
      data Tph2evUni( 0, 8)   / 0d0 /
      data Tph2evUni( 0, 9)   / 0d0 /
      data Tph2evUni( 0,10)   / 0d0 /
      data Tph2evUni( 0,11)   / 0d0 /
      data Tph2evUni( 0,12)   / 0d0 /
      data Tph2evUni( 0,13)   / 0d0 /
      data Tph2evUni( 1, 0)   / 1d0 /
      data Tph2evUni( 1, 1)   / 0d0 /
      data Tph2evUni( 1, 2)   / 0d0 /
      data Tph2evUni( 1, 3)   / 0d0 /
      data Tph2evUni( 1, 4)   / 0d0 /
      data Tph2evUni( 1, 5)   / 0d0 /
      data Tph2evUni( 1, 6)   / 0d0 /
      data Tph2evUni( 1, 7)   / 0d0 /
      data Tph2evUni( 1, 8)   / 0d0 /
      data Tph2evUni( 1, 9)   / 0d0 /
      data Tph2evUni( 1,10)   / 0d0 /
      data Tph2evUni( 1,11)   / 0d0 /
      data Tph2evUni( 1,12)   / 0d0 /
      data Tph2evUni( 1,13)   / 0d0 /
      data Tph2evUni( 2, 0)   / 0d0 /
      data Tph2evUni( 2, 1)   / 1d0 /
      data Tph2evUni( 2, 2)   / 1d0 /
      data Tph2evUni( 2, 3)   / 1d0 /
      data Tph2evUni( 2, 4)   / 1d0 /
      data Tph2evUni( 2, 5)   / 1d0 /
      data Tph2evUni( 2, 6)   / 1d0 /
      data Tph2evUni( 2, 7)   / 0d0 /
      data Tph2evUni( 2, 8)   / 1d0 /
      data Tph2evUni( 2, 9)   / 1d0 /
      data Tph2evUni( 2,10)   / 1d0 /
      data Tph2evUni( 2,11)   / 1d0 /
      data Tph2evUni( 2,12)   / 1d0 /
      data Tph2evUni( 2,13)   / 1d0 /
      data Tph2evUni( 3, 0)   / 0d0 /
      data Tph2evUni( 3, 1)   / 1d0 /
      data Tph2evUni( 3, 2)   /-1d0 /
      data Tph2evUni( 3, 3)   / 1d0 /
      data Tph2evUni( 3, 4)   /-1d0 /
      data Tph2evUni( 3, 5)   / 1d0 /
      data Tph2evUni( 3, 6)   /-1d0 /
      data Tph2evUni( 3, 7)   / 0d0 /
      data Tph2evUni( 3, 8)   /-1d0 /
      data Tph2evUni( 3, 9)   / 1d0 /
      data Tph2evUni( 3,10)   /-1d0 /
      data Tph2evUni( 3,11)   / 1d0 /
      data Tph2evUni( 3,12)   /-1d0 /
      data Tph2evUni( 3,13)   / 1d0 /
      data Tph2evUni( 4, 0)   / 0d0 /
      data Tph2evUni( 4, 1)   / 0d0 /
      data Tph2evUni( 4, 2)   / 0d0 /
      data Tph2evUni( 4, 3)   /-1d0 /
      data Tph2evUni( 4, 4)   / 0d0 /
      data Tph2evUni( 4, 5)   / 1d0 /
      data Tph2evUni( 4, 6)   / 0d0 /
      data Tph2evUni( 4, 7)   / 0d0 /
      data Tph2evUni( 4, 8)   / 0d0 /
      data Tph2evUni( 4, 9)   / 1d0 /
      data Tph2evUni( 4,10)   / 0d0 /
      data Tph2evUni( 4,11)   /-1d0 /
      data Tph2evUni( 4,12)   / 0d0 /
      data Tph2evUni( 4,13)   / 0d0 /
      data Tph2evUni( 5, 0)   / 0d0 /
      data Tph2evUni( 5, 1)   /-2d0 /
      data Tph2evUni( 5, 2)   / 0d0 /
      data Tph2evUni( 5, 3)   / 1d0 /
      data Tph2evUni( 5, 4)   / 0d0 /
      data Tph2evUni( 5, 5)   / 1d0 /
      data Tph2evUni( 5, 6)   / 0d0 /
      data Tph2evUni( 5, 7)   / 0d0 /
      data Tph2evUni( 5, 8)   / 0d0 /
      data Tph2evUni( 5, 9)   / 1d0 /
      data Tph2evUni( 5,10)   / 0d0 /
      data Tph2evUni( 5,11)   / 1d0 /
      data Tph2evUni( 5,12)   / 0d0 /
      data Tph2evUni( 5,13)   /-2d0 /
      data Tph2evUni( 6, 0)   / 0d0 /
      data Tph2evUni( 6, 1)   / 0d0 /
      data Tph2evUni( 6, 2)   / 0d0 /
      data Tph2evUni( 6, 3)   / 0d0 /
      data Tph2evUni( 6, 4)   /-1d0 /
      data Tph2evUni( 6, 5)   / 0d0 /
      data Tph2evUni( 6, 6)   / 1d0 /
      data Tph2evUni( 6, 7)   / 0d0 /
      data Tph2evUni( 6, 8)   / 1d0 /
      data Tph2evUni( 6, 9)   / 0d0 /
      data Tph2evUni( 6,10)   /-1d0 /
      data Tph2evUni( 6,11)   / 0d0 /
      data Tph2evUni( 6,12)   / 0d0 /
      data Tph2evUni( 6,13)   / 0d0 /
      data Tph2evUni( 7, 0)   / 0d0 /
      data Tph2evUni( 7, 1)   / 0d0 /
      data Tph2evUni( 7, 2)   /-2d0 /
      data Tph2evUni( 7, 3)   / 0d0 /
      data Tph2evUni( 7, 4)   / 1d0 /
      data Tph2evUni( 7, 5)   / 0d0 /
      data Tph2evUni( 7, 6)   / 1d0 /
      data Tph2evUni( 7, 7)   / 0d0 /
      data Tph2evUni( 7, 8)   / 1d0 /
      data Tph2evUni( 7, 9)   / 0d0 /
      data Tph2evUni( 7,10)   / 1d0 /
      data Tph2evUni( 7,11)   / 0d0 /
      data Tph2evUni( 7,12)   /-2d0 /
      data Tph2evUni( 7,13)   / 0d0 /
      data Tph2evUni( 8, 0)   / 0d0 /
      data Tph2evUni( 8, 1)   /-1d0 /
      data Tph2evUni( 8, 2)   /-1d0 /
      data Tph2evUni( 8, 3)   /-1d0 /
      data Tph2evUni( 8, 4)   /-1d0 /
      data Tph2evUni( 8, 5)   /-1d0 /
      data Tph2evUni( 8, 6)   /-1d0 /
      data Tph2evUni( 8, 7)   / 0d0 /
      data Tph2evUni( 8, 8)   / 1d0 /
      data Tph2evUni( 8, 9)   / 1d0 /
      data Tph2evUni( 8,10)   / 1d0 /
      data Tph2evUni( 8,11)   / 1d0 /
      data Tph2evUni( 8,12)   / 1d0 /
      data Tph2evUni( 8,13)   / 1d0 /
      data Tph2evUni( 9, 0)   / 0d0 /
      data Tph2evUni( 9, 1)   /-1d0 /
      data Tph2evUni( 9, 2)   / 1d0 /
      data Tph2evUni( 9, 3)   /-1d0 /
      data Tph2evUni( 9, 4)   / 1d0 /
      data Tph2evUni( 9, 5)   /-1d0 /
      data Tph2evUni( 9, 6)   / 1d0 /
      data Tph2evUni( 9, 7)   / 0d0 /
      data Tph2evUni( 9, 8)   /-1d0 /
      data Tph2evUni( 9, 9)   / 1d0 /
      data Tph2evUni( 9,10)   /-1d0 /
      data Tph2evUni( 9,11)   / 1d0 /
      data Tph2evUni( 9,12)   /-1d0 /
      data Tph2evUni( 9,13)   / 1d0 /
      data Tph2evUni(10, 0)   / 0d0 /
      data Tph2evUni(10, 1)   / 0d0 /
      data Tph2evUni(10, 2)   / 0d0 /
      data Tph2evUni(10, 3)   / 1d0 /
      data Tph2evUni(10, 4)   / 0d0 /
      data Tph2evUni(10, 5)   /-1d0 /
      data Tph2evUni(10, 6)   / 0d0 /
      data Tph2evUni(10, 7)   / 0d0 /
      data Tph2evUni(10, 8)   / 0d0 /
      data Tph2evUni(10, 9)   / 1d0 /
      data Tph2evUni(10,10)   / 0d0 /
      data Tph2evUni(10,11)   /-1d0 /
      data Tph2evUni(10,12)   / 0d0 /
      data Tph2evUni(10,13)   / 0d0 /
      data Tph2evUni(11, 0)   / 0d0 /
      data Tph2evUni(11, 1)   / 2d0 /
      data Tph2evUni(11, 2)   / 0d0 /
      data Tph2evUni(11, 3)   /-1d0 /
      data Tph2evUni(11, 4)   / 0d0 /
      data Tph2evUni(11, 5)   /-1d0 /
      data Tph2evUni(11, 6)   / 0d0 /
      data Tph2evUni(11, 7)   / 0d0 /
      data Tph2evUni(11, 8)   / 0d0 /
      data Tph2evUni(11, 9)   / 1d0 /
      data Tph2evUni(11,10)   / 0d0 /
      data Tph2evUni(11,11)   / 1d0 /
      data Tph2evUni(11,12)   / 0d0 /
      data Tph2evUni(11,13)   /-2d0 /
      data Tph2evUni(12, 0)   / 0d0 /
      data Tph2evUni(12, 1)   / 0d0 /
      data Tph2evUni(12, 2)   / 0d0 /
      data Tph2evUni(12, 3)   / 0d0 /
      data Tph2evUni(12, 4)   / 1d0 /
      data Tph2evUni(12, 5)   / 0d0 /
      data Tph2evUni(12, 6)   /-1d0 /
      data Tph2evUni(12, 7)   / 0d0 /
      data Tph2evUni(12, 8)   / 1d0 /
      data Tph2evUni(12, 9)   / 0d0 /
      data Tph2evUni(12,10)   /-1d0 /
      data Tph2evUni(12,11)   / 0d0 /
      data Tph2evUni(12,12)   / 0d0 /
      data Tph2evUni(12,13)   / 0d0 /
      data Tph2evUni(13, 0)   / 0d0 /
      data Tph2evUni(13, 1)   / 0d0 /
      data Tph2evUni(13, 2)   / 2d0 /
      data Tph2evUni(13, 3)   / 0d0 /
      data Tph2evUni(13, 4)   /-1d0 /
      data Tph2evUni(13, 5)   / 0d0 /
      data Tph2evUni(13, 6)   /-1d0 /
      data Tph2evUni(13, 7)   / 0d0 /
      data Tph2evUni(13, 8)   / 1d0 /
      data Tph2evUni(13, 9)   / 0d0 /
      data Tph2evUni(13,10)   / 1d0 /
      data Tph2evUni(13,11)   / 0d0 /
      data Tph2evUni(13,12)   /-2d0 /
      data Tph2evUni(13,13)   / 0d0 /
*
      data Tev2phUni( 0, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 1)   /  1.00000000000000d0 /
      data Tev2phUni( 0, 2)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 3)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 5)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 7)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 8)   /  0.00000000000000d0 /
      data Tev2phUni( 0, 9)   /  0.00000000000000d0 /
      data Tev2phUni( 0,10)   /  0.00000000000000d0 /
      data Tev2phUni( 0,11)   /  0.00000000000000d0 /
      data Tev2phUni( 0,12)   /  0.00000000000000d0 /
      data Tev2phUni( 0,13)   /  0.00000000000000d0 /
      data Tev2phUni( 1, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 1, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 1, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 1, 3)   /  0.08333333333333d0 /
      data Tev2phUni( 1, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 1, 5)   / -0.16666666666667d0 /
      data Tev2phUni( 1, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 1, 7)   /  0.00000000000000d0 /
      data Tev2phUni( 1, 8)   / -0.08333333333333d0 /
      data Tev2phUni( 1, 9)   / -0.08333333333333d0 /
      data Tev2phUni( 1,10)   /  0.00000000000000d0 /
      data Tev2phUni( 1,11)   /  0.16666666666667d0 /
      data Tev2phUni( 1,12)   /  0.00000000000000d0 /
      data Tev2phUni( 1,13)   /  0.00000000000000d0 /
      data Tev2phUni( 2, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 2, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 2, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 2, 3)   / -0.08333333333333d0 /
      data Tev2phUni( 2, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 2, 5)   /  0.00000000000000d0 /
      data Tev2phUni( 2, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 2, 7)   / -0.16666666666667d0 /
      data Tev2phUni( 2, 8)   / -0.08333333333333d0 /
      data Tev2phUni( 2, 9)   /  0.08333333333333d0 /
      data Tev2phUni( 2,10)   /  0.00000000000000d0 /
      data Tev2phUni( 2,11)   /  0.00000000000000d0 /
      data Tev2phUni( 2,12)   /  0.00000000000000d0 /
      data Tev2phUni( 2,13)   /  0.16666666666667d0 /
      data Tev2phUni( 3, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 3, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 3, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 3, 3)   /  0.08333333333333d0 /
      data Tev2phUni( 3, 4)   / -0.25000000000000d0 /
      data Tev2phUni( 3, 5)   /  0.08333333333333d0 /
      data Tev2phUni( 3, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 3, 7)   /  0.00000000000000d0 /
      data Tev2phUni( 3, 8)   / -0.08333333333333d0 /
      data Tev2phUni( 3, 9)   / -0.08333333333333d0 /
      data Tev2phUni( 3,10)   /  0.25000000000000d0 /
      data Tev2phUni( 3,11)   / -0.08333333333333d0 /
      data Tev2phUni( 3,12)   /  0.00000000000000d0 /
      data Tev2phUni( 3,13)   /  0.00000000000000d0 /
      data Tev2phUni( 4, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 4, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 4, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 4, 3)   / -0.08333333333333d0 /
      data Tev2phUni( 4, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 4, 5)   /  0.00000000000000d0 /
      data Tev2phUni( 4, 6)   / -0.25000000000000d0 /
      data Tev2phUni( 4, 7)   /  0.08333333333333d0 /
      data Tev2phUni( 4, 8)   / -0.08333333333333d0 /
      data Tev2phUni( 4, 9)   /  0.08333333333333d0 /
      data Tev2phUni( 4,10)   /  0.00000000000000d0 /
      data Tev2phUni( 4,11)   /  0.00000000000000d0 /
      data Tev2phUni( 4,12)   /  0.25000000000000d0 /
      data Tev2phUni( 4,13)   / -0.08333333333333d0 /
      data Tev2phUni( 5, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 5, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 5, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 5, 3)   /  0.08333333333333d0 /
      data Tev2phUni( 5, 4)   /  0.25000000000000d0 /
      data Tev2phUni( 5, 5)   /  0.08333333333333d0 /
      data Tev2phUni( 5, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 5, 7)   /  0.00000000000000d0 /
      data Tev2phUni( 5, 8)   / -0.08333333333333d0 /
      data Tev2phUni( 5, 9)   / -0.08333333333333d0 /
      data Tev2phUni( 5,10)   / -0.25000000000000d0 /
      data Tev2phUni( 5,11)   / -0.08333333333333d0 /
      data Tev2phUni( 5,12)   /  0.00000000000000d0 /
      data Tev2phUni( 5,13)   /  0.00000000000000d0 /
      data Tev2phUni( 6, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 6, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 6, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 6, 3)   / -0.08333333333333d0 /
      data Tev2phUni( 6, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 6, 5)   /  0.00000000000000d0 /
      data Tev2phUni( 6, 6)   /  0.25000000000000d0 /
      data Tev2phUni( 6, 7)   /  0.08333333333333d0 /
      data Tev2phUni( 6, 8)   / -0.08333333333333d0 /
      data Tev2phUni( 6, 9)   /  0.08333333333333d0 /
      data Tev2phUni( 6,10)   /  0.00000000000000d0 /
      data Tev2phUni( 6,11)   /  0.00000000000000d0 /
      data Tev2phUni( 6,12)   / -0.25000000000000d0 /
      data Tev2phUni( 6,13)   / -0.08333333333333d0 /
      data Tev2phUni( 7, 0)   /  1.00000000000000d0 /
      data Tev2phUni( 7, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 2)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 3)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 5)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 7)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 8)   /  0.00000000000000d0 /
      data Tev2phUni( 7, 9)   /  0.00000000000000d0 /
      data Tev2phUni( 7,10)   /  0.00000000000000d0 /
      data Tev2phUni( 7,11)   /  0.00000000000000d0 /
      data Tev2phUni( 7,12)   /  0.00000000000000d0 /
      data Tev2phUni( 7,13)   /  0.00000000000000d0 /
      data Tev2phUni( 8, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 8, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 8, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 8, 3)   / -0.08333333333333d0 /
      data Tev2phUni( 8, 4)   /  0.00000000000000d0 /
      data Tev2phUni( 8, 5)   /  0.00000000000000d0 /
      data Tev2phUni( 8, 6)   /  0.25000000000000d0 /
      data Tev2phUni( 8, 7)   /  0.08333333333333d0 /
      data Tev2phUni( 8, 8)   /  0.08333333333333d0 /
      data Tev2phUni( 8, 9)   / -0.08333333333333d0 /
      data Tev2phUni( 8,10)   /  0.00000000000000d0 /
      data Tev2phUni( 8,11)   /  0.00000000000000d0 /
      data Tev2phUni( 8,12)   /  0.25000000000000d0 /
      data Tev2phUni( 8,13)   /  0.08333333333333d0 /
      data Tev2phUni( 9, 0)   /  0.00000000000000d0 /
      data Tev2phUni( 9, 1)   /  0.00000000000000d0 /
      data Tev2phUni( 9, 2)   /  0.08333333333333d0 /
      data Tev2phUni( 9, 3)   /  0.08333333333333d0 /
      data Tev2phUni( 9, 4)   /  0.25000000000000d0 /
      data Tev2phUni( 9, 5)   /  0.08333333333333d0 /
      data Tev2phUni( 9, 6)   /  0.00000000000000d0 /
      data Tev2phUni( 9, 7)   /  0.00000000000000d0 /
      data Tev2phUni( 9, 8)   /  0.08333333333333d0 /
      data Tev2phUni( 9, 9)   /  0.08333333333333d0 /
      data Tev2phUni( 9,10)   /  0.25000000000000d0 /
      data Tev2phUni( 9,11)   /  0.08333333333333d0 /
      data Tev2phUni( 9,12)   /  0.00000000000000d0 /
      data Tev2phUni( 9,13)   /  0.00000000000000d0 /
      data Tev2phUni(10, 0)   /  0.00000000000000d0 /
      data Tev2phUni(10, 1)   /  0.00000000000000d0 /
      data Tev2phUni(10, 2)   /  0.08333333333333d0 /
      data Tev2phUni(10, 3)   / -0.08333333333333d0 /
      data Tev2phUni(10, 4)   /  0.00000000000000d0 /
      data Tev2phUni(10, 5)   /  0.00000000000000d0 /
      data Tev2phUni(10, 6)   / -0.25000000000000d0 /
      data Tev2phUni(10, 7)   /  0.08333333333333d0 /
      data Tev2phUni(10, 8)   /  0.08333333333333d0 /
      data Tev2phUni(10, 9)   / -0.08333333333333d0 /
      data Tev2phUni(10,10)   /  0.00000000000000d0 /
      data Tev2phUni(10,11)   /  0.00000000000000d0 /
      data Tev2phUni(10,12)   / -0.25000000000000d0 /
      data Tev2phUni(10,13)   /  0.08333333333333d0 /
      data Tev2phUni(11, 0)   /  0.00000000000000d0 /
      data Tev2phUni(11, 1)   /  0.00000000000000d0 /
      data Tev2phUni(11, 2)   /  0.08333333333333d0 /
      data Tev2phUni(11, 3)   /  0.08333333333333d0 /
      data Tev2phUni(11, 4)   / -0.25000000000000d0 /
      data Tev2phUni(11, 5)   /  0.08333333333333d0 /
      data Tev2phUni(11, 6)   /  0.00000000000000d0 /
      data Tev2phUni(11, 7)   /  0.00000000000000d0 /
      data Tev2phUni(11, 8)   /  0.08333333333333d0 /
      data Tev2phUni(11, 9)   /  0.08333333333333d0 /
      data Tev2phUni(11,10)   / -0.25000000000000d0 /
      data Tev2phUni(11,11)   /  0.08333333333333d0 /
      data Tev2phUni(11,12)   /  0.00000000000000d0 /
      data Tev2phUni(11,13)   /  0.00000000000000d0 /
      data Tev2phUni(12, 0)   /  0.00000000000000d0 /
      data Tev2phUni(12, 1)   /  0.00000000000000d0 /
      data Tev2phUni(12, 2)   /  0.08333333333333d0 /
      data Tev2phUni(12, 3)   / -0.08333333333333d0 /
      data Tev2phUni(12, 4)   /  0.00000000000000d0 /
      data Tev2phUni(12, 5)   /  0.00000000000000d0 /
      data Tev2phUni(12, 6)   /  0.00000000000000d0 /
      data Tev2phUni(12, 7)   / -0.16666666666667d0 /
      data Tev2phUni(12, 8)   /  0.08333333333333d0 /
      data Tev2phUni(12, 9)   / -0.08333333333333d0 /
      data Tev2phUni(12,10)   /  0.00000000000000d0 /
      data Tev2phUni(12,11)   /  0.00000000000000d0 /
      data Tev2phUni(12,12)   /  0.00000000000000d0 /
      data Tev2phUni(12,13)   / -0.16666666666667d0 /
      data Tev2phUni(13, 0)   /  0.00000000000000d0 /
      data Tev2phUni(13, 1)   /  0.00000000000000d0 /
      data Tev2phUni(13, 2)   /  0.08333333333333d0 /
      data Tev2phUni(13, 3)   /  0.08333333333333d0 /
      data Tev2phUni(13, 4)   /  0.00000000000000d0 /
      data Tev2phUni(13, 5)   / -0.16666666666667d0 /
      data Tev2phUni(13, 6)   /  0.00000000000000d0 /
      data Tev2phUni(13, 7)   /  0.00000000000000d0 /
      data Tev2phUni(13, 8)   /  0.08333333333333d0 /
      data Tev2phUni(13, 9)   /  0.08333333333333d0 /
      data Tev2phUni(13,10)   /  0.00000000000000d0 /
      data Tev2phUni(13,11)   / -0.16666666666667d0 /
      data Tev2phUni(13,12)   /  0.00000000000000d0 /
      data Tev2phUni(13,13)   /  0.00000000000000d0 /
*
      data TevQCD2evUni( 0, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 2)   /  1.00000000000000d0 /
      data TevQCD2evUni( 0, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0,10)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0,11)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 0,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 0)   /  1.00000000000000d0 /
      data TevQCD2evUni( 1, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1,10)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1,11)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 1,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 1)   /  1.00000000000000d0 /
      data TevQCD2evUni( 2, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2,10)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2,11)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 2,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 3, 9)   /  1.00000000000000d0 /
      data TevQCD2evUni( 3,10)   /  0.33333333333333d0 /
      data TevQCD2evUni( 3,11)   / -0.33333333333333d0 /
      data TevQCD2evUni( 3,12)   /  0.20000000000000d0 /
      data TevQCD2evUni( 3,13)   / -0.20000000000000d0 /
      data TevQCD2evUni( 4, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4, 9)   /  0.50000000000000d0 /
      data TevQCD2evUni( 4,10)   /  0.16666666666667d0 /
      data TevQCD2evUni( 4,11)   /  0.33333333333333d0 /
      data TevQCD2evUni( 4,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 4,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 5, 9)   /  0.50000000000000d0 /
      data TevQCD2evUni( 5,10)   /  0.16666666666667d0 /
      data TevQCD2evUni( 5,11)   / -0.16666666666667d0 /
      data TevQCD2evUni( 5,12)   /  0.10000000000000d0 /
      data TevQCD2evUni( 5,13)   /  0.40000000000000d0 /
      data TevQCD2evUni( 6, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6, 9)   / -0.50000000000000d0 /
      data TevQCD2evUni( 6,10)   /  0.50000000000000d0 /
      data TevQCD2evUni( 6,11)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 6,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 7, 9)   / -0.50000000000000d0 /
      data TevQCD2evUni( 7,10)   / -0.16666666666667d0 /
      data TevQCD2evUni( 7,11)   /  0.16666666666667d0 /
      data TevQCD2evUni( 7,12)   /  0.50000000000000d0 /
      data TevQCD2evUni( 7,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 3)   /  1.00000000000000d0 /
      data TevQCD2evUni( 8, 4)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 5)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8,10)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8,11)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 8,13)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9, 4)   /  1.00000000000000d0 /
      data TevQCD2evUni( 9, 5)   /  0.33333333333333d0 /
      data TevQCD2evUni( 9, 6)   / -0.33333333333333d0 /
      data TevQCD2evUni( 9, 7)   /  0.20000000000000d0 /
      data TevQCD2evUni( 9, 8)   / -0.20000000000000d0 /
      data TevQCD2evUni( 9, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9,10)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9,11)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9,12)   /  0.00000000000000d0 /
      data TevQCD2evUni( 9,13)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 4)   /  0.50000000000000d0 /
      data TevQCD2evUni(10, 5)   /  0.16666666666667d0 /
      data TevQCD2evUni(10, 6)   /  0.33333333333333d0 /
      data TevQCD2evUni(10, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni(10, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni(10,10)   /  0.00000000000000d0 /
      data TevQCD2evUni(10,11)   /  0.00000000000000d0 /
      data TevQCD2evUni(10,12)   /  0.00000000000000d0 /
      data TevQCD2evUni(10,13)   /  0.00000000000000d0 /
      data TevQCD2evUni(11, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni(11, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni(11, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni(11, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni(11, 4)   /  0.50000000000000d0 /
      data TevQCD2evUni(11, 5)   /  0.16666666666667d0 /
      data TevQCD2evUni(11, 6)   / -0.16666666666667d0 /
      data TevQCD2evUni(11, 7)   /  0.10000000000000d0 /
      data TevQCD2evUni(11, 8)   /  0.40000000000000d0 /
      data TevQCD2evUni(11, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni(11,10)   /  0.00000000000000d0 /
      data TevQCD2evUni(11,11)   /  0.00000000000000d0 /
      data TevQCD2evUni(11,12)   /  0.00000000000000d0 /
      data TevQCD2evUni(11,13)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 4)   / -0.50000000000000d0 /
      data TevQCD2evUni(12, 5)   /  0.50000000000000d0 /
      data TevQCD2evUni(12, 6)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 7)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni(12, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni(12,10)   /  0.00000000000000d0 /
      data TevQCD2evUni(12,11)   /  0.00000000000000d0 /
      data TevQCD2evUni(12,12)   /  0.00000000000000d0 /
      data TevQCD2evUni(12,13)   /  0.00000000000000d0 /
      data TevQCD2evUni(13, 0)   /  0.00000000000000d0 /
      data TevQCD2evUni(13, 1)   /  0.00000000000000d0 /
      data TevQCD2evUni(13, 2)   /  0.00000000000000d0 /
      data TevQCD2evUni(13, 3)   /  0.00000000000000d0 /
      data TevQCD2evUni(13, 4)   / -0.50000000000000d0 /
      data TevQCD2evUni(13, 5)   / -0.16666666666667d0 /
      data TevQCD2evUni(13, 6)   /  0.16666666666667d0 /
      data TevQCD2evUni(13, 7)   /  0.50000000000000d0 /
      data TevQCD2evUni(13, 8)   /  0.00000000000000d0 /
      data TevQCD2evUni(13, 9)   /  0.00000000000000d0 /
      data TevQCD2evUni(13,10)   /  0.00000000000000d0 /
      data TevQCD2evUni(13,11)   /  0.00000000000000d0 /
      data TevQCD2evUni(13,12)   /  0.00000000000000d0 /
      data TevQCD2evUni(13,13)   /  0.00000000000000d0 /
*
      data TevUni2evQCD( 0, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 1)   /  1.00000000000000d0/
      data TevUni2evQCD( 0, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 0, 9)   /  0.00000000000000d0/
      data TevUni2evQCD( 0,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 0,11)   /  0.00000000000000d0/
      data TevUni2evQCD( 0,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 0,13)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 2)   /  1.00000000000000d0/
      data TevUni2evQCD( 1, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 1, 9)   /  0.00000000000000d0/
      data TevUni2evQCD( 1,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 1,11)   /  0.00000000000000d0/
      data TevUni2evQCD( 1,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 1,13)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 0)   /  1.00000000000000d0/
      data TevUni2evQCD( 2, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 2, 9)   /  0.00000000000000d0/
      data TevUni2evQCD( 2,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 2,11)   /  0.00000000000000d0/
      data TevUni2evQCD( 2,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 2,13)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 3, 8)   /  1.00000000000000d0/
      data TevUni2evQCD( 3, 9)   /  0.00000000000000d0/
      data TevUni2evQCD( 3,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 3,11)   /  0.00000000000000d0/
      data TevUni2evQCD( 3,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 3,13)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 4, 9)   /  0.33333333333333d0/
      data TevUni2evQCD( 4,10)   /  0.50000000000000d0/
      data TevUni2evQCD( 4,11)   /  0.16666666666667d0/
      data TevUni2evQCD( 4,12)   / -0.50000000000000d0/
      data TevUni2evQCD( 4,13)   / -0.16666666666667d0/
      data TevUni2evQCD( 5, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 5, 9)   /  0.33333333333333d0/
      data TevUni2evQCD( 5,10)   /  0.50000000000000d0/
      data TevUni2evQCD( 5,11)   /  0.16666666666667d0/
      data TevUni2evQCD( 5,12)   /  1.50000000000000d0/
      data TevUni2evQCD( 5,13)   / -0.16666666666667d0/
      data TevUni2evQCD( 6, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 6, 9)   / -0.66666666666667d0/
      data TevUni2evQCD( 6,10)   /  2.00000000000000d0/
      data TevUni2evQCD( 6,11)   / -0.33333333333333d0/
      data TevUni2evQCD( 6,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 6,13)   /  0.33333333333333d0/
      data TevUni2evQCD( 7, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 7, 9)   /  0.66666666666667d0/
      data TevUni2evQCD( 7,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 7,11)   /  0.33333333333333d0/
      data TevUni2evQCD( 7,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 7,13)   /  1.66666666666667d0/
      data TevUni2evQCD( 8, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 3)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 4)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 5)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 6)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 7)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 8, 9)   / -1.00000000000000d0/
      data TevUni2evQCD( 8,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 8,11)   /  2.00000000000000d0/
      data TevUni2evQCD( 8,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 8,13)   /  0.00000000000000d0/
      data TevUni2evQCD( 9, 0)   /  0.00000000000000d0/
      data TevUni2evQCD( 9, 1)   /  0.00000000000000d0/
      data TevUni2evQCD( 9, 2)   /  0.00000000000000d0/
      data TevUni2evQCD( 9, 3)   /  0.33333333333333d0/
      data TevUni2evQCD( 9, 4)   /  0.50000000000000d0/
      data TevUni2evQCD( 9, 5)   /  0.16666666666667d0/
      data TevUni2evQCD( 9, 6)   / -0.50000000000000d0/
      data TevUni2evQCD( 9, 7)   / -0.16666666666667d0/
      data TevUni2evQCD( 9, 8)   /  0.00000000000000d0/
      data TevUni2evQCD( 9, 9)   /  0.00000000000000d0/
      data TevUni2evQCD( 9,10)   /  0.00000000000000d0/
      data TevUni2evQCD( 9,11)   /  0.00000000000000d0/
      data TevUni2evQCD( 9,12)   /  0.00000000000000d0/
      data TevUni2evQCD( 9,13)   /  0.00000000000000d0/
      data TevUni2evQCD(10, 0)   /  0.00000000000000d0/
      data TevUni2evQCD(10, 1)   /  0.00000000000000d0/
      data TevUni2evQCD(10, 2)   /  0.00000000000000d0/
      data TevUni2evQCD(10, 3)   /  0.33333333333333d0/
      data TevUni2evQCD(10, 4)   /  0.50000000000000d0/
      data TevUni2evQCD(10, 5)   /  0.16666666666667d0/
      data TevUni2evQCD(10, 6)   /  1.50000000000000d0/
      data TevUni2evQCD(10, 7)   / -0.16666666666667d0/
      data TevUni2evQCD(10, 8)   /  0.00000000000000d0/
      data TevUni2evQCD(10, 9)   /  0.00000000000000d0/
      data TevUni2evQCD(10,10)   /  0.00000000000000d0/
      data TevUni2evQCD(10,11)   /  0.00000000000000d0/
      data TevUni2evQCD(10,12)   /  0.00000000000000d0/
      data TevUni2evQCD(10,13)   /  0.00000000000000d0/
      data TevUni2evQCD(11, 0)   /  0.00000000000000d0/
      data TevUni2evQCD(11, 1)   /  0.00000000000000d0/
      data TevUni2evQCD(11, 2)   /  0.00000000000000d0/
      data TevUni2evQCD(11, 3)   / -0.66666666666667d0/
      data TevUni2evQCD(11, 4)   /  2.00000000000000d0/
      data TevUni2evQCD(11, 5)   / -0.33333333333333d0/
      data TevUni2evQCD(11, 6)   /  0.00000000000000d0/
      data TevUni2evQCD(11, 7)   /  0.33333333333333d0/
      data TevUni2evQCD(11, 8)   /  0.00000000000000d0/
      data TevUni2evQCD(11, 9)   /  0.00000000000000d0/
      data TevUni2evQCD(11,10)   /  0.00000000000000d0/
      data TevUni2evQCD(11,11)   /  0.00000000000000d0/
      data TevUni2evQCD(11,12)   /  0.00000000000000d0/
      data TevUni2evQCD(11,13)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 0)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 1)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 2)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 3)   /  0.66666666666667d0/
      data TevUni2evQCD(12, 4)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 5)   /  0.33333333333333d0/
      data TevUni2evQCD(12, 6)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 7)   /  1.66666666666667d0/
      data TevUni2evQCD(12, 8)   /  0.00000000000000d0/
      data TevUni2evQCD(12, 9)   /  0.00000000000000d0/
      data TevUni2evQCD(12,10)   /  0.00000000000000d0/
      data TevUni2evQCD(12,11)   /  0.00000000000000d0/
      data TevUni2evQCD(12,12)   /  0.00000000000000d0/
      data TevUni2evQCD(12,13)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 0)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 1)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 2)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 3)   / -1.00000000000000d0/
      data TevUni2evQCD(13, 4)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 5)   /  2.00000000000000d0/
      data TevUni2evQCD(13, 6)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 7)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 8)   /  0.00000000000000d0/
      data TevUni2evQCD(13, 9)   /  0.00000000000000d0/
      data TevUni2evQCD(13,10)   /  0.00000000000000d0/
      data TevUni2evQCD(13,11)   /  0.00000000000000d0/
      data TevUni2evQCD(13,12)   /  0.00000000000000d0/
      data TevUni2evQCD(13,13)   /  0.00000000000000d0/
