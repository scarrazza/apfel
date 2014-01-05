*     -*-fortran-*-
*
*     Tranformations between QCD evolution and physical basis
*
      double precision Tev2phQCD(0:13,0:13)
      double precision Tph2evQCD(0:13,0:13)
*
      data Tev2phQCD(0,0)   / 1.000000000000000d0 /
      data Tev2phQCD(0,1)   / 0.000000000000000d0 /
      data Tev2phQCD(0,2)   / 0.000000000000000d0 /
      data Tev2phQCD(0,3)   / 0.000000000000000d0 /
      data Tev2phQCD(0,4)   / 0.000000000000000d0 /
      data Tev2phQCD(0,5)   / 0.000000000000000d0 /
      data Tev2phQCD(0,6)   / 0.000000000000000d0 /
      data Tev2phQCD(0,7)   / 0.000000000000000d0 /
      data Tev2phQCD(0,8)   / 0.000000000000000d0 /
      data Tev2phQCD(0,9)   / 0.000000000000000d0 /
      data Tev2phQCD(0,10)  / 0.000000000000000d0 /
      data Tev2phQCD(0,11)  / 0.000000000000000d0 /
      data Tev2phQCD(0,12)  / 0.000000000000000d0 /
      data Tev2phQCD(0,13)  / 0.000000000000000d0 /
      data Tev2phQCD(1,0)   / 0.000000000000000d0 /
      data Tev2phQCD(1,1)   / 0.083333333333333d0 /
      data Tev2phQCD(1,2)   / 0.000000000000000d0 /
      data Tev2phQCD(1,3)   /-0.083333333333333d0 /
      data Tev2phQCD(1,4)   / 0.000000000000000d0 /
      data Tev2phQCD(1,5)   / 0.000000000000000d0 /
      data Tev2phQCD(1,6)   / 0.000000000000000d0 /
      data Tev2phQCD(1,7)   / 0.000000000000000d0 /
      data Tev2phQCD(1,8)   / 0.083333333333333d0 /
      data Tev2phQCD(1,9)   / 0.000000000000000d0 /
      data Tev2phQCD(1,10)  / 0.000000000000000d0 /
      data Tev2phQCD(1,11)  / 0.000000000000000d0 /
      data Tev2phQCD(1,12)  / 0.000000000000000d0 /
      data Tev2phQCD(1,13)  /-0.083333333333333d0 /
      data Tev2phQCD(2,0)   / 0.000000000000000d0 /
      data Tev2phQCD(2,1)   / 0.083333333333333d0 /
      data Tev2phQCD(2,2)   / 0.000000000000000d0 /
      data Tev2phQCD(2,3)   /-0.083333333333333d0 /
      data Tev2phQCD(2,4)   / 0.000000000000000d0 /
      data Tev2phQCD(2,5)   / 0.000000000000000d0 /
      data Tev2phQCD(2,6)   / 0.000000000000000d0 /
      data Tev2phQCD(2,7)   / 0.100000000000000d0 /
      data Tev2phQCD(2,8)   /-0.016666666666667d0 /
      data Tev2phQCD(2,9)   / 0.000000000000000d0 /
      data Tev2phQCD(2,10)  / 0.000000000000000d0 /
      data Tev2phQCD(2,11)  / 0.000000000000000d0 /
      data Tev2phQCD(2,12)  /-0.100000000000000d0 /
      data Tev2phQCD(2,13)  / 0.016666666666667d0 /
      data Tev2phQCD(3,0)   / 0.000000000000000d0 /
      data Tev2phQCD(3,1)   / 0.083333333333333d0 /
      data Tev2phQCD(3,2)   / 0.000000000000000d0 /
      data Tev2phQCD(3,3)   /-0.083333333333333d0 /
      data Tev2phQCD(3,4)   / 0.000000000000000d0 /
      data Tev2phQCD(3,5)   / 0.000000000000000d0 /
      data Tev2phQCD(3,6)   / 0.125000000000000d0 /
      data Tev2phQCD(3,7)   /-0.025000000000000d0 /
      data Tev2phQCD(3,8)   /-0.016666666666667d0 /
      data Tev2phQCD(3,9)   / 0.000000000000000d0 /
      data Tev2phQCD(3,10)  / 0.000000000000000d0 /
      data Tev2phQCD(3,11)  /-0.125000000000000d0 /
      data Tev2phQCD(3,12)  / 0.025000000000000d0 /
      data Tev2phQCD(3,13)  / 0.016666666666667d0 /
      data Tev2phQCD(4,0)   / 0.000000000000000d0 /
      data Tev2phQCD(4,1)   / 0.083333333333333d0 /
      data Tev2phQCD(4,2)   / 0.000000000000000d0 /
      data Tev2phQCD(4,3)   /-0.083333333333333d0 /
      data Tev2phQCD(4,4)   / 0.000000000000000d0 /
      data Tev2phQCD(4,5)   / 0.166666666666667d0 /
      data Tev2phQCD(4,6)   /-0.041666666666667d0 /
      data Tev2phQCD(4,7)   /-0.025000000000000d0 /
      data Tev2phQCD(4,8)   /-0.016666666666667d0 /
      data Tev2phQCD(4,9)   / 0.000000000000000d0 /
      data Tev2phQCD(4,10)  /-0.166666666666667d0 /
      data Tev2phQCD(4,11)  / 0.041666666666667d0 /
      data Tev2phQCD(4,12)  / 0.025000000000000d0 /
      data Tev2phQCD(4,13)  / 0.016666666666667d0 /
      data Tev2phQCD(5,0)   / 0.000000000000000d0 /
      data Tev2phQCD(5,1)   / 0.083333333333333d0 /
      data Tev2phQCD(5,2)   / 0.000000000000000d0 /
      data Tev2phQCD(5,3)   /-0.083333333333333d0 /
      data Tev2phQCD(5,4)   /-0.250000000000000d0 /
      data Tev2phQCD(5,5)   /-0.083333333333333d0 /
      data Tev2phQCD(5,6)   /-0.041666666666667d0 /
      data Tev2phQCD(5,7)   /-0.025000000000000d0 /
      data Tev2phQCD(5,8)   /-0.016666666666667d0 /
      data Tev2phQCD(5,9)   / 0.250000000000000d0 /
      data Tev2phQCD(5,10)  / 0.083333333333333d0 /
      data Tev2phQCD(5,11)  / 0.041666666666667d0 /
      data Tev2phQCD(5,12)  / 0.025000000000000d0 /
      data Tev2phQCD(5,13)  / 0.016666666666667d0 /
      data Tev2phQCD(6,0)   / 0.000000000000000d0 /
      data Tev2phQCD(6,1)   / 0.083333333333333d0 /
      data Tev2phQCD(6,2)   / 0.000000000000000d0 /
      data Tev2phQCD(6,3)   /-0.083333333333333d0 /
      data Tev2phQCD(6,4)   / 0.250000000000000d0 /
      data Tev2phQCD(6,5)   /-0.083333333333333d0 /
      data Tev2phQCD(6,6)   /-0.041666666666667d0 /
      data Tev2phQCD(6,7)   /-0.025000000000000d0 /
      data Tev2phQCD(6,8)   /-0.016666666666667d0 /
      data Tev2phQCD(6,9)   /-0.250000000000000d0 /
      data Tev2phQCD(6,10)  / 0.083333333333333d0 /
      data Tev2phQCD(6,11)  / 0.041666666666667d0 /
      data Tev2phQCD(6,12)  / 0.025000000000000d0 /
      data Tev2phQCD(6,13)  / 0.016666666666667d0 /
      data Tev2phQCD(7,0)   / 0.000000000000000d0 /
      data Tev2phQCD(7,1)   / 0.000000000000000d0 /
      data Tev2phQCD(7,2)   / 1.000000000000000d0 /
      data Tev2phQCD(7,3)   / 0.000000000000000d0 /
      data Tev2phQCD(7,4)   / 0.000000000000000d0 /
      data Tev2phQCD(7,5)   / 0.000000000000000d0 /
      data Tev2phQCD(7,6)   / 0.000000000000000d0 /
      data Tev2phQCD(7,7)   / 0.000000000000000d0 /
      data Tev2phQCD(7,8)   / 0.000000000000000d0 /
      data Tev2phQCD(7,9)   / 0.000000000000000d0 /
      data Tev2phQCD(7,10)  / 0.000000000000000d0 /
      data Tev2phQCD(7,11)  / 0.000000000000000d0 /
      data Tev2phQCD(7,12)  / 0.000000000000000d0 /
      data Tev2phQCD(7,13)  / 0.000000000000000d0 /
      data Tev2phQCD(8,0)   / 0.000000000000000d0 /
      data Tev2phQCD(8,1)   / 0.083333333333333d0 /
      data Tev2phQCD(8,2)   / 0.000000000000000d0 /
      data Tev2phQCD(8,3)   / 0.083333333333333d0 /
      data Tev2phQCD(8,4)   /-0.250000000000000d0 /
      data Tev2phQCD(8,5)   / 0.083333333333333d0 /
      data Tev2phQCD(8,6)   / 0.041666666666667d0 /
      data Tev2phQCD(8,7)   / 0.025000000000000d0 /
      data Tev2phQCD(8,8)   / 0.016666666666667d0 /
      data Tev2phQCD(8,9)   /-0.250000000000000d0 /
      data Tev2phQCD(8,10)  / 0.083333333333333d0 /
      data Tev2phQCD(8,11)  / 0.041666666666667d0 /
      data Tev2phQCD(8,12)  / 0.025000000000000d0 /
      data Tev2phQCD(8,13)  / 0.016666666666667d0 /
      data Tev2phQCD(9,0)   / 0.000000000000000d0 /
      data Tev2phQCD(9,1)   / 0.083333333333333d0 /
      data Tev2phQCD(9,2)   / 0.000000000000000d0 /
      data Tev2phQCD(9,3)   / 0.083333333333333d0 /
      data Tev2phQCD(9,4)   / 0.250000000000000d0 /
      data Tev2phQCD(9,5)   / 0.083333333333333d0 /
      data Tev2phQCD(9,6)   / 0.041666666666667d0 /
      data Tev2phQCD(9,7)   / 0.025000000000000d0 /
      data Tev2phQCD(9,8)   / 0.016666666666667d0 /
      data Tev2phQCD(9,9)   / 0.250000000000000d0 /
      data Tev2phQCD(9,10)  / 0.083333333333333d0 /
      data Tev2phQCD(9,11)  / 0.041666666666667d0 /
      data Tev2phQCD(9,12)  / 0.025000000000000d0 /
      data Tev2phQCD(9,13)  / 0.016666666666667d0 /
      data Tev2phQCD(10,0)  / 0.000000000000000d0 /
      data Tev2phQCD(10,1)  / 0.083333333333333d0 /
      data Tev2phQCD(10,2)  / 0.000000000000000d0 /
      data Tev2phQCD(10,3)  / 0.083333333333333d0 /
      data Tev2phQCD(10,4)  / 0.000000000000000d0 /
      data Tev2phQCD(10,5)  /-0.166666666666667d0 /
      data Tev2phQCD(10,6)  / 0.041666666666667d0 /
      data Tev2phQCD(10,7)  / 0.025000000000000d0 /
      data Tev2phQCD(10,8)  / 0.016666666666667d0 /
      data Tev2phQCD(10,9)  / 0.000000000000000d0 /
      data Tev2phQCD(10,10) /-0.166666666666667d0 /
      data Tev2phQCD(10,11) / 0.041666666666667d0 /
      data Tev2phQCD(10,12) / 0.025000000000000d0 /
      data Tev2phQCD(10,13) / 0.016666666666667d0 /
      data Tev2phQCD(11,0)  / 0.000000000000000d0 /
      data Tev2phQCD(11,1)  / 0.083333333333333d0 /
      data Tev2phQCD(11,2)  / 0.000000000000000d0 /
      data Tev2phQCD(11,3)  / 0.083333333333333d0 /
      data Tev2phQCD(11,4)  / 0.000000000000000d0 /
      data Tev2phQCD(11,5)  / 0.000000000000000d0 /
      data Tev2phQCD(11,6)  /-0.125000000000000d0 /
      data Tev2phQCD(11,7)  / 0.025000000000000d0 /
      data Tev2phQCD(11,8)  / 0.016666666666667d0 /
      data Tev2phQCD(11,9)  / 0.000000000000000d0 /
      data Tev2phQCD(11,10) / 0.000000000000000d0 /
      data Tev2phQCD(11,11) /-0.125000000000000d0 /
      data Tev2phQCD(11,12) / 0.025000000000000d0 /
      data Tev2phQCD(11,13) / 0.016666666666667d0 /
      data Tev2phQCD(12,0)  / 0.000000000000000d0 /
      data Tev2phQCD(12,1)  / 0.083333333333333d0 /
      data Tev2phQCD(12,2)  / 0.000000000000000d0 /
      data Tev2phQCD(12,3)  / 0.083333333333333d0 /
      data Tev2phQCD(12,4)  / 0.000000000000000d0 /
      data Tev2phQCD(12,5)  / 0.000000000000000d0 /
      data Tev2phQCD(12,6)  / 0.000000000000000d0 /
      data Tev2phQCD(12,7)  /-0.100000000000000d0 /
      data Tev2phQCD(12,8)  / 0.016666666666667d0 /
      data Tev2phQCD(12,9)  / 0.000000000000000d0 /
      data Tev2phQCD(12,10) / 0.000000000000000d0 /
      data Tev2phQCD(12,11) / 0.000000000000000d0 /
      data Tev2phQCD(12,12) /-0.100000000000000d0 /
      data Tev2phQCD(12,13) / 0.016666666666667d0 /
      data Tev2phQCD(13,0)  / 0.000000000000000d0 /
      data Tev2phQCD(13,1)  / 0.083333333333333d0 /
      data Tev2phQCD(13,2)  / 0.000000000000000d0 /
      data Tev2phQCD(13,3)  / 0.083333333333333d0 /
      data Tev2phQCD(13,4)  / 0.000000000000000d0 /
      data Tev2phQCD(13,5)  / 0.000000000000000d0 /
      data Tev2phQCD(13,6)  / 0.000000000000000d0 /
      data Tev2phQCD(13,7)  / 0.000000000000000d0 /
      data Tev2phQCD(13,8)  /-0.083333333333333d0 /
      data Tev2phQCD(13,9)  / 0.000000000000000d0 /
      data Tev2phQCD(13,10) / 0.000000000000000d0 /
      data Tev2phQCD(13,11) / 0.000000000000000d0 /
      data Tev2phQCD(13,12) / 0.000000000000000d0 /
      data Tev2phQCD(13,13) /-0.083333333333333d0 /
*
      data Tph2evQCD(0,0)   / 1d0 /
      data Tph2evQCD(0,1)   / 0d0 /
      data Tph2evQCD(0,2)   / 0d0 /
      data Tph2evQCD(0,3)   / 0d0 /
      data Tph2evQCD(0,4)   / 0d0 /
      data Tph2evQCD(0,5)   / 0d0 /
      data Tph2evQCD(0,6)   / 0d0 /
      data Tph2evQCD(0,7)   / 0d0 /
      data Tph2evQCD(0,8)   / 0d0 /
      data Tph2evQCD(0,9)   / 0d0 /
      data Tph2evQCD(0,10)  / 0d0 /
      data Tph2evQCD(0,11)  / 0d0 /
      data Tph2evQCD(0,12)  / 0d0 /
      data Tph2evQCD(0,13)  / 0d0 /
      data Tph2evQCD(1,0)   / 0d0 /
      data Tph2evQCD(1,1)   / 1d0 /
      data Tph2evQCD(1,2)   / 1d0 /
      data Tph2evQCD(1,3)   / 1d0 /
      data Tph2evQCD(1,4)   / 1d0 /
      data Tph2evQCD(1,5)   / 1d0 /
      data Tph2evQCD(1,6)   / 1d0 /
      data Tph2evQCD(1,7)   / 0d0 /
      data Tph2evQCD(1,8)   / 1d0 /
      data Tph2evQCD(1,9)   / 1d0 /
      data Tph2evQCD(1,10)  / 1d0 /
      data Tph2evQCD(1,11)  / 1d0 /
      data Tph2evQCD(1,12)  / 1d0 /
      data Tph2evQCD(1,13)  / 1d0 /
      data Tph2evQCD(2,0)   / 0d0 /
      data Tph2evQCD(2,1)   / 0d0 /
      data Tph2evQCD(2,2)   / 0d0 /
      data Tph2evQCD(2,3)   / 0d0 /
      data Tph2evQCD(2,4)   / 0d0 /
      data Tph2evQCD(2,5)   / 0d0 /
      data Tph2evQCD(2,6)   / 0d0 /
      data Tph2evQCD(2,7)   / 1d0 /
      data Tph2evQCD(2,8)   / 0d0 /
      data Tph2evQCD(2,9)   / 0d0 /
      data Tph2evQCD(2,10)  / 0d0 /
      data Tph2evQCD(2,11)  / 0d0 /
      data Tph2evQCD(2,12)  / 0d0 /
      data Tph2evQCD(2,13)  / 0d0 /
      data Tph2evQCD(3,0)   / 0d0 /
      data Tph2evQCD(3,1)   /-1d0 /
      data Tph2evQCD(3,2)   /-1d0 /
      data Tph2evQCD(3,3)   /-1d0 /
      data Tph2evQCD(3,4)   /-1d0 /
      data Tph2evQCD(3,5)   /-1d0 /
      data Tph2evQCD(3,6)   /-1d0 /
      data Tph2evQCD(3,7)   / 0d0 /
      data Tph2evQCD(3,8)   / 1d0 /
      data Tph2evQCD(3,9)   / 1d0 /
      data Tph2evQCD(3,10)  / 1d0 /
      data Tph2evQCD(3,11)  / 1d0 /
      data Tph2evQCD(3,12)  / 1d0 /
      data Tph2evQCD(3,13)  / 1d0 /
      data Tph2evQCD(4,0)   / 0d0 /
      data Tph2evQCD(4,1)   / 0d0 /
      data Tph2evQCD(4,2)   / 0d0 /
      data Tph2evQCD(4,3)   / 0d0 /
      data Tph2evQCD(4,4)   / 0d0 /
      data Tph2evQCD(4,5)   /-1d0 /
      data Tph2evQCD(4,6)   / 1d0 /
      data Tph2evQCD(4,7)   / 0d0 /
      data Tph2evQCD(4,8)   /-1d0 /
      data Tph2evQCD(4,9)   / 1d0 /
      data Tph2evQCD(4,10)  / 0d0 /
      data Tph2evQCD(4,11)  / 0d0 /
      data Tph2evQCD(4,12)  / 0d0 /
      data Tph2evQCD(4,13)  / 0d0 /
      data Tph2evQCD(5,0)   / 0d0 /
      data Tph2evQCD(5,1)   / 0d0 /
      data Tph2evQCD(5,2)   / 0d0 /
      data Tph2evQCD(5,3)   / 0d0 /
      data Tph2evQCD(5,4)   / 2d0 /
      data Tph2evQCD(5,5)   /-1d0 /
      data Tph2evQCD(5,6)   /-1d0 /
      data Tph2evQCD(5,7)   / 0d0 /
      data Tph2evQCD(5,8)   / 1d0 /
      data Tph2evQCD(5,9)   / 1d0 /
      data Tph2evQCD(5,10)  /-2d0 /
      data Tph2evQCD(5,11)  / 0d0 /
      data Tph2evQCD(5,12)  / 0d0 /
      data Tph2evQCD(5,13)  / 0d0 /
      data Tph2evQCD(6,0)   / 0d0 /
      data Tph2evQCD(6,1)   / 0d0 /
      data Tph2evQCD(6,2)   / 0d0 /
      data Tph2evQCD(6,3)   / 3d0 /
      data Tph2evQCD(6,4)   /-1d0 /
      data Tph2evQCD(6,5)   /-1d0 /
      data Tph2evQCD(6,6)   /-1d0 /
      data Tph2evQCD(6,7)   / 0d0 /
      data Tph2evQCD(6,8)   / 1d0 /
      data Tph2evQCD(6,9)   / 1d0 /
      data Tph2evQCD(6,10)  / 1d0 /
      data Tph2evQCD(6,11)  /-3d0 /
      data Tph2evQCD(6,12)  / 0d0 /
      data Tph2evQCD(6,13)  / 0d0 /
      data Tph2evQCD(7,0)   / 0d0 /
      data Tph2evQCD(7,1)   / 0d0 /
      data Tph2evQCD(7,2)   / 4d0 /
      data Tph2evQCD(7,3)   /-1d0 /
      data Tph2evQCD(7,4)   /-1d0 /
      data Tph2evQCD(7,5)   /-1d0 /
      data Tph2evQCD(7,6)   /-1d0 /
      data Tph2evQCD(7,7)   / 0d0 /
      data Tph2evQCD(7,8)   / 1d0 /
      data Tph2evQCD(7,9)   / 1d0 /
      data Tph2evQCD(7,10)  / 1d0 /
      data Tph2evQCD(7,11)  / 1d0 /
      data Tph2evQCD(7,12)  /-4d0 /
      data Tph2evQCD(7,13)  / 0d0 /
      data Tph2evQCD(8,0)   / 0d0 /
      data Tph2evQCD(8,1)   / 5d0 /
      data Tph2evQCD(8,2)   /-1d0 /
      data Tph2evQCD(8,3)   /-1d0 /
      data Tph2evQCD(8,4)   /-1d0 /
      data Tph2evQCD(8,5)   /-1d0 /
      data Tph2evQCD(8,6)   /-1d0 /
      data Tph2evQCD(8,7)   / 0d0 /
      data Tph2evQCD(8,8)   / 1d0 /
      data Tph2evQCD(8,9)   / 1d0 /
      data Tph2evQCD(8,10)  / 1d0 /
      data Tph2evQCD(8,11)  / 1d0 /
      data Tph2evQCD(8,12)  / 1d0 /
      data Tph2evQCD(8,13)  /-5d0 /
      data Tph2evQCD(9,0)   / 0d0 /
      data Tph2evQCD(9,1)   / 0d0 /
      data Tph2evQCD(9,2)   / 0d0 /
      data Tph2evQCD(9,3)   / 0d0 /
      data Tph2evQCD(9,4)   / 0d0 /
      data Tph2evQCD(9,5)   / 1d0 /
      data Tph2evQCD(9,6)   /-1d0 /
      data Tph2evQCD(9,7)   / 0d0 /
      data Tph2evQCD(9,8)   /-1d0 /
      data Tph2evQCD(9,9)   / 1d0 /
      data Tph2evQCD(9,10)  / 0d0 /
      data Tph2evQCD(9,11)  / 0d0 /
      data Tph2evQCD(9,12)  / 0d0 /
      data Tph2evQCD(9,13)  / 0d0 /
      data Tph2evQCD(10,0)  / 0d0 /
      data Tph2evQCD(10,1)  / 0d0 /
      data Tph2evQCD(10,2)  / 0d0 /
      data Tph2evQCD(10,3)  / 0d0 /
      data Tph2evQCD(10,4)  /-2d0 /
      data Tph2evQCD(10,5)  / 1d0 /
      data Tph2evQCD(10,6)  / 1d0 /
      data Tph2evQCD(10,7)  / 0d0 /
      data Tph2evQCD(10,8)  / 1d0 /
      data Tph2evQCD(10,9)  / 1d0 /
      data Tph2evQCD(10,10) /-2d0 /
      data Tph2evQCD(10,11) / 0d0 /
      data Tph2evQCD(10,12) / 0d0 /
      data Tph2evQCD(10,13) / 0d0 /
      data Tph2evQCD(11,0)  / 0d0 /
      data Tph2evQCD(11,1)  / 0d0 /
      data Tph2evQCD(11,2)  / 0d0 /
      data Tph2evQCD(11,3)  /-3d0 /
      data Tph2evQCD(11,4)  / 1d0 /
      data Tph2evQCD(11,5)  / 1d0 /
      data Tph2evQCD(11,6)  / 1d0 /
      data Tph2evQCD(11,7)  / 0d0 /
      data Tph2evQCD(11,8)  / 1d0 /
      data Tph2evQCD(11,9)  / 1d0 /
      data Tph2evQCD(11,10) / 1d0 /
      data Tph2evQCD(11,11) /-3d0 /
      data Tph2evQCD(11,12) / 0d0 /
      data Tph2evQCD(11,13) / 0d0 /
      data Tph2evQCD(12,0)  / 0d0 /
      data Tph2evQCD(12,1)  / 0d0 /
      data Tph2evQCD(12,2)  /-4d0 /
      data Tph2evQCD(12,3)  / 1d0 /
      data Tph2evQCD(12,4)  / 1d0 /
      data Tph2evQCD(12,5)  / 1d0 /
      data Tph2evQCD(12,6)  / 1d0 /
      data Tph2evQCD(12,7)  / 0d0 /
      data Tph2evQCD(12,8)  / 1d0 /
      data Tph2evQCD(12,9)  / 1d0 /
      data Tph2evQCD(12,10) / 1d0 /
      data Tph2evQCD(12,11) / 1d0 /
      data Tph2evQCD(12,12) /-4d0 /
      data Tph2evQCD(12,13) / 0d0 /
      data Tph2evQCD(13,0)  / 0d0 /
      data Tph2evQCD(13,1)  /-5d0 /
      data Tph2evQCD(13,2)  / 1d0 /
      data Tph2evQCD(13,3)  / 1d0 /
      data Tph2evQCD(13,4)  / 1d0 /
      data Tph2evQCD(13,5)  / 1d0 /
      data Tph2evQCD(13,6)  / 1d0 /
      data Tph2evQCD(13,7)  / 0d0 /
      data Tph2evQCD(13,8)  / 1d0 /
      data Tph2evQCD(13,9)  / 1d0 /
      data Tph2evQCD(13,10) / 1d0 /
      data Tph2evQCD(13,11) / 1d0 /
      data Tph2evQCD(13,12) / 1d0 /
      data Tph2evQCD(13,13) /-5d0 /


