*     -*-fortran-*-
*
*     Tranformations between QCD evolution and physical basis
*
*                   {{  120,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 },
*                    {    0,    10,     0,   -10,     0,     0,     0,     0,    10,     0,     0,     0,     0,   -10 },
*                    {    0,    10,     0,   -10,     0,     0,     0,    12,    -2,     0,     0,     0,   -12,     2 },
*                    {    0,    10,     0,   -10,     0,     0,    15,    -3,    -2,     0,     0,   -15,     3,     2 },
*                    {    0,    10,     0,   -10,     0,    20,    -5,    -3,    -2,     0,   -20,     5,     3,     2 },
*                    {    0,    10,     0,   -10,   -30,   -10,    -5,    -3,    -2,    30,    10,     5,     3,     2 },
*                    {    0,    10,     0,   -10,    30,   -10,    -5,    -3,    -2,   -30,    10,     5,     3,     2 },
* 120 * Tev2phQCD6 = {    0,     0,   120,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0 },
*                    {    0,    10,     0,    10,   -30,    10,     5,     3,     2,   -30,    10,     5,     3,     2 },
*                    {    0,    10,     0,    10,    30,    10,     5,     3,     2,    30,    10,     5,     3,     2 },
*                    {    0,    10,     0,    10,     0,   -20,     5,     3,     2,     0,   -20,     5,     3,     2 },
*                    {    0,    10,     0,    10,     0,     0,   -15,     3,     2,     0,     0,   -15,     3,     2 },
*                    {    0,    10,     0,    10,     0,     0,     0,   -12,     2,     0,     0,     0,   -12,     2 },
*                    {    0,    10,     0,    10,     0,     0,     0,     0,   -10,     0,     0,     0,     0,   -10 }}
*
*             {{   1,     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 },
*              {   0,     1,    1,    1,    1,    1,    1,    0,    1,    1,    1,    1,    1,    1 },
*              {   0,     0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0 },
*              {   0,    -1,   -1,   -1,   -1,   -1,   -1,    0,    1,    1,    1,    1,    1,    1 },
*              {   0,     0,    0,    0,    0,   -1,    1,    0,   -1,    1,    0,    0,    0,    0 },
*              {   0,     0,    0,    0,    2,   -1,   -1,    0,    1,    1,   -2,    0,    0,    0 },
*              {   0,     0,    0,    3,   -1,   -1,   -1,    0,    1,    1,    1,   -3,    0,    0 },
* Tph2evQCD6 = {   0,     0,    4,   -1,   -1,   -1,   -1,    0,    1,    1,    1,    1,   -4,    0 },
*              {   0,     5,   -1,   -1,   -1,   -1,   -1,    0,    1,    1,    1,    1,    1,   -5 },
*              {   0,     0,    0,    0,    0,    1,   -1,    0,   -1,    1,    0,    0,    0,    0 },
*              {   0,     0,    0,    0,   -2,    1,    1,    0,    1,    1,   -2,    0,    0,    0 },
*              {   0,     0,    0,   -3,    1,    1,    1,    0,    1,    1,    1,   -3,    0,    0 },
*              {   0,     0,   -4,    1,    1,    1,    1,    0,    1,    1,    1,    1,   -4,    0 },
*              {   0,    -5,    1,    1,    1,    1,    1,    0,    1,    1,    1,    1,    1,   -5 }}
*
      double precision Tev2phQCD(3:6,0:13,0:13)
      double precision Tph2evQCD(3:6,0:13,0:13)
*
      data Tev2phQCD(3, 0, 0)   /  1.00000000000000d0 /
      data Tev2phQCD(3, 0, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 0,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 1,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 2,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 3,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 1)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 4, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 3)   / -0.16666666666667d0 /
      data Tev2phQCD(3, 4, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 5)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 4, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4,10)   / -0.16666666666667d0 /
      data Tev2phQCD(3, 4,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 4,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5, 1)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 5, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5, 3)   / -0.16666666666667d0 /
      data Tev2phQCD(3, 5, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(3, 5, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(3, 5, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(3, 5,10)   /  0.08333333333333d0 /
      data Tev2phQCD(3, 5,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 5,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6, 1)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 6, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6, 3)   / -0.16666666666667d0 /
      data Tev2phQCD(3, 6, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(3, 6, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(3, 6, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(3, 6,10)   /  0.08333333333333d0 /
      data Tev2phQCD(3, 6,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 6,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 2)   /  1.00000000000000d0 /
      data Tev2phQCD(3, 7, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 7,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8, 1)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 8, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8, 3)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 8, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(3, 8, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(3, 8, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(3, 8,10)   /  0.08333333333333d0 /
      data Tev2phQCD(3, 8,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 8,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9, 1)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 9, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9, 3)   /  0.16666666666667d0 /
      data Tev2phQCD(3, 9, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(3, 9, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(3, 9, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(3, 9,10)   /  0.08333333333333d0 /
      data Tev2phQCD(3, 9,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3, 9,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 1)   /  0.16666666666667d0 /
      data Tev2phQCD(3,10, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 3)   /  0.16666666666667d0 /
      data Tev2phQCD(3,10, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 5)   / -0.16666666666667d0 /
      data Tev2phQCD(3,10, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10,10)   / -0.16666666666667d0 /
      data Tev2phQCD(3,10,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3,10,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3,11,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3,12,13)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13,10)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13,11)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13,12)   /  0.00000000000000d0 /
      data Tev2phQCD(3,13,13)   /  0.00000000000000d0 /
   
      data Tev2phQCD(4, 0, 0)   /  1.00000000000000d0 /
      data Tev2phQCD(4, 0, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0,11)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 0,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1,11)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 1,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2,11)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 2,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 3, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 3)   / -0.12500000000000d0 /
      data Tev2phQCD(4, 3, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 6)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 3, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3,11)   / -0.12500000000000d0 /
      data Tev2phQCD(4, 3,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 3,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 4, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4, 3)   / -0.12500000000000d0 /
      data Tev2phQCD(4, 4, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4, 5)   /  0.16666666666667d0 /
      data Tev2phQCD(4, 4, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(4, 4, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4,10)   / -0.16666666666667d0 /
      data Tev2phQCD(4, 4,11)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 4,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 4,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 5, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 5, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 5, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 5, 3)   / -0.12500000000000d0 /
      data Tev2phQCD(4, 5, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(4, 5, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(4, 5, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(4, 5, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 5, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 5, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(4, 5,10)   /  0.08333333333333d0 /
      data Tev2phQCD(4, 5,11)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 5,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 5,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 6, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 6, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 6, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 6, 3)   / -0.12500000000000d0 /
      data Tev2phQCD(4, 6, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(4, 6, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(4, 6, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(4, 6, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 6, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 6, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(4, 6,10)   /  0.08333333333333d0 /
      data Tev2phQCD(4, 6,11)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 6,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 6,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 2)   /  1.00000000000000d0 /
      data Tev2phQCD(4, 7, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7,11)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 7,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 8, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 8, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 8, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 8, 3)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 8, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(4, 8, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(4, 8, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 8, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 8, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 8, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(4, 8,10)   /  0.08333333333333d0 /
      data Tev2phQCD(4, 8,11)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 8,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 8,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 9, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 9, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 9, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 9, 3)   /  0.12500000000000d0 /
      data Tev2phQCD(4, 9, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(4, 9, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(4, 9, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 9, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 9, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 9, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(4, 9,10)   /  0.08333333333333d0 /
      data Tev2phQCD(4, 9,11)   /  0.04166666666667d0 /
      data Tev2phQCD(4, 9,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4, 9,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4,10, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10, 3)   /  0.12500000000000d0 /
      data Tev2phQCD(4,10, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10, 5)   / -0.16666666666667d0 /
      data Tev2phQCD(4,10, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(4,10, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10,10)   / -0.16666666666667d0 /
      data Tev2phQCD(4,10,11)   /  0.04166666666667d0 /
      data Tev2phQCD(4,10,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4,10,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 1)   /  0.12500000000000d0 /
      data Tev2phQCD(4,11, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 3)   /  0.12500000000000d0 /
      data Tev2phQCD(4,11, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 6)   / -0.12500000000000d0 /
      data Tev2phQCD(4,11, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11,11)   / -0.12500000000000d0 /
      data Tev2phQCD(4,11,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4,11,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12,11)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4,12,13)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13,10)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13,11)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13,12)   /  0.00000000000000d0 /
      data Tev2phQCD(4,13,13)   /  0.00000000000000d0 /
   
      data Tev2phQCD(5, 0, 0)   /  1.00000000000000d0 /
      data Tev2phQCD(5, 0, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0,11)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0,12)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 0,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1,11)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1,12)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 1,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 2, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 3)   / -0.10000000000000d0 /
      data Tev2phQCD(5, 2, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 7)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 2, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2,11)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 2,12)   / -0.10000000000000d0 /
      data Tev2phQCD(5, 2,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 3, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3, 3)   / -0.10000000000000d0 /
      data Tev2phQCD(5, 3, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3, 6)   /  0.12500000000000d0 /
      data Tev2phQCD(5, 3, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(5, 3, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 3,11)   / -0.12500000000000d0 /
      data Tev2phQCD(5, 3,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 3,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 4, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 4, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 4, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 4, 3)   / -0.10000000000000d0 /
      data Tev2phQCD(5, 4, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 4, 5)   /  0.16666666666667d0 /
      data Tev2phQCD(5, 4, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(5, 4, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(5, 4, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 4, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 4,10)   / -0.16666666666667d0 /
      data Tev2phQCD(5, 4,11)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 4,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 4,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 5, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 5, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 5, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 5, 3)   / -0.10000000000000d0 /
      data Tev2phQCD(5, 5, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(5, 5, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(5, 5, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(5, 5, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(5, 5, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 5, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(5, 5,10)   /  0.08333333333333d0 /
      data Tev2phQCD(5, 5,11)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 5,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 5,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 6, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 6, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 6, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 6, 3)   / -0.10000000000000d0 /
      data Tev2phQCD(5, 6, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(5, 6, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(5, 6, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(5, 6, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(5, 6, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 6, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(5, 6,10)   /  0.08333333333333d0 /
      data Tev2phQCD(5, 6,11)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 6,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 6,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 2)   /  1.00000000000000d0 /
      data Tev2phQCD(5, 7, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7,11)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7,12)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 7,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 8, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 8, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 8, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 8, 3)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 8, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(5, 8, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(5, 8, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 8, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 8, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 8, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(5, 8,10)   /  0.08333333333333d0 /
      data Tev2phQCD(5, 8,11)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 8,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 8,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 9, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 9, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 9, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 9, 3)   /  0.10000000000000d0 /
      data Tev2phQCD(5, 9, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(5, 9, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(5, 9, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 9, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 9, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5, 9, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(5, 9,10)   /  0.08333333333333d0 /
      data Tev2phQCD(5, 9,11)   /  0.04166666666667d0 /
      data Tev2phQCD(5, 9,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5, 9,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5,10, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5,10, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5,10, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5,10, 3)   /  0.10000000000000d0 /
      data Tev2phQCD(5,10, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5,10, 5)   / -0.16666666666667d0 /
      data Tev2phQCD(5,10, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(5,10, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(5,10, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5,10, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5,10,10)   / -0.16666666666667d0 /
      data Tev2phQCD(5,10,11)   /  0.04166666666667d0 /
      data Tev2phQCD(5,10,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5,10,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5,11, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11, 3)   /  0.10000000000000d0 /
      data Tev2phQCD(5,11, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11, 6)   / -0.12500000000000d0 /
      data Tev2phQCD(5,11, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(5,11, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5,11,11)   / -0.12500000000000d0 /
      data Tev2phQCD(5,11,12)   /  0.02500000000000d0 /
      data Tev2phQCD(5,11,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 1)   /  0.10000000000000d0 /
      data Tev2phQCD(5,12, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 3)   /  0.10000000000000d0 /
      data Tev2phQCD(5,12, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 7)   / -0.10000000000000d0 /
      data Tev2phQCD(5,12, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12,11)   /  0.00000000000000d0 /
      data Tev2phQCD(5,12,12)   / -0.10000000000000d0 /
      data Tev2phQCD(5,12,13)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13,10)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13,11)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13,12)   /  0.00000000000000d0 /
      data Tev2phQCD(5,13,13)   /  0.00000000000000d0 /
   
      data Tev2phQCD(6, 0, 0)   /  1.00000000000000d0 /
      data Tev2phQCD(6, 0, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0,11)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0,12)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 0,13)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 1, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 3)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 1, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1, 8)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 1, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1,11)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1,12)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 1,13)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 2, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 2, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2, 3)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 2, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2, 7)   /  0.10000000000000d0 /
      data Tev2phQCD(6, 2, 8)   / -0.01666666666667d0 /
      data Tev2phQCD(6, 2, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2,11)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 2,12)   / -0.10000000000000d0 /
      data Tev2phQCD(6, 2,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 3, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 3, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 3, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 3, 3)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 3, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 3, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 3, 6)   /  0.12500000000000d0 /
      data Tev2phQCD(6, 3, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(6, 3, 8)   / -0.01666666666667d0 /
      data Tev2phQCD(6, 3, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 3,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 3,11)   / -0.12500000000000d0 /
      data Tev2phQCD(6, 3,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 3,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 4, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 4, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 4, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 4, 3)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 4, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 4, 5)   /  0.16666666666667d0 /
      data Tev2phQCD(6, 4, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(6, 4, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(6, 4, 8)   / -0.01666666666667d0 /
      data Tev2phQCD(6, 4, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 4,10)   / -0.16666666666667d0 /
      data Tev2phQCD(6, 4,11)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 4,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 4,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 5, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 5, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 5, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 5, 3)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 5, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(6, 5, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 5, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(6, 5, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(6, 5, 8)   / -0.01666666666667d0 /
      data Tev2phQCD(6, 5, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(6, 5,10)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 5,11)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 5,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 5,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 6, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 6, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 6, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 6, 3)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 6, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(6, 6, 5)   / -0.08333333333333d0 /
      data Tev2phQCD(6, 6, 6)   / -0.04166666666667d0 /
      data Tev2phQCD(6, 6, 7)   / -0.02500000000000d0 /
      data Tev2phQCD(6, 6, 8)   / -0.01666666666667d0 /
      data Tev2phQCD(6, 6, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(6, 6,10)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 6,11)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 6,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 6,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 7, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 1)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 2)   /  1.00000000000000d0 /
      data Tev2phQCD(6, 7, 3)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 8)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7,11)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7,12)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 7,13)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 8, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 8, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 8, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 8, 3)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 8, 4)   / -0.25000000000000d0 /
      data Tev2phQCD(6, 8, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 8, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 8, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 8, 8)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 8, 9)   / -0.25000000000000d0 /
      data Tev2phQCD(6, 8,10)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 8,11)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 8,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 8,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 9, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 9, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 9, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6, 9, 3)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 9, 4)   /  0.25000000000000d0 /
      data Tev2phQCD(6, 9, 5)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 9, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 9, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 9, 8)   /  0.01666666666667d0 /
      data Tev2phQCD(6, 9, 9)   /  0.25000000000000d0 /
      data Tev2phQCD(6, 9,10)   /  0.08333333333333d0 /
      data Tev2phQCD(6, 9,11)   /  0.04166666666667d0 /
      data Tev2phQCD(6, 9,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6, 9,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6,10, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6,10, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6,10, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6,10, 3)   /  0.08333333333333d0 /
      data Tev2phQCD(6,10, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6,10, 5)   / -0.16666666666667d0 /
      data Tev2phQCD(6,10, 6)   /  0.04166666666667d0 /
      data Tev2phQCD(6,10, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(6,10, 8)   /  0.01666666666667d0 /
      data Tev2phQCD(6,10, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6,10,10)   / -0.16666666666667d0 /
      data Tev2phQCD(6,10,11)   /  0.04166666666667d0 /
      data Tev2phQCD(6,10,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6,10,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6,11, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6,11, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6,11, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6,11, 3)   /  0.08333333333333d0 /
      data Tev2phQCD(6,11, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6,11, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6,11, 6)   / -0.12500000000000d0 /
      data Tev2phQCD(6,11, 7)   /  0.02500000000000d0 /
      data Tev2phQCD(6,11, 8)   /  0.01666666666667d0 /
      data Tev2phQCD(6,11, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6,11,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6,11,11)   / -0.12500000000000d0 /
      data Tev2phQCD(6,11,12)   /  0.02500000000000d0 /
      data Tev2phQCD(6,11,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6,12, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6,12, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12, 3)   /  0.08333333333333d0 /
      data Tev2phQCD(6,12, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12, 7)   / -0.10000000000000d0 /
      data Tev2phQCD(6,12, 8)   /  0.01666666666667d0 /
      data Tev2phQCD(6,12, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12,11)   /  0.00000000000000d0 /
      data Tev2phQCD(6,12,12)   / -0.10000000000000d0 /
      data Tev2phQCD(6,12,13)   /  0.01666666666667d0 /
      data Tev2phQCD(6,13, 0)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13, 1)   /  0.08333333333333d0 /
      data Tev2phQCD(6,13, 2)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13, 3)   /  0.08333333333333d0 /
      data Tev2phQCD(6,13, 4)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13, 5)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13, 6)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13, 7)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13, 8)   / -0.08333333333333d0 /
      data Tev2phQCD(6,13, 9)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13,10)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13,11)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13,12)   /  0.00000000000000d0 /
      data Tev2phQCD(6,13,13)   / -0.08333333333333d0 /
   
      data Tph2evQCD(3, 0, 0)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 0, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0,10)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 0,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 1, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 1, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 1, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 1, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 1,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 1,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 1,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 7)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 2, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2,10)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 2,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 3, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 3, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 3, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 3, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 3,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 3,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 3,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 4, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 4, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 4, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 4,10)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 4,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5, 4)   /  2.00000000000000d0 /
      data Tph2evQCD(3, 5, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 5, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 5, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 5, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 5,10)   / -2.00000000000000d0 /
      data Tph2evQCD(3, 5,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 5,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 6, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 6, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 6, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 6, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 6,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 6,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 6,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 7, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 7, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 7, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 7, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 7,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 7,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 7,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 8, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 8, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 8, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 8, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 8,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 8,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 8,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 9, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 9, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(3, 9, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3, 9,10)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3, 9,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10, 4)   / -2.00000000000000d0 /
      data Tph2evQCD(3,10, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(3,10, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(3,10, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3,10, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3,10,10)   / -2.00000000000000d0 /
      data Tph2evQCD(3,10,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3,10,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(3,11, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(3,11, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(3,11, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3,11, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3,11,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3,11,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3,11,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(3,12, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(3,12, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(3,12, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3,12, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3,12,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3,12,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3,12,13)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(3,13, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(3,13, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(3,13, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(3,13, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(3,13,10)   /  1.00000000000000d0 /
      data Tph2evQCD(3,13,11)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13,12)   /  0.00000000000000d0 /
      data Tph2evQCD(3,13,13)   /  0.00000000000000d0 /
   
      data Tph2evQCD(4, 0, 0)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 0, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0,10)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0,11)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 0,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 1, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 1, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 1, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 1, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 1, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1,11)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 1,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 1,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 7)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 2, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2,10)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2,11)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 2,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 3, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 3, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 3, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 3, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 3, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 3, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 3, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 3, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 3, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 3, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 3,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 3,11)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 3,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 3,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 4, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 4, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 4, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 4,10)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4,11)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 4,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5, 4)   /  2.00000000000000d0 /
      data Tph2evQCD(4, 5, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 5, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 5, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 5, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 5,10)   / -2.00000000000000d0 /
      data Tph2evQCD(4, 5,11)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 5,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 6, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 6, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 6, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 6, 3)   /  3.00000000000000d0 /
      data Tph2evQCD(4, 6, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 6, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 6, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 6, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 6, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 6, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 6,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 6,11)   / -3.00000000000000d0 /
      data Tph2evQCD(4, 6,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 6,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 7, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 7, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 7, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 7, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 7, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 7, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 7, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 7, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 7, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 7, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 7,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 7,11)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 7,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 7,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 8, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 8, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 8, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 8, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 8, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 8, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 8, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 8, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 8, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 8, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 8,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 8,11)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 8,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 8,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 9, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 9, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(4, 9, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4, 9,10)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9,11)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4, 9,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10, 4)   / -2.00000000000000d0 /
      data Tph2evQCD(4,10, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(4,10, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(4,10, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4,10, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4,10,10)   / -2.00000000000000d0 /
      data Tph2evQCD(4,10,11)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4,10,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4,11, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4,11, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4,11, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4,11, 3)   / -3.00000000000000d0 /
      data Tph2evQCD(4,11, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(4,11, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(4,11, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(4,11, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4,11, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4,11, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4,11,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4,11,11)   / -3.00000000000000d0 /
      data Tph2evQCD(4,11,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4,11,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4,12, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4,12, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4,12, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4,12, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4,12, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12,11)   /  1.00000000000000d0 /
      data Tph2evQCD(4,12,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4,12,13)   /  0.00000000000000d0 /
      data Tph2evQCD(4,13, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(4,13, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(4,13, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(4,13, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(4,13, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13,10)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13,11)   /  1.00000000000000d0 /
      data Tph2evQCD(4,13,12)   /  0.00000000000000d0 /
      data Tph2evQCD(4,13,13)   /  0.00000000000000d0 /
   
      data Tph2evQCD(5, 0, 0)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 0, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0,10)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0,11)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 0,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 1, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 1, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 1, 2)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 1, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1,11)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1,12)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 1,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 7)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 2, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2,10)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2,11)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 2,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 3, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 3, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 3, 2)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 3, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 3, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 3, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 3, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 3, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 3, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 3, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 3,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 3,11)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 3,12)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 3,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 4, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 4, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 4, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 4,10)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4,11)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 4,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5, 4)   /  2.00000000000000d0 /
      data Tph2evQCD(5, 5, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 5, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 5, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 5, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 5,10)   / -2.00000000000000d0 /
      data Tph2evQCD(5, 5,11)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 5,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 6, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 6, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 6, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 6, 3)   /  3.00000000000000d0 /
      data Tph2evQCD(5, 6, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 6, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 6, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 6, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 6, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 6, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 6,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 6,11)   / -3.00000000000000d0 /
      data Tph2evQCD(5, 6,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 6,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 7, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 7, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 7, 2)   /  4.00000000000000d0 /
      data Tph2evQCD(5, 7, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 7, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 7, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 7, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 7, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 7, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 7, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 7,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 7,11)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 7,12)   / -4.00000000000000d0 /
      data Tph2evQCD(5, 7,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 8, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 8, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 8, 2)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 8, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 8, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 8, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 8, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 8, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 8, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 8, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 8,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 8,11)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 8,12)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 8,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 9, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 9, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(5, 9, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5, 9,10)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9,11)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5, 9,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10, 4)   / -2.00000000000000d0 /
      data Tph2evQCD(5,10, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(5,10, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(5,10, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5,10, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5,10,10)   / -2.00000000000000d0 /
      data Tph2evQCD(5,10,11)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5,10,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5,11, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5,11, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5,11, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(5,11, 3)   / -3.00000000000000d0 /
      data Tph2evQCD(5,11, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(5,11, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(5,11, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(5,11, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5,11, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5,11, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5,11,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5,11,11)   / -3.00000000000000d0 /
      data Tph2evQCD(5,11,12)   /  0.00000000000000d0 /
      data Tph2evQCD(5,11,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5,12, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5,12, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5,12, 2)   / -4.00000000000000d0 /
      data Tph2evQCD(5,12, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5,12, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12,11)   /  1.00000000000000d0 /
      data Tph2evQCD(5,12,12)   / -4.00000000000000d0 /
      data Tph2evQCD(5,12,13)   /  0.00000000000000d0 /
      data Tph2evQCD(5,13, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(5,13, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(5,13, 2)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(5,13, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13,10)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13,11)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13,12)   /  1.00000000000000d0 /
      data Tph2evQCD(5,13,13)   /  0.00000000000000d0 /
   
      data Tph2evQCD(6, 0, 0)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 0, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0,10)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0,11)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 0,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 1, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 1, 1)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 2)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 1, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1,11)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1,12)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 1,13)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 2, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 5)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 6)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 7)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 2, 8)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2, 9)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2,10)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2,11)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 2,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 3, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 3, 1)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 3, 2)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 3, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 3, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 3, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 3, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 3, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 3, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 3, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 3,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 3,11)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 3,12)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 3,13)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 4, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 4, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 4, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 4, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 4,10)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4,11)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 4,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5, 4)   /  2.00000000000000d0 /
      data Tph2evQCD(6, 5, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 5, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 5, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 5, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 5,10)   / -2.00000000000000d0 /
      data Tph2evQCD(6, 5,11)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 5,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 6, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 6, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 6, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 6, 3)   /  3.00000000000000d0 /
      data Tph2evQCD(6, 6, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 6, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 6, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 6, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 6, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 6, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 6,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 6,11)   / -3.00000000000000d0 /
      data Tph2evQCD(6, 6,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 6,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 7, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 7, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 7, 2)   /  4.00000000000000d0 /
      data Tph2evQCD(6, 7, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 7, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 7, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 7, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 7, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 7, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 7, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 7,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 7,11)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 7,12)   / -4.00000000000000d0 /
      data Tph2evQCD(6, 7,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 8, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 8, 1)   /  5.00000000000000d0 /
      data Tph2evQCD(6, 8, 2)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 8, 3)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 8, 4)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 8, 5)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 8, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 8, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 8, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 8, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 8,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 8,11)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 8,12)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 8,13)   / -5.00000000000000d0 /
      data Tph2evQCD(6, 9, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9, 4)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 9, 6)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 9, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9, 8)   / -1.00000000000000d0 /
      data Tph2evQCD(6, 9, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6, 9,10)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9,11)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6, 9,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10, 3)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10, 4)   / -2.00000000000000d0 /
      data Tph2evQCD(6,10, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(6,10, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(6,10, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6,10, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6,10,10)   / -2.00000000000000d0 /
      data Tph2evQCD(6,10,11)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6,10,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6,11, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6,11, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6,11, 2)   /  0.00000000000000d0 /
      data Tph2evQCD(6,11, 3)   / -3.00000000000000d0 /
      data Tph2evQCD(6,11, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(6,11, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(6,11, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(6,11, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6,11, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6,11, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6,11,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6,11,11)   / -3.00000000000000d0 /
      data Tph2evQCD(6,11,12)   /  0.00000000000000d0 /
      data Tph2evQCD(6,11,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6,12, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6,12, 1)   /  0.00000000000000d0 /
      data Tph2evQCD(6,12, 2)   / -4.00000000000000d0 /
      data Tph2evQCD(6,12, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6,12, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12,11)   /  1.00000000000000d0 /
      data Tph2evQCD(6,12,12)   / -4.00000000000000d0 /
      data Tph2evQCD(6,12,13)   /  0.00000000000000d0 /
      data Tph2evQCD(6,13, 0)   /  0.00000000000000d0 /
      data Tph2evQCD(6,13, 1)   / -5.00000000000000d0 /
      data Tph2evQCD(6,13, 2)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13, 3)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13, 4)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13, 5)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13, 6)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13, 7)   /  0.00000000000000d0 /
      data Tph2evQCD(6,13, 8)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13, 9)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13,10)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13,11)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13,12)   /  1.00000000000000d0 /
      data Tph2evQCD(6,13,13)   / -5.00000000000000d0 /




