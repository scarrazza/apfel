*     -*-fortran-*-
*
*     Switch for NLO QED corrections
*
      logical NLOQED
      character*4 InNLOQED

      logical SFNLOQED
      character*4 InSFNLOQED
*
      common / NLOQEDCorrectionsAPFEL / NLOQED,InNLOQED
      common / SFNLOQEDCorrectionsAPFEL / SFNLOQED,InSFNLOQED
