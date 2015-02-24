*     -*-fortran-*-
*
*     Switch for target mass corrections
*
      logical TMC
      character*4 InTMC
*
      double precision rhop
*
      common / TargetMassCorrectionsAPFEL / TMC,InTMC
      common / ProtonMassOverQAPFEL / rhop
