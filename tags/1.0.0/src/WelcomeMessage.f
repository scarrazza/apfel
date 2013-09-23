************************************************************************
*
*     WelcomeMessage.f:
*
*     This routine returns the welcome message.
*
************************************************************************
      subroutine WelcomeMessage
*
      character*6 apfelversion

      call getapfelversion(apfelversion)

      write(6,*) " "
      write(6,*) "Welcome to "
      write(6,*) "     _/_/_/    _/_/_/_/   _/_/_/_/   _/_/_/_/   _/"
      write(6,*) "   _/    _/   _/    _/   _/         _/         _/"
      write(6,*) "  _/_/_/_/   _/_/_/_/   _/_/_/     _/_/_/     _/"
      write(6,*) " _/    _/   _/         _/         _/         _/"
      write(6,*) "_/    _/   _/         _/         _/_/_/_/   _/_/_/_/"
      write(6,*) "_____v", apfelversion,
     1     "A Partonic Function Evolution Library"      
      write(6,*) " "
*
      return
      end
