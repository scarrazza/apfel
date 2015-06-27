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
*
      call getapfelversion(apfelversion)
*
      write(6,*) achar(27)//"[1;32m"
      write(6,*) "Welcome to "
      write(6,*) "     _/_/_/    _/_/_/_/   _/_/_/_/   _/_/_/_/   _/"
      write(6,*) "   _/    _/   _/    _/   _/         _/         _/"
      write(6,*) "  _/_/_/_/   _/_/_/_/   _/_/_/     _/_/_/     _/"
      write(6,*) " _/    _/   _/         _/         _/         _/"
      write(6,*) "_/    _/   _/         _/         _/_/_/_/   _/_/_/_/"
      write(6,*) "_____v", apfelversion,
     1     "A PDF Evolution Library, arXiv:1310.1394"      
      write(6,*) "     Authors: V. Bertone, S. Carrazza, J. Rojo"
      write(6,*) achar(27)//"[0m"
*
      return
      end
