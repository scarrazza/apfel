#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

// Define APFEL I/O LHAPDF
#ifndef NOLHAPDF

// LHAPDF implementation
#include "LHAPDF/LHAPDF.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6

// global objects
LHAPDF::PDF* _pdfs;
extern "C" {
  void   mkpdfs_(int *mem, char* setname, int len)
  {
    LHAPDF::setVerbosity(0);
    std::string str(setname, len);

    if (_pdfs)
      delete _pdfs;
    try
      {
	_pdfs = LHAPDF::mkPDF(str, *mem);
      }
    catch(LHAPDF::Exception e)
      {
	std::cout << e.what() << std::endl;
      }
  }

  double xfxq_(int *fl, double *x, double *Q)
  {
    return _pdfs->xfxQ(*fl, *x, *Q);
  }

  bool islhapdf6_()
  {
    return true;
  }
}
#else

extern "C" {
  void   mkpdfs_(char* setname, int len)
  {
    std::string str(setname, len);
    LHAPDF::initPDFSetByName(str);
  }

  double xfxq_(int *mem, int *fl, double *x, double *Q)
  {
    LHAPDF::initPDF(*mem);
    std::vector<double> pdf = LHAPDF::xfx(*x, *Q);
    return pdf[*fl+6];
  }

  bool   islhapdf6_() { return false; }
}

#endif

#else
extern "C" void stop() { printf(" [Error] LHAPDF disabled.\n"); exit(-1); }
extern "C" void mkpdfs_(){ stop(); return ; }
extern "C" double xfxq_(){ stop(); return 0;}
extern "C" int numberpdf_(){ stop(); return 0;}
extern "C" bool islhapdf6_() { return true; }
#endif

// function to create folder for LHAPDFgrid.f
extern "C" void mkdir_(char* name, int) { mkdir(name, 0777); }
