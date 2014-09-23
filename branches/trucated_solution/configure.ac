AC_PREREQ(2.61)
AC_INIT([APFEL],[3.0.0],[stefano.carrazza@mi.infn.it],[apfel])
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported:\
 use a different --prefix argument])
fi
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
AM_MAINTAINER_MODE([disable])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


## Library version flags (at least update the version comment with each new release version)
VERSIONFLAGS="-version-info 3.0.0" ## for 3.0.0
AC_SUBST(VERSIONFLAGS)
AC_DEFINE_UNQUOTED(APFEL_VERSION, "$PACKAGE_VERSION", "APFEL version string")

## Checks and config headers
AC_LANG(C++)
#AC_CONFIG_SRCDIR([src/LHpdflib.F])
#AC_CONFIG_HEADERS([config/config.h include/LHAPDF/FortranWrappers.h include/LHAPDF/LHAPDFConfig.h])
AC_CONFIG_HEADERS([config/config.h include/APFEL/FortranWrappers.h])


## Set Fortran compiler behaviour
if test "x$FCFLAGS" == "x"; then
  FCFLAGS="-O3"
fi
# Try to respect users' Fortran compiler variables
if test "x$FC" == "x"; then
  if test "x$F77" == "x"; then
    FC="$GFORTRAN"
  else
    FC="$F77"
  fi
fi

#AC_LANG_PUSH(Fortran)
AC_PROG_FC
#AC_FC_SRCEXT(f90)
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
#AC_FC_FREEFORM
#AC_LANG_POP(Fortran)

AC_DEFUN(AC_PROG_LHAPDF, [AC_CHECK_PROG(LHAPDF,lhapdf-config,yes)])
AC_PROG_LHAPDF
if test x"${LHAPDF}" != x"yes" ; then
   AC_MSG_ERROR([Please install LHAPDF before trying to build.])
fi

## Enable binreloc
AM_BINRELOC

## Hack to make libtool use the FC Fortran compiler via F77, since
## libtool < 2.x doesn't understand --tag=FC
F77=$FC
AM_FFLAGS=$AM_FCFLAGS
FFLAGS="$FFLAGS $FCFLAGS"
AC_SUBST(FCFLAGS)
AC_SUBST(FFLAGS)
AC_SUBST(F77)
AC_SUBST(AM_FFLAGS)

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

## Enable debugging symbols.
AC_ARG_ENABLE([debug],
  [AC_HELP_STRING(--enable-debug, [build with debugging symbols (default=no)])],
  [], [enable_debug=no])
if test x$enable_debug == xyes; then
  [AM_FCFLAGS="$AM_FCFLAGS -g"]
fi


## Enable C++ wrapper
#AC_ARG_ENABLE([ccwrap],
#  [AC_HELP_STRING(--disable-ccwrap, [build without C++ wrapper])],
#  [], [enable_ccwrap=yes])
#if test x$enable_ccwrap == xyes; then
#   AC_MSG_NOTICE(Building C++ wrapper)
#else
#   AC_MSG_NOTICE(Not building C++ wrapper)
#fi
#AM_CONDITIONAL(ENABLE_CCWRAP, [test x$enable_ccwrap == xyes])


## Build Doxygen if possible
AC_ARG_ENABLE([doxygen],
  [AC_HELP_STRING(--disable-doxygen, [don't try to make Doxygen documentation])],
  [], [enable_doxygen=yes])
if test x$enable_doxygen == xyes; then
   AC_PATH_PROG(DOXYGEN, doxygen)
fi
AM_CONDITIONAL(WITH_DOXYGEN, [test x$DOXYGEN != x])


## Python extension
AC_ARG_ENABLE(pywrap, [AC_HELP_STRING(--disable-pywrap,
  [don't build Python module (default=build)])],
  [], [enable_pywrap=yes])   

enable_pywrap=yes
AZ_PYTHON_DEFAULT
## Basic Python checks
if test x$enable_pyext == xyes; then
  AZ_PYTHON_PATH
  AZ_PYTHON_VERSION_ENSURE([2.3])
 ## Test for Python header
  if test -x "$PYTHON"; then
    AC_MSG_CHECKING([for Python include path])
    python_incpath=`$PYTHON -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc();"`
    AC_MSG_RESULT([$python_incpath])
    python_header="$python_incpath/Python.h"
    if test -z "$python_incpath"; then
      AC_MSG_WARN([Can't build Python extension since include directory cannot be determined from distutils])
      enable_pywrap=no
    elif test ! -e "$python_header"; then
      AC_MSG_WARN([Can't build Python extension since header file $python_header cannot be found])
      enable_pywrap=no
    fi
  fi
fi

## SWIG version checks
if test x$enable_pywrap == xyes -a ! -e pywrap/apfel_wrap.cc; then
  AC_PATH_PROG(SWIG, swig)
  if test -x "$SWIG"; then
    swig_version=`$SWIG -version 2>&1 | grep Version | cut -f3 -d" "`
    swig_major=`echo $swig_version | cut -f1 -d.`
    swig_minor=`echo $swig_version | cut -f2 -d.`
    swig_micro=`echo $swig_version | cut -f3 -d.`
    let swig_version_num="1000*$swig_major+100*$swig_minor+$swig_micro"
    if test "$swig_version_num" -lt 1333; then
      AC_MSG_WARN([Can't build Python extension since SWIG version is too old to build wrapper files: please use a tarball release])
      enable_pywrap=no
    fi
  fi
fi
## Finish
AM_CONDITIONAL(ENABLE_PYWRAP, [test x$enable_pywrap == xyes])
if test x$enable_pywrap == xyes; then
  AC_MSG_NOTICE([All Python build checks successful: 'apfel' Python extension will be built])
fi


## Set final FCFLAGS, CXXFLAGS and CPPFLAGS
AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/include -I\$(top_builddir)/include"
AM_CPPFLAGS=["$AM_CPPFLAGS -DDATA_PATH="$datadir" -DAPFEL_VERSION="$PACKAGE_VERSION" "]
AC_SUBST(AM_CPPFLAGS)

#
#AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic "])
#AC_CEDAR_CHECKCXXFLAG([-ansi], [AM_CXXFLAGS="$AM_CXXFLAGS -ansi "])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall "])
AC_SUBST(AM_CXXFLAGS)
#
#AC_CEDAR_CHECKFCFLAG([-Wunused], [AM_FCFLAGS="$AM_FCFLAGS -Wunused "])
#AC_CEDAR_CHECKFCFLAG([-Wall], [AM_FCFLAGS="$AM_FCFLAGS -Wall "])
#AC_CEDAR_CHECKFCFLAG([-ffree-form], [AM_FCFLAGS="$AM_FCFLAGS -ffree-form "])
## For libtool F77 mode:
#AM_FFLAGS=$AM_FCFLAGS
#AC_SUBST(AM_FCFLAGS)
#AC_SUBST(AM_FFLAGS)


## Output
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([include/Makefile include/APFEL/Makefile])
AC_CONFIG_FILES([ccwrap/Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([examples/Makefile])
AC_CONFIG_FILES([pywrap/Makefile pywrap/setup.py])
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([bin/Makefile bin/apfel-config bin/apfel])
AC_CONFIG_FILES([apfelGUI/Makefile])
AC_OUTPUT
