AC_PREREQ([2.69])
AC_INIT([APFEL],[3.0.7],[valerio.bertone@cern.ch,stefano.carrazza@cern.ch],[apfel])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE

AC_DEFINE_UNQUOTED(APFEL_VERSION, "$PACKAGE_VERSION", "APFEL version string")

AC_PROG_CXX
AC_CONFIG_HEADERS([config/config.h include/APFEL/FortranWrappers.h])

## Set Fortran compiler behaviour
if test "x$FCFLAGS" == "x"; then
#  FFLAGS="-fopenmp -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -fopenmp -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fvar-tracking-assignments -pipe -fsanitize=address"
# CXXFLAGS="-fsanitize=address"
# CCFLAGS="-fsanitize=address"
  FCFLAGS="-g -O3 -Wunused -cpp"
fi
# Try to respect users' Fortran compiler variables
if test "x$FC" == "x"; then
  if test "x$F77" == "x"; then
    FC="$GFORTRAN"
  else
    FC="$F77"
  fi
fi

AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported:\
 use a different --prefix argument])
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

LT_INIT

## Enable debugging symbols.
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],[build with debugging symbols (default=no)])],
  [], [enable_debug=no])
if test x$enable_debug == xyes; then
  [AM_FCFLAGS="$AM_FCFLAGS -g"]
fi

# Check host
AC_CANONICAL_HOST
case $host_os in
  darwin*)
	AM_CPPFLAGS=["$AM_CPPFLAGS -DDARWIN"]
        ;;
  linux*)
	#AM_CPPFLAGS=["$AM_CPPFLAGS -DLINUX"]
	;;
esac

## Build Doxygen if possible
AC_ARG_ENABLE([lhapdf],
  [AS_HELP_STRING([--disable-lhapdf],[don't link to LHAPDF])],
  [], [enable_lhapdf=yes])
if test x$enable_lhapdf == xyes; then
  AC_SEARCH_LHAPDF  
else
  AM_CPPFLAGS=["$AM_CPPFLAGS -DNOLHAPDF"]  
fi
AM_CONDITIONAL(ENABLE_LHAPDF, [test x$enable_lhapdf == xyes])

## Python extension
AC_ARG_ENABLE(pywrap, [AS_HELP_STRING([--disable-pywrap],[don't build Python module (default=build)])],
  [], [enable_pywrap=yes])   

enable_pywrap=yes
AZ_PYTHON_DEFAULT
## Basic Python checks
if test x$enable_pywrap == xyes; then
  AZ_PYTHON_PATH
  #AZ_PYTHON_VERSION_ENSURE([2.3])
 ## Test for Python header
  if test -x "$PYTHON"; then
    AC_MSG_CHECKING([for Python include path])
    python_incpath=`$PYTHON -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"`
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

AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall -std=c++11 "])
AC_SUBST(AM_CXXFLAGS)

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
AC_OUTPUT
