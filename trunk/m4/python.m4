# ===========================================================================
#               http://autoconf-archive.cryp.to/az_python.html
# ===========================================================================
#
# SYNOPSIS
#
#   AZ_PYTHON_DEFAULT
#   AZ_PYTHON_ENABLE
#   AZ_PYTHON_WITH
#   AZ_PYTHON_PATH
#   AZ_PYTHON_VERSION_ENSURE( [2.2] )
#   AZ_PYTHON_CSPEC
#   AZ_PYTHON_LSPEC
#
# DESCRIPTION
#
#   This file provides autoconf support for those applications that want to
#   embed python. It supports all pythons >= 2.2 which is the first official
#   release containing distutils. Version 2.2 of python was released
#   December 21, 2001. Since it actually executes the python, cross platform
#   configuration will probably not work. Also, most of the platforms
#   supported are consistent until you look into MacOSX. The python included
#   with it is installed as a framework which is a very different
#   environment to set up the normal tools such as gcc and libtool to deal
#   with. Therefore, once we establish which python that we are going to
#   use, we use its distutils to actually compile and link our modules or
#   applications.
#
#   At this time, it does NOT support linking with Python statically. It
#   does support dynamic linking.
#
#   This set of macros help define $PYTHON, $PYTHON_USE, $PYTHON_CSPEC and
#   $PYTHON_LSPEC. $PYTHON defines the full executable path for the Python
#   being linked to and is used within these macros to determine if that has
#   been specified or found. These macros do execute this python version so
#   it must be present on the system at configure time.
#
#   $PYTHON_USE is an automake variable that defines whether Python support
#   should be included or not in your application. $PYTHON_CSPEC is a
#   variable that supplies additional CFLAGS for the compilation of the
#   application/shared library. $PYTHON_LSPEC is a variable that supplies
#   additional LDFLAGS for linking the application/shared library.
#
#   The following is an example of how to set up for python usage within
#   your application in your configure.in:
#
#     AZ_PYTHON_DEFAULT( )
#     AZ_PYTHON_ENABLE( )             # Optional
#     AZ_PYTHON_WITH( )               # Optional
#     AZ_PYTHON_PATH( )               # or AZ_PYTHON_INSIST( )
#     # if $PYTHON is not defined, then the following do nothing.
#     AZ_PYTHON_VERSION_ENSURE( [2.2] )
#     AZ_PYTHON_CSPEC
#     AZ_PYTHON_LSPEC
#
#   The AZ_PYTHON_DEFAULT sets the $PYTHON_USE to false. Thereby, excluding
#   it if it was optional.
#
#   The AZ_PYTHON_ENABLE looks for the optional configure parameters of
#   --enable-python/--disable-python and establishes the $PYTHON and
#   $PYTHON_USE variables accordingly.
#
#   The AZ_PYTHON_WITH looks for the optional configure parameters of
#   --with-python/--without-python and establishes the $PYTHON and
#   $PYTHON_USE variables accordingly.
#
#   The AZ_PYTHON_PATH looks for python assuming that none has been
#   previously found or defined and issues an error if it does not find it.
#   If it does find it, it establishes the $PYTHON and $PYTHON_USE variables
#   accordingly. AZ_PYTHON_INSIST could be used here instead if you want to
#   insist that Python support be included using the --enable-python or
#   --with-python checks previously done.
#
#   The AZ_PYTHON_VERSION_ENSURE issues an error if the Python previously
#   found is not of version 2.2 or greater.
#
#   Once that these macros have be run, we can use PYTHON_USE within the
#   makefile.am file to conditionally add the Python support such as:
#
#   Makefile.am example showing optional inclusion of directories:
#
#    if PYTHON_USE
#    plugins = plugins
#    src = src
#    else
#    plugins =
#    src =
#    endif
#
#    SUBDIRS = . $(plugins) $(src)
#
#   Makefile.am example showing optional shared library build:
#
#    if PYTHON_USE
#    lib_LTLIBRARIES        = libElemList.la
#    libElemList_la_SOURCES = libElemList.c
#    libElemList_la_CFLAGS  = @PYTHON_CSPEC@
#    libElemList_la_LDFLAGS = @PYTHON_LSPEC@
#    endif
#
#   Makefile.am example showing optional program build:
#
#    if PYTHON_USE
#    bin_PROGRAMS    = runFunc
#    runFunc_SOURCES = runFunc.c
#    runFunc_CFLAGS  = @PYTHON_CSPEC@
#    runFunc_LDFLAGS = @PYTHON_LSPEC@
#    endif
#
#   The above compiles the modules only if PYTHON_USE was specified as true.
#   Also, the else portion of the if was optional.
#
# LAST MODIFICATION
#
#   2008-04-12
#
# COPYLEFT
#
#   Copyright (c) 2008 Robert White <kranki@mac.com>
#   Copyright (c) 2008 Dustin J. Mitchell <dustin@cs.uchicago.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

# AZ_PYTHON_DEFAULT( )
# -----------------
# Sets the default to not include Python support.

AC_DEFUN([AZ_PYTHON_DEFAULT],
[
    az_python_use=false
    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
])



# AZ_PYTHON_ENABLE( [path] )
# -----------------------------------------------------------------
# Handles the various --enable-python commands.
# Input:
#   $1 is the optional search path for the python executable if needed
# Ouput:
#   PYTHON_USE (AM_CONDITIONAL) is true if python executable found
#   and --enable-python was requested; otherwise false.
#   $PYTHON contains the full executable path to python if PYTHON_ENABLE_USE
#   is true.
#
# Example:
#   AZ_PYTHON_ENABLE( )
#   or
#   AZ_PYTHON_ENABLE( "/usr/bin" )

AC_DEFUN([AZ_PYTHON_ENABLE],
[
    AC_ARG_VAR([PYTHON],[Python Executable Path])

    # unless PYTHON was supplied to us (as a precious variable),
    # see if --enable-python[=PythonExecutablePath], --enable-python,
    # --disable-python or --enable-python=no was given.
    if test -z "$PYTHON"
    then
        AC_MSG_CHECKING(for --enable-python)
        AC_ARG_ENABLE(
            python,
            AC_HELP_STRING([--enable-python@<:@=PYTHON@:>@],
                [absolute path name of Python executable]
            ),
            [
                if test "$enableval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG(
                        [PYTHON],
                        python,
                        [],
                        $1
                    )
                    if test -z "$PYTHON"
                    then
                        AC_MSG_ERROR(no path to python found)
                    fi
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                elif test "$enableval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_python_use=false
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                else
                    # $enableval must be the executable path then.
                    AC_SUBST([PYTHON], ["${enableval}"])
                    AC_MSG_RESULT($withval)
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                fi
            ],
            [
                # --with-python was not specified.
                AC_MSG_RESULT(no)
                az_python_use=false
                AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
            ]
        )
    fi

])



# AZ_PYTHON_CSPEC( )
# -----------------
# Set up the c compiler options to compile Python
# embedded programs/libraries in $PYTHON_CSPEC if
# $PYTHON has been defined.

AC_DEFUN([AZ_PYTHON_CSPEC],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        az_python_prefix=`${PYTHON} -c "import sys; print sys.prefix"`
        if test -z "$az_python_prefix"
        then
            AC_MSG_ERROR([Python Prefix is not known])
        fi
        az_python_execprefix=`${PYTHON} -c "import sys; print sys.exec_prefix"`
        az_python_version=`$PYTHON -c "import sys; print sys.version[[:3]]"`
        az_python_includespec="-I${az_python_prefix}/include/python${az_python_version}"
        if test x"$python_prefix" != x"$python_execprefix"; then
            az_python_execspec="-I${az_python_execprefix}/include/python${az_python_version}"
            az_python_includespec="${az_python_includespec} $az_python_execspec"
        fi
        az_python_ccshared=`${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_config_var('CFLAGSFORSHARED')"`
        az_python_cspec="${az_python_ccshared} ${az_python_includespec}"
        AC_SUBST([PYTHON_CSPEC], [${az_python_cspec}])
        AC_MSG_NOTICE([PYTHON_CSPEC=${az_python_cspec}])
    fi
])



# AZ_PYTHON_INSIST( )
# -----------------
# Look for Python and set the output variable 'PYTHON'
# to 'python' if found, empty otherwise.

AC_DEFUN([AZ_PYTHON_PATH],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    fi
])



# AZ_PYTHON_LSPEC( )
# -----------------
# Set up the linker options to link Python embedded
# programs/libraries in $PYTHON_LSPEC if $PYTHON
# has been defined.

AC_DEFUN([AZ_PYTHON_LSPEC],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        AZ_PYTHON_RUN([
import sys
import distutils.sysconfig
strUseFrameWork = "--enable-framework"
dictConfig = distutils.sysconfig.get_config_vars( )
strConfigArgs = dictConfig.get("CONFIG_ARGS")
strLinkSpec =  dictConfig.get('LDFLAGS')
if -1 ==  strConfigArgs.find(strUseFrameWork):
    strLibPL = dictConfig.get("LIBPL")
    if strLibPL and (strLibPL != ""):
        strLinkSpec += " -L%s" % (strLibPL)
    strSys = dictConfig.get("SYSLIBS")
    if strSys and (strSys != ""):
        strLinkSpec += " %s" % (strSys)
    strSHL = dictConfig.get("SHLIBS")
    if strSHL and (strSHL != ""):
        strLinkSpec += " %s" % (strSHL)
    # Construct the Python Library Name.
    strTmplte = " -lpython%d.%d"
    if (sys.platform == "win32") or (sys.platform == "os2emx"):
        strTmplte = " -lpython%d%d"
    strWrk = strTmplte % ( (sys.hexversion >> 24),
                            ((sys.hexversion >> 16) & 0xff))
    strLinkSpec += strWrk
else:
    # This is not ideal since it changes the search path
    # for Frameworks which could have side-effects on
    # other included Frameworks.  However, it is necessary
    # where someone has installed more than one frameworked
    # Python.  Frameworks are really only used in MacOSX.
    strLibFW = dictConfig.get("PYTHONFRAMEWORKPREFIX")
    if strLibFW and (strLibFW != ""):
        strLinkSpec += " -F%s" % (strLibFW)
strLinkSpec += " %s" % (dictConfig.get('LINKFORSHARED'))
print strLinkSpec
        ])
        AC_SUBST([PYTHON_LSPEC], [${az_python_output}])
        AC_MSG_NOTICE([PYTHON_LSPEC=${az_python_output}])
    fi
])



# AZ_PYTHON_PATH( )
# -----------------
# Look for Python and set the output variable 'PYTHON'
# to 'python' if found, empty otherwise.

AC_DEFUN([AZ_PYTHON_PATH],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    AC_PATH_PROG( PYTHON, python, [], $1 )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        az_python_use=true
    fi
    AM_CONDITIONAL(PYTHON_USE, test "$az_python_use" = "true")
])



# AZ_PYTHON_PREFIX( )
# -------------------
# Use the values of $prefix and $exec_prefix for the corresponding
# values of PYTHON_PREFIX and PYTHON_EXEC_PREFIX.

AC_DEFUN([AZ_PYTHON_PREFIX],
[
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable Path is not known])
    fi
    ax_python_prefix=`${PYTHON} -c "import sys; print sys.prefix"`
    ax_python_execprefix=`${PYTHON} -c "import sys; print sys.exec_prefix"`
    AC_SUBST([PYTHON_PREFIX], ["${ax_python_prefix}"])
    AC_SUBST([PYTHON_EXECPREFIX], ["${ax_python_execprefix}"])
])



# AZ_PYTHON_RUN( PYTHON_PROGRAM )
# -----------------
# Run a Python Test Program saving its output
# in az_python_output and its condition code
# in az_python_cc.

AC_DEFUN([AZ_PYTHON_RUN],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        cat >conftest.py <<_ACEOF
$1
_ACEOF
        az_python_output=`$PYTHON conftest.py`
        az_python_cc=$?
        rm conftest.py
        if test -f "conftest.pyc"
        then
            rm conftest.pyc
        fi
    fi
])



# AZ_PYTHON_VERSION_CHECK( VERSION, [ACTION-IF-TRUE], [ACTION-IF-FALSE] )
# -----------------------------------------------------------------------------
# Run ACTION-IF-TRUE if the Python interpreter has version >= VERSION.
# Run ACTION-IF-FALSE otherwise.
# This test uses sys.hexversion instead of the string equivalant (first
# word of sys.version), in order to cope with versions such as 2.2c1.
# hexversion has been introduced in Python 1.5.2; it's probably not
# worth to support older versions (1.5.1 was released on October 31, 1998).

AC_DEFUN([AZ_PYTHON_VERSION_CHECK],
 [
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        AC_MSG_CHECKING([whether $PYTHON version >= $1])
        AZ_PYTHON_RUN([
import sys, string
# split strings by '.' and convert to numeric.  Append some zeros
# because we need at least 4 digits for the hex conversion.
minver = map(int, string.split('$1', '.')) + [[0, 0, 0]]
minverhex = 0
for i in xrange(0, 4): minverhex = (minverhex << 8) + minver[[i]]
if sys.hexversion >= minverhex:
    sys.exit( 0 )
else:
    sys.exit( 1 )
        ])
        if test $az_python_cc -eq 0
        then
            $2
        m4_ifvaln(
            [$3],
            [else $3]
        )
        fi
    fi
])



# AZ_PYTHON_VERSION_ENSURE( VERSION )
# -----------------
# Insure that the Python Interpreter Version
# is greater than or equal to the VERSION
# parameter.

AC_DEFUN([AZ_PYTHON_VERSION_ENSURE],
[
    AZ_PYTHON_VERSION_CHECK(
        [$1],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR(too old)]
    )
])



# AZ_PYTHON_WITH( [path] )
# -----------------------------------------------------------------
# Handles the various --with-python commands.
# Input:
#   $1 is the optional search path for the python executable if needed
# Ouput:
#   PYTHON_USE (AM_CONDITIONAL) is true if python executable found
#   and --with-python was requested; otherwise false.
#   $PYTHON contains the full executable path to python if PYTHON_USE
#   is true.
#
# Example:
#   AZ_PYTHON_WITH( )
#   or
#   AZ_PYTHON_WITH("/usr/bin")

AC_DEFUN([AZ_PYTHON_WITH],
[
    AC_ARG_VAR([PYTHON],[Python Executable Path])

    # unless PYTHON was supplied to us (as a precious variable),
    # see if --with-python[=PythonExecutablePath], --with-python,
    # --without-python or --with-python=no was given.
    if test -z "$PYTHON"
    then
        AC_MSG_CHECKING(for --with-python)
        AC_ARG_WITH(
            python,
            AC_HELP_STRING([--with-python@<:@=PYTHON@:>@],
                [absolute path name of Python executable]
            ),
            [
                if test "$withval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG(
                        [PYTHON],
                        python,
                        [],
                        $1
                    )
                    if test -z "$PYTHON"
                    then
                        AC_MSG_ERROR(no path to python found)
                    fi
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                elif test "$withval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_python_use=false
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                else
                    # $withval must be the executable path then.
                    AC_SUBST([PYTHON], ["${withval}"])
                    AC_MSG_RESULT($withval)
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                fi
            ],
            [
                # --with-python was not specified.
                AC_MSG_RESULT(no)
                az_python_use=false
                AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
            ]
        )
    fi

])

# ===========================================================================
#            http://autoconf-archive.cryp.to/ac_python_devel.html
# ===========================================================================
#
# SYNOPSIS
#
#   AC_PYTHON_DEVEL([version])
#
# DESCRIPTION
#
#   Note: Defines as a precious variable "PYTHON_VERSION". Don't override it
#   in your configure.ac.
#
#   This macro checks for Python and tries to get the include path to
#   'Python.h'. It provides the $(PYTHON_CPPFLAGS) and $(PYTHON_LDFLAGS)
#   output variables. It also exports $(PYTHON_EXTRA_LIBS) and
#   $(PYTHON_EXTRA_LDFLAGS) for embedding Python in your code.
#
#   You can search for some particular version of Python by passing a
#   parameter to this macro, for example ">= '2.3.1'", or "== '2.4'". Please
#   note that you *have* to pass also an operator along with the version to
#   match, and pay special attention to the single quotes surrounding the
#   version number. Don't use "PYTHON_VERSION" for this: that environment
#   variable is declared as precious and thus reserved for the end-user.
#
#   This macro should work for all versions of Python >= 2.1.0. As an end
#   user, you can disable the check for the python version by setting the
#   PYTHON_NOVERSIONCHECK environment variable to something else than the
#   empty string.
#
#   If you need to use this macro for an older Python version, please
#   contact the authors. We're always open for feedback.
#
# LAST MODIFICATION
#
#   2008-04-12
#
# COPYLEFT
#
#   Copyright (c) 2008 Sebastian Huber <sebastian-huber@web.de>
#   Copyright (c) 2008 Alan W. Irwin <irwin@beluga.phys.uvic.ca>
#   Copyright (c) 2008 Rafael Laboissiere <rafael@laboissiere.net>
#   Copyright (c) 2008 Andrew Collier <colliera@ukzn.ac.za>
#   Copyright (c) 2008 Matteo Settenvini <matteo@member.fsf.org>
#   Copyright (c) 2008 Horst Knorr <hk_classes@knoda.org>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AC_PYTHON_DEVEL],[
	#
	# Allow the use of a (user set) custom python version
	#
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string
		will be appended to the Python interpreter
		canonical name.])

	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test -z "$PYTHON"; then
	   AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	fi

	#
	# Check for a version of Python >= 2.1.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys, string; \
		ver = string.split(sys.version)[[0]]; \
		print ver >= '2.1.0'"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
		else
			AC_MSG_RESULT([skip at user request])
		fi
	else
		AC_MSG_RESULT([yes])
	fi

	#
	# if the macro parameter ``version'' is set, honour it
	#
	if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys, string; \
			ver = string.split(sys.version)[[0]]; \
			print ver $1"`
		if test "$ac_supports_python_ver" = "True"; then
	   	   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_ERROR([this package requires Python $1.
If you have it installed, but it isn't the default Python
interpreter in your system path, please pass the PYTHON_VERSION
variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	fi

	#
	# Check if you have distutils, else fail
	#
	AC_MSG_CHECKING([for the distutils Python package])
	ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
	if test -z "$ac_distutils_result"; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		AC_MSG_ERROR([cannot import Python module "distutils".
Please check your Python installation. The error was:
$ac_distutils_result])
		PYTHON_VERSION=""
	fi

	#
	# Check for Python include path
	#
	AC_MSG_CHECKING([for Python include path])
	if test -z "$PYTHON_CPPFLAGS"; then
		python_path=`$PYTHON -c "import distutils.sysconfig; \
           		print distutils.sysconfig.get_python_inc();"`
		if test -n "${python_path}"; then
		   	python_path="-I$python_path"
		fi
		PYTHON_CPPFLAGS=$python_path
	fi
	AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	AC_SUBST([PYTHON_CPPFLAGS])

	#
	# Check for Python library path
	#
	AC_MSG_CHECKING([for Python library path])
	if test -z "$PYTHON_LDFLAGS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		py_version=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print join(get_config_vars('VERSION'))"`
		if test "$py_version" == "[None]"; then
			if test -n "$PYTHON_VERSION"; then
				py_version=$PYTHON_VERSION
			else
				py_version=`$PYTHON -c "import sys; \
					print sys.version[[:3]]"`
			fi
		fi

		PYTHON_LDFLAGS=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print '-L' + get_python_lib(0,1), \
		      	'-lpython';"`$py_version
	fi
	AC_MSG_RESULT([$PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_LDFLAGS])

	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
		        print distutils.sysconfig.get_python_lib(0,0);"`
	fi
	AC_MSG_RESULT([$PYTHON_SITE_PKG])
	AC_SUBST([PYTHON_SITE_PKG])

	#
	# libraries which must be linked in when embedding
	#
	AC_MSG_CHECKING(python extra libraries)
	if test -z "$PYTHON_EXTRA_LIBS"; then
	   PYTHON_EXTRA_LIBS=`$PYTHON -c "import distutils.sysconfig; \
                conf = distutils.sysconfig.get_config_var; \
                print conf('LOCALMODLIBS'), conf('LIBS')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LIBS])
	AC_SUBST(PYTHON_EXTRA_LIBS)

	#
	# linking flags needed when embedding
	#
	AC_MSG_CHECKING(python extra linking flags)
	if test -z "$PYTHON_EXTRA_LDFLAGS"; then
		PYTHON_EXTRA_LDFLAGS=`$PYTHON -c "import distutils.sysconfig; \
			conf = distutils.sysconfig.get_config_var; \
			print conf('LINKFORSHARED')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	AC_LANG_PUSH([C])
	# save current global flags
	LIBS="$ac_save_LIBS $PYTHON_LDFLAGS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_TRY_LINK([
		#include <Python.h>
	],[
		Py_Initialize();
	],[pythonexists=yes],[pythonexists=no])

	AC_MSG_RESULT([$pythonexists])

        if test ! "$pythonexists" = "yes"; then
	   AC_MSG_ERROR([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
	  PYTHON_VERSION=""
	fi
	AC_LANG_POP
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"

	#
	# all done!
	#
])
