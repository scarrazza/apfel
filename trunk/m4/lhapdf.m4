# AC_SEARCH_LHAPDF(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_LHAPDF], [

AC_PATH_PROG(LHAPDFCONFIG, lhapdf-config, [], [$PATH])
if test -f "$LHAPDFCONFIG"; then
  LHAPDF_CPPFLAGS=`$LHAPDFCONFIG --cxxflags`
  LHAPDF_LDFLAGS=`$LHAPDFCONFIG --ldflags`
else
  AC_MSG_ERROR([LHAPDF cannot be found!])
  exit 1
fi
AC_SUBST(LHAPDF_CPPFLAGS)
AC_SUBST(LHAPDF_LDFLAGS)
$1
])
