#ifndef APFEL_FORTRANWRAPPERS_H
#define APFEL_FORTRANWRAPPERS_H

#define STDCALLBULL

/* #include "FCMangle.h" */

#cmakedefine HAVE_INTTYPES_H
/* #cmakedefine FC_FUNC ${FC_FUNC} */
/* #cmakedefine FC_FUNC_ ${FC_FUNC_} */
/* #cmakedefine FC_CHAR_PTR${FC_CHAR_PTR} */

#define FC_FUNC(name, NAME) name ## _
#define FC_FUNC_(name, NAME) name ## _

#endif
