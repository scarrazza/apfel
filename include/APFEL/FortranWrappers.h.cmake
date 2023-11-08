#ifndef APFEL_FORTRANWRAPPERS_H
#define APFEL_FORTRANWRAPPERS_H

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) @FortranCInterface_GLOBAL_PREFIX@##@FortranCInterface_GLOBAL_SUFFIX@

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_


#endif
