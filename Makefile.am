ACLOCAL_AMFLAGS = -I m4

SUBDIRS = include ccwrap src lib examples pywrap doc bin
dist_doc_DATA = README
dist_pkgdata_DATA = src/HELL/data/*

mostlyclean-local:
	rm -rf a.out
	rm -rf autom4te.cache
	rm -rf src/DIS/.deps src/core/.deps
	rm -rf src/DIS_new/.deps

