ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src include ccwrap lib examples pywrap doc bin apfelGUI
dist_doc_DATA = README
dist_pkgdata_DATA = src/HELL/data/*

mostlyclean-local:
	rm -rf a.out
	rm -rf autom4te.cache
	rm -rf src/DIS/.deps src/core/.deps

