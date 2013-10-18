#-------------------------------------------------
#
# Project created by QtCreator 2013-10-16T19:54:08
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = apfelgui
TEMPLATE = app


INCLUDEPATH += $$system(root-config --incdir)
INCLUDEPATH += $$system(apfel-config --incdir)
INCLUDEPATH += include

LIBS += $$system(apfel-config --ldflags)
LIBS += $$system(lhapdf-config --ldflags)
LIBS += $$system(root-config --glibs)

DEPENDPATH += . forms resources src include

CONFIG += release
OBJECTS_DIR = src
UI_DIR = include/
MOC_DIR = src/
RCC_DIR = resources/
DISTFILES += apfelgui

SOURCES += src/main.cxx \
           src/apfelmainwindow.cxx \
           src/apfelthread.cxx

HEADERS  += include/apfelmainwindow.h \
            include/apfelthread.h

FORMS    += forms/apfelmainwindow.ui

RESOURCES += resources/resource.qrc
