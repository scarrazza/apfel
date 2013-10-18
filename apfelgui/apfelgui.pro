#-------------------------------------------------
#
# Project created by QtCreator 2013-10-16T19:54:08
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = apfelgui
TEMPLATE = app

# Add the following lines if apfel is installed in a custom location
#INCLUDEPATH += /opt/apfel/include
#LIBS += -L/opt/apfel/lib -lAPFEL

INCLUDEPATH += $$system(root-config --incdir)
INCLUDEPATH += include

LIBS += -lAPFEL
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
