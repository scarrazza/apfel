#-------------------------------------------------
#
# Project created by QtCreator 2013-10-16T19:54:08
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = apfelgui
TEMPLATE = app

INCLUDEPATH += /opt/apfel/include
INCLUDEPATH += $$system(root-config --incdir)

LIBS += -L/opt/apfel/lib -lAPFEL
LIBS += $$system(lhapdf-config --ldflags)
LIBS += $$system(root-config --glibs)


SOURCES += main.cxx \
           apfelmainwindow.cxx \
    apfelthread.cxx

HEADERS  += apfelmainwindow.h \
    apfelthread.h

FORMS    += apfelmainwindow.ui

RESOURCES += \
    resource.qrc
