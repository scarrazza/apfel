#-------------------------------------------------
#
# Project created by QtCreator 2013-10-25T23:16:49
#
#-------------------------------------------------

QT       += core gui svg

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = apfelgui
TEMPLATE = app

INCLUDEPATH += $$system(root-config --incdir)
INCLUDEPATH += $$system(apfel-config --incdir)
INCLUDEPATH += $$system(lhapdf-config --ldflags)
INCLUDEPATH += inc

INCLUDEPATH += /opt/apfel/include
INCLUDEPATH += /opt/lhapdf/include
#LIBS += -L/opt/apfel/lib -lAPFEL -L/opt/lhapdf/lib -lLHAPDF

DEPENDPATH += . forms src inc

LIBS += $$system(apfel-config --ldflags)
LIBS += $$system(lhapdf-config --ldflags)
LIBS += $$system(root-config --glibs)

OBJECTS_DIR = src
UI_DIR = inc/
MOC_DIR = src/
RCC_DIR = resources/
DISTFILES += apfelgui

SOURCES += src/main.cxx\
           src/apfelmainwindow.cxx \
           src/pdfdialog.cxx \
           src/plotmembers.cxx \
           src/plotall.cxx \
           src/common.cxx \
           src/plotcomparison.cxx \
           src/plotlumi.cxx \
           src/utils.cxx

HEADERS  += inc/apfelmainwindow.h \
            inc/pdfdialog.h \
            inc/plotmembers.h \
            inc/plotall.h \
            inc/common.h \
            inc/plotcomparison.h \
            inc/plotlumi.h \
            inc/utils.h

FORMS    += forms/apfelmainwindow.ui \
            forms/pdfdialog.ui \
            forms/plotmembers.ui \
            forms/plotall.ui \
            forms/plotcomparison.ui \
            forms/plotlumi.ui

RESOURCES += \
    resources/resource.qrc
