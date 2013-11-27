#/****************************************************************************
#*   APFEL GUI
#*   Copyright (C) 2013- Stefano Carrazza. All rights reserved.
#*									    *
#*   This program is free software: you can redistribute it and/or modify
#*   it under the terms of the GNU General Public License as published by
#*   the Free Software Foundation, either version 3 of the License, or
#*   (at your option) any later version.
#*
#*   This program is distributed in the hope that it will be useful,
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*   GNU General Public License for more details.
#*
#*   You should have received a copy of the GNU General Public License
#*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************/

QMAKE_MACOSX_DEPLOYMENT_TARGET=10.7

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
