/****************************************************************************
*   APFEL GUI                                                               *
*   Copyright (C) 2013- Stefano Carrazza. All rights reserved.              *
*									    *
*   This program is free software: you can redistribute it and/or modify    *
*   it under the terms of the GNU General Public License as published by    *
*   the Free Software Foundation, either version 3 of the License, or	    *
*   (at your option) any later version.				            *
*									    *
*   This program is distributed in the hope that it will be useful,	    *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of	    *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	    *
*   GNU General Public License for more details.			    *
*									    *
*   You should have received a copy of the GNU General Public License	    *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/

#ifndef APFELTHREAD_H
#define APFELTHREAD_H

#include <string>
#include <QThread>
using std::string;

namespace Ui {
  class MainWindow;
}

class QProgressDialog;
class TCanvas;

class apfelthread : public QThread
{
  Q_OBJECT

public:
  apfelthread(QObject *parent = 0, Ui::MainWindow *fui = 0, int fmod = 0);
  ~apfelthread();
  void run();      //start calculation
  void savecanvas(string file);
  QString getresult(int i);

signals:
  void progress(int i);
  void progress();

private:
  Ui::MainWindow *ui;
  TCanvas *C;
  int mod;
  double momsr;
  double uvsr;
  double dvsr;
  double svsr;
};

#endif // APFELTHREAD_H
