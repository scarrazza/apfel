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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
class apfelthread;
class QProgressDialog;

namespace Ui {
  class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();
  
private slots:
  void on_pushButton_clicked();
  void ThreadFinished();
  void ThreadProgress(int);
  void Thread2Progress();
  void Thread2Finished();
  void on_comboBox_3_currentIndexChanged(int index);
  void on_comboBox_currentIndexChanged(int index);
  void on_pushButton_2_clicked();
  void on_pushButton_3_clicked();
  void on_pushButton_4_clicked();

private:
  Ui::MainWindow *ui;
  apfelthread *thread;
  apfelthread *thread2;
  QProgressDialog *d;
};

#endif // MAINWINDOW_H
