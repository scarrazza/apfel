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
  QProgressDialog *d2;
};

#endif // MAINWINDOW_H
