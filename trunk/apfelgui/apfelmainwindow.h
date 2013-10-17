#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
class apfelthread;

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

private:
  Ui::MainWindow *ui;
  apfelthread *thread;
};

#endif // MAINWINDOW_H
