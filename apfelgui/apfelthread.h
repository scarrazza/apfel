#ifndef APFELTHREAD_H
#define APFELTHREAD_H

#include <QThread>

namespace Ui {
  class MainWindow;
}

class QProgressDialog;

class apfelthread : public QThread
{
  Q_OBJECT

public:
  apfelthread(QObject *parent = 0, Ui::MainWindow *fui = 0, QProgressDialog *fd = 0);
  ~apfelthread();
  void run();      //start calculation

private:
  Ui::MainWindow *ui;
  QProgressDialog *d;
};

#endif // APFELTHREAD_H
