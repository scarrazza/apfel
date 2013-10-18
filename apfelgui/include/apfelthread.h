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
  apfelthread(QObject *parent = 0, Ui::MainWindow *fui = 0, QProgressDialog *fd = 0, int fmod = 0);
  ~apfelthread();
  void run();      //start calculation
  void savecanvas(string file);

private:
  Ui::MainWindow *ui;
  QProgressDialog *d;
  TCanvas *C;
  int mod;
};

#endif // APFELTHREAD_H
