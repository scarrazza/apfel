#ifndef PLOTDISQ_H
#define PLOTDISQ_H

#include <QWidget>
#include <QThread>

class PDFDialog;
class TCanvas;
class PlotDISQ;

namespace Ui {
  class PlotDISQ;
}

class disqthread: public QThread
{
  Q_OBJECT

public:
  disqthread(QObject *parent,QString filename);
  ~disqthread();
  void run();      //start calculation
  void stop();
  void SaveCanvas(QString filename);
  QString getFileName() { return fFileName; }

signals:
  void progress(int i);

private:
  TCanvas *fC;
  PlotDISQ *fp;
  QString fFileName;
  bool fIsTerminated;

};

class PlotDISQ : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotDISQ(QWidget *parent = 0, PDFDialog *pdf = NULL);
  ~PlotDISQ();
  
private slots:
  void on_process_currentIndexChanged(int index);
  void on_output_currentIndexChanged(int index);
  void on_automaticrange_toggled(bool checked);
  void on_playButton_clicked();
  void on_saveButton_clicked();
  void ThreadFinished();
  void ThreadProgress(int);

private:
  Ui::PlotDISQ *ui;
  disqthread *thread;
  PDFDialog *fPDF;
  QString fPlotName;
  bool fIsRunning;

  friend class disqthread;

};

#endif // PLOTDISQ_H
