#ifndef PLOTDIS_H
#define PLOTDIS_H

#include <QWidget>
#include <QThread>

class PDFDialog;
class TCanvas;
class PlotDIS;

namespace Ui {
  class PlotDIS;
}

class disthread: public QThread
{
  Q_OBJECT

public:
  disthread(QObject *parent,QString filename);
  ~disthread();
  void run();      //start calculation
  void stop();
  void SaveCanvas(QString filename);
  QString getFileName() { return fFileName; }

signals:
  void progress(int i);

private:
  TCanvas *fC;
  PlotDIS *fp;
  QString fFileName;
  bool fIsTerminated;

};


class PlotDIS : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotDIS(QWidget *parent = 0, PDFDialog *pdf = NULL);
  ~PlotDIS();
  
private slots:
  void on_process_currentIndexChanged(int index);
  void on_output_currentIndexChanged(int index);
  void on_automaticrange_toggled(bool checked);
  void on_playButton_clicked();
  void on_saveButton_clicked();
  void ThreadFinished();
  void ThreadProgress(int);

private:
  Ui::PlotDIS *ui;
  disthread *thread;
  PDFDialog *fPDF;
  QString fPlotName;
  bool fIsRunning;

  friend class disthread;
};

#endif // PLOTDIS_H
