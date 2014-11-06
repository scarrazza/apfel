#ifndef PLOTALL_H
#define PLOTALL_H

#include <QWidget>
#include <QThread>
class PDFDialog;
class TCanvas;
class PlotAll;

namespace Ui {
  class PlotAll;
}


class plotthread: public QThread
{
  Q_OBJECT

public:
  plotthread(QObject *parent,QString filename);
  ~plotthread();
  void run();      //start calculation
  void stop();
  void SaveCanvas(QString filename);
  QString getFileName() { return fFileName; }

signals:
  void progress(int i);

private:
  TCanvas *fC;
  PlotAll *fp;
  QString fFileName;
  bool fIsTerminated;

};

class PlotAll : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotAll(QWidget *parent = 0, PDFDialog *pdf = NULL);
  ~PlotAll();
  
private slots:
  void on_playButton_clicked();
  void on_checkBox_toggled(bool checked);
  void on_automaticrange_toggled(bool checked);
  void ThreadFinished();
  void ThreadProgress(int);
  void on_saveButton_clicked();

private:
  Ui::PlotAll *ui;
  plotthread *thread;
  PDFDialog *fPDF;
  QString fPlotName;
  bool fIsRunning;

  friend class plotthread;

};

#endif // PLOTALL_H
