#ifndef PLOTLUMI_H
#define PLOTLUMI_H

#include <QWidget>
#include <QThread>

class PDFDialog;
class PlotLumi;
class TCanvas;

namespace Ui {
  class PlotLumi;
}

class lumithread: public QThread
{
  Q_OBJECT

public:
  lumithread(QObject *parent,QString filename);
  ~lumithread();
  void run();      //start calculation
  void SaveCanvas(QString filename);
  QString getFileName() { return fFileName; }

signals:
  void progress(int i);

private:
  TCanvas *fC;
  PlotLumi *fp;
  QString fFileName;
};

class PlotLumi : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotLumi(QWidget *parent = 0,std::vector<PDFDialog*> pdf = std::vector<PDFDialog*>());
  ~PlotLumi();
  
private slots:
  void on_referenceSet_currentIndexChanged(int index);
  void on_automaticrange_toggled(bool checked);
  void on_playButton_clicked();
  void on_saveButton_clicked();
  void ThreadFinished();
  void ThreadProgress(int);
  void on_PDFflavor_currentIndexChanged(int index);

private:
  Ui::PlotLumi *ui;
  lumithread *thread;
  std::vector<PDFDialog*> fPDF;
  QString fPlotName;
  int fRef;

  friend class lumithread;

};

#endif // PLOTLUMI_H
