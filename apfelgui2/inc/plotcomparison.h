#ifndef PLOTCOMPARISON_H
#define PLOTCOMPARISON_H

#include <QWidget>
#include <QThread>
class PDFDialog;
class PlotComparison;
class TCanvas;

namespace Ui {
  class PlotComparison;
}

class comparisonthread: public QThread
{
  Q_OBJECT

public:
  comparisonthread(QObject *parent,QString filename);
  ~comparisonthread();
  void run();      //start calculation
  void SaveCanvas(QString filename);
  QString getFileName() { return fFileName; }

signals:
  void progress(int i);

private:
  TCanvas *fC;
  PlotComparison *fp;
  QString fFileName;
};

class PlotComparison : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotComparison(QWidget *parent = 0,std::vector<PDFDialog*> pdf = std::vector<PDFDialog*>());
  ~PlotComparison();
  
private slots:
  void on_playButton_clicked();
  void on_saveButton_clicked();
  void ThreadFinished();
  void ThreadProgress(int);
  void on_PDFflavor_currentIndexChanged(int index);
  void on_referenceSet_currentIndexChanged(int index);
  void on_ratio_toggled(bool checked);
  void on_automaticrange_toggled(bool checked);

private:
  Ui::PlotComparison *ui;
  comparisonthread *thread;
  std::vector<PDFDialog*> fPDF;
  QString fPlotName;
  int fRef;

  friend class comparisonthread;
};

#endif // PLOTCOMPARISON_H
