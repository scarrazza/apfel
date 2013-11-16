#ifndef PLOTMEMBERS_H
#define PLOTMEMBERS_H

#include <QWidget>
#include <QThread>

class PDFDialog;
class TCanvas;
class PlotMembers;

namespace Ui {
  class PlotMembers;
}

class memberthread: public QThread
{
  Q_OBJECT

public:
  memberthread(QObject *parent,QString filename);
  ~memberthread();
  void run();      //start calculation
  void SaveCanvas(QString filename);
  QString getFileName() { return fFileName; }

signals:
  void progress(int i);

private:
  TCanvas *fC;
  PlotMembers *fp;
  QString fFileName;
};


class PlotMembers : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotMembers(QWidget *parent = 0, PDFDialog *pdf = NULL);
  ~PlotMembers();
  
private slots:
  void on_playButton_clicked();
  void on_checkBox_toggled(bool checked);
  void ThreadFinished();
  void ThreadProgress(int);
  void on_automaticrange_toggled(bool checked);
  void on_saveButton_clicked();
  void on_PDFflavor_currentIndexChanged(int index);

private:
  Ui::PlotMembers *ui;
  memberthread *thread;
  PDFDialog *fPDF;
  QString fPlotName;

  friend class memberthread;
};

#endif // PLOTMEMBERS_H
