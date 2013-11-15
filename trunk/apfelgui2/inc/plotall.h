#ifndef PLOTALL_H
#define PLOTALL_H

#include <QWidget>
class PDFDialog;

namespace Ui {
  class PlotAll;
}

class PlotAll : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotAll(QWidget *parent = 0, PDFDialog *pdf = NULL);
  ~PlotAll();
  
private:
  Ui::PlotAll *ui;
  PDFDialog *fPDF;
  QString fPlotName;
};

#endif // PLOTALL_H
