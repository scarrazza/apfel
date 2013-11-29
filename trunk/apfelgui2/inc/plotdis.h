#ifndef PLOTDIS_H
#define PLOTDIS_H

#include <QWidget>

namespace Ui {
  class PlotDIS;
}

class PlotDIS : public QWidget
{
  Q_OBJECT
  
public:
  explicit PlotDIS(QWidget *parent = 0);
  ~PlotDIS();
  
private:
  Ui::PlotDIS *ui;
};

#endif // PLOTDIS_H
