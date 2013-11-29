#include "plotdis.h"
#include "ui_plotdis.h"

#include <QGraphicsScene>
#include <QGraphicsItemGroup>
#include <QGraphicsSvgItem>
#include <QFile>
#include <QtSvg/QSvgWidget>
#include <QDebug>
#include <QFileDialog>
#include <QDesktopWidget>

PlotDIS::PlotDIS(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::PlotDIS)
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());
}

PlotDIS::~PlotDIS()
{
  delete ui;
}
