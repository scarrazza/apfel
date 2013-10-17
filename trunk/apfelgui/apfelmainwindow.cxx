#include <QWidget>
#include <QMessageBox>
#include <QAction>
#include <QString>
#include <QGraphicsPixmapItem>
#include <QDesktopWidget>
#include <QProgressDialog>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QComboBox>
#include <QPushButton>
#include <QFileDialog>

#include "apfelmainwindow.h"
#include "ui_apfelmainwindow.h"
#include "apfelthread.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow),
  thread(0)
{
  ui->setupUi(this);

  //center mainwindow position on desktop
  QDesktopWidget *desktop = QApplication::desktop();

  int screenWidth, width;
  int screenHeight, height;
  int x, y;
  QSize windowSize;

  screenWidth = desktop->width(); // get width of screen
  screenHeight = desktop->height(); // get height of screen

  windowSize = size(); // size of our application window
  width = windowSize.width();
  height = windowSize.height();

  // little computations
  x = (screenWidth - width) / 2;
  y = (screenHeight - height) / 2;

  // move window to desired coordinates
  move ( x, y );

  d = new QProgressDialog("Computing, please wait...", "", 0, 10, this);
  d->setModal(true);
  d->setCancelButton(0);
  thread = new apfelthread(this,ui,d);
  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));

}

MainWindow::~MainWindow()
{
  delete ui;
  delete thread;
}

void MainWindow::on_pushButton_clicked()
{
  ui->pushButton->setEnabled(false);

  d->setValue(1);
  d->show();
  QApplication::processEvents();

  thread->start();
}

void MainWindow::ThreadFinished()
{
  ui->pushButton->setEnabled(true);

  QImage image("plot.png");

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);
  QGraphicsPixmapItem* item = new QGraphicsPixmapItem(QPixmap::fromImage(image));
  scene->addItem(item);
  ui->graphicsView->setScene(scene);
  ui->graphicsView->show();

  d->setValue(10);
}

void MainWindow::on_comboBox_3_currentIndexChanged(int index)
{
  if (index == 1)
    {
      ui->spinBox->setEnabled(true);
      ui->spinBox_2->setEnabled(false);
      ui->spinBox_3->setEnabled(false);
      ui->comboBox_4->setEnabled(false);
      ui->lineEdit_10->setEnabled(false);
      ui->lineEdit_11->setEnabled(false);
      ui->lineEdit_12->setEnabled(false);
    }
  else
    {
      ui->spinBox->setEnabled(false);
      ui->spinBox_2->setEnabled(true);
      ui->spinBox_3->setEnabled(true);
      ui->comboBox_4->setEnabled(true);
      ui->lineEdit_10->setEnabled(true);
      ui->lineEdit_11->setEnabled(true);
      ui->lineEdit_12->setEnabled(true);
    }
}

void MainWindow::on_comboBox_currentIndexChanged(int index)
{
  if (index == 0)
    {
      ui->lineEdit_13->setEnabled(true);
      ui->lineEdit_14->setEnabled(true);
      ui->lineEdit_15->setEnabled(false);
      ui->lineEdit_16->setEnabled(false);

      ui->comboBox_2->setEnabled(true);
    }
  else if (index == 1)
    {
      ui->lineEdit_13->setEnabled(false);
      ui->lineEdit_14->setEnabled(false);
      ui->lineEdit_15->setEnabled(true);
      ui->lineEdit_16->setEnabled(true);

      ui->comboBox_2->setEnabled(false);
    }
  else
    {
      ui->lineEdit_13->setEnabled(true);
      ui->lineEdit_14->setEnabled(true);
      ui->lineEdit_15->setEnabled(true);
      ui->lineEdit_16->setEnabled(true);

      ui->comboBox_2->setEnabled(true);
    }

}
