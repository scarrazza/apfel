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
#include <QStatusBar>
#include <QDebug>

#include "apfelmainwindow.h"
#include "ui_apfelmainwindow.h"
#include "apfelthread.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow),
  thread(0),
  thread2(0),
  d(0),
  d2(0)
{
  ui->setupUi(this);

  statusBar()->showMessage("APFEL: V. Bertone, S. Carrazza and J. Rojo (arXiv:1310.1394)");

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
  d->setWindowTitle("Please wait");
  d->setModal(true);
  d->setCancelButton(0);
  thread  = new apfelthread(this,ui,d,0);

  d2 = new QProgressDialog("Computing, please wait...", "", 0, 500, this);
  d2->setWindowTitle("Please wait");
  d2->setModal(true);
  d2->setCancelButton(0);
  thread2 = new apfelthread(this,ui,d2,1);
  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread2, SIGNAL(finished()), this, SLOT(Thread2Finished()));

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
  ui->pushButton_2->setEnabled(true);

  QImage image("apfelplot.svg");

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);
  QGraphicsPixmapItem* item = new QGraphicsPixmapItem(QPixmap::fromImage(image));
  scene->addItem(item);
  ui->graphicsView->setScene(scene);
  ui->graphicsView->show();

  d->setValue(10);

  QFile::remove("apfelplot.svg") ;
}

void MainWindow::Thread2Finished()
{
  ui->pushButton_3->setEnabled(true);
  ui->pushButton_4->setEnabled(true);

  QImage image("apfelplot2.svg");

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView_2);
  QGraphicsPixmapItem* item = new QGraphicsPixmapItem(QPixmap::fromImage(image));
  scene->addItem(item);
  ui->graphicsView_2->setScene(scene);
  ui->graphicsView_2->show();  

  d2->setValue(4 + 10*ui->spinBox_5->value());

  QFile::remove("apfelplot2.svg") ;
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


void MainWindow::on_pushButton_2_clicked()
{
  QString path;
  path = QFileDialog::getSaveFileName(this,tr("Save as"),"",tr("*.eps (*.eps);;All files (*.*)"));

  if(path != 0) thread->savecanvas(path.toStdString() + ".eps");
}

void MainWindow::on_pushButton_3_clicked()
{
  ui->pushButton_3->setEnabled(false);

  d2->setMaximum(4 + 10*ui->spinBox_5->value());

  d2->setValue(1);
  d2->show();
  QApplication::processEvents();

  thread2->start();
}

void MainWindow::on_pushButton_4_clicked()
{
  QString path;
  path = QFileDialog::getSaveFileName(this,tr("Save as"),"",tr("*.eps (*.eps);;All files (*.*)"));

  if(path != 0) thread2->savecanvas(path.toStdString() + ".eps");
}
