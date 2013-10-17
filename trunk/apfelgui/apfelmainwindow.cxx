#include <QMessageBox>
#include <QAction>
#include <QString>
#include <QGraphicsPixmapItem>
#include <QDesktopWidget>
#include <QProgressDialog>

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

  thread = new apfelthread(this,ui);
  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));

}

MainWindow::~MainWindow()
{
  delete ui;
  delete thread;
}

void MainWindow::on_pushButton_clicked()
{
  QProgressDialog progress("Copying files...", "Abort Copy",0,10, this);
  progress.setWindowModality(Qt::WindowModal);

  thread->start();
}

void MainWindow::ThreadFinished()
{
  QImage image("plot.png");

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);
  QGraphicsPixmapItem* item = new QGraphicsPixmapItem(QPixmap::fromImage(image));
  scene->addItem(item);
  ui->graphicsView->setScene(scene);
  ui->graphicsView->show();
}
