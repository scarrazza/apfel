#include <QMessageBox>
#include <QAction>
#include <QString>
#include <QGraphicsPixmapItem>
#include <QDesktopWidget>
#include <QApplication>
#include <QLineEdit>

#include "apfelthread.h"
#include "apfelmainwindow.h"
#include "ui_apfelmainwindow.h"

#include "APFEL/APFEL.h"
#include "TCanvas.h"
#include "TF1.h"

apfelthread::apfelthread(QObject *parent, Ui::MainWindow *fui):
  QThread(parent),
  ui(fui)
{
}

apfelthread::~apfelthread()
{

}

void apfelthread::run()
{
  // Initialize apfel
  APFEL::InitializeAPFEL();
  APFEL::EvolveAPFEL(ui->doubleSpinBox->value(), ui->doubleSpinBox_2->value());

  ui->lineEdit_3->setText(QString::number(APFEL::AlphaQCD(ui->doubleSpinBox_2->value())));
  ui->lineEdit_4->setText(QString::number(APFEL::AlphaQED(ui->doubleSpinBox_2->value())));

  // Plot
  TCanvas *c = new TCanvas("c","",466.66,400);
  TF1 *f = new TF1("f","sin(x)",-5,5);
  f->Draw();
  c->SaveAs("plot.png");

}
