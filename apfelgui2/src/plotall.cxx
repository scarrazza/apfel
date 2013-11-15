#include "plotall.h"
#include "ui_plotall.h"
#include "pdfdialog.h"
#include "common.h"

#include <QDesktopWidget>

PlotAll::PlotAll(QWidget *parent, PDFDialog *pdf) :
  QWidget(parent),
  ui(new Ui::PlotAll),
  fPDF(pdf),
  fPlotName("allplot.svg")
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  long int t = static_cast<long int> (time(NULL));
  QString str;
  str.append(QString("%1").arg(t));

  fPlotName = "memberplot_" + str + ".svg";
  //thread  = new memberthread(this,fPlotName);

  //connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  //connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));
  /*
  ui->graphicsView->scale(1.2,1.2);

  if (fPDF->isLHAPDF()) ui->Qi->setEnabled(false);

  ui->xtitle->setText("x");
  ui->ytitle->setText("");
  ui->title->setText(name[ui->PDFflavor->currentIndex()]+ ", " + fPDF->PDFname()+ " members");
  */
}

PlotAll::~PlotAll()
{
  delete ui;
}
