#include "plotall.h"
#include "ui_plotall.h"
#include "pdfdialog.h"
#include "common.h"
#include <cmath>

#include <QDesktopWidget>
#include <QGraphicsScene>
#include <QGraphicsItemGroup>
#include <QGraphicsSvgItem>
#include <QFile>
#include <QtSvg/QSvgWidget>
#include <QDebug>
#include <QFileDialog>

#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"

PlotAll::PlotAll(QWidget *parent, PDFDialog *pdf) :
  QWidget(parent),
  ui(new Ui::PlotAll),
  fPDF(pdf),
  fPlotName("allplot.svg"),
  fIsRunning(false)
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  long int t = static_cast<long int> (time(NULL));
  QString str;
  str.append(QString("%1").arg(t));

  fPlotName = "allplot_" + str + ".svg";
  thread  = new plotthread(this,fPlotName);

  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));

  ui->graphicsView->scale(1.2,1.2);

  if (fPDF->isLHAPDF()) ui->Qi->setEnabled(false);
  if (fPDF->numberPDF() == 1)
    {
      ui->checkBox->setEnabled(false);
      ui->stddev->setEnabled(false);
      ui->stddev->setChecked(false);
    }

  ui->xtitle->setText("x");
  ui->ytitle->setText("xf(x,Q)");
  ui->title->setText(fPDF->PDFname()+ " PDFs");

  if (fPDF->PDFname() == "ToyLH (APFEL internal)") {
      ui->stddev->setEnabled(false);
      ui->checkBox->setEnabled(false);
    }
}

PlotAll::~PlotAll()
{
  delete ui;
  if (thread)
    {
      thread->terminate();
      delete thread;
    }
}

void PlotAll::on_playButton_clicked()
{
  if (!fIsRunning)
    {
      fIsRunning = true;
      ui->playButton->setIcon(QIcon(":/images/Stop1NormalRed.png"));
      ui->playButton->setText("Stop");
      if (ui->graphicsView->scene()) ui->graphicsView->scene()->clear();
      QApplication::processEvents();
      thread->start();
    }
  else
    {
      fIsRunning = false;
      thread->stop();
      ui->playButton->setEnabled(false);
    }
}

void PlotAll::on_checkBox_toggled(bool checked)
{
  if (checked == true)
    {
      ui->setmember->setEnabled(false);
      ui->stddev->setEnabled(true);
    }
  else
    {
      ui->setmember->setEnabled(true);
      ui->stddev->setEnabled(false);
      ui->stddev->setChecked(false);
    }
}

void PlotAll::on_automaticrange_toggled(bool checked)
{
  if (checked == true)
    {
      ui->xmin->setEnabled(false);
      ui->xmax->setEnabled(false);
      ui->ymin->setEnabled(false);
      ui->ymax->setEnabled(false);
    }
  else
    {
      ui->xmin->setEnabled(true);
      ui->xmax->setEnabled(true);
      ui->ymin->setEnabled(true);
      ui->ymax->setEnabled(true);
    }
}

void PlotAll::ThreadFinished()
{
  ui->playButton->setIcon(QIcon(":/images/StepForwardNormalBlue.png"));
  ui->playButton->setText("Compute");
  fIsRunning = false;

  ui->playButton->setEnabled(true);
  ui->saveButton->setEnabled(true);

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);
  QGraphicsSvgItem * item = new QGraphicsSvgItem(fPlotName);
  scene->addItem(item);
  ui->graphicsView->setScene(scene);
  ui->graphicsView->show();

  ui->progressBar->setValue(0);

  QFile::remove(fPlotName) ;
}

void PlotAll::ThreadProgress(int i)
{
  ui->progressBar->setValue(i);
}

void PlotAll::on_saveButton_clicked()
{
  QString selectedFilter;
  QString path = QFileDialog::getSaveFileName(this,
                                              tr("Save as"),"",
                                              tr(".eps;;.ps;;.pdf;;.png;;.root;;.C"),&selectedFilter);
  if (path != 0) thread->SaveCanvas(path + selectedFilter);
}

plotthread::plotthread(QObject *parent, QString filename):
  QThread(parent),
  fp((PlotAll*)parent),
  fFileName(filename),
  fIsTerminated(false)

{
}

plotthread::~plotthread()
{
}

void plotthread::run()
{
  fC = new TCanvas();
  fC->SetTickx();
  fC->SetTicky();

  if (fp->ui->logx->isChecked())
    fC->SetLogx();

  if (fp->ui->logy->isChecked())
    fC->SetLogy();

  // Initialize PDFs

  const int N = fp->ui->xpoints->value();
  double xmin = 1e-3;
  double xmax = 1;
  double ymin = 0;
  double ymax = 1;;
  int legindex = 0;

  if (fp->ui->logy->isChecked()) ymin = 1e-5;

  if (!fp->ui->automaticrange->isChecked())
    {
      xmin = fp->ui->xmin->text().toDouble();
      xmax = fp->ui->xmax->text().toDouble();
      ymin = fp->ui->ymin->text().toDouble();
      ymax = fp->ui->ymax->text().toDouble();
    }

  const double Qi = fp->ui->Qi->text().toDouble();
  const double Qf = fp->ui->Qf->text().toDouble();
  fp->fPDF->InitPDFset(Qi,Qf);

  TLegend *leg = NULL;
  if (fp->ui->logx->isChecked())
    leg = new TLegend(0.130747,0.491525,0.300287,0.883475);
  else
    leg = new TLegend(0.716954,0.491525,0.886494,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TMultiGraph *mg = new TMultiGraph();

  int nf = fp->ui->maxFlavors->value();
  double *x = new double[N];
  for (int i = 0; i < N; i++)
    {
      if (fp->ui->logx->isChecked()) x[i] = exp(log(xmin)+i*(log(xmax)-log(xmin))/(N-1));
      else x[i] = xmin+i*(xmax-xmin)/(N-1);
    }

  int scales [] = { fp->ui->gluon->value(),
                    fp->ui->down->value(),fp->ui->up->value(),fp->ui->strange->value(),
                    fp->ui->charm->value(),fp->ui->bottom->value(),fp->ui->top->value()
                  };

  for (int fl = -nf; fl <= nf; fl++)
    {
      if(fIsTerminated){ fIsTerminated = false; return; }

      emit progress( (fl+nf)*100./(2*nf+1));

      TGraphErrors *g = new TGraphErrors(N);
      g->SetLineWidth(2);
      g->SetLineColor(colors2[fl+6]);
      g->SetFillColor(colors2[fl+6]);      

      double *xPDF = new double[N];
      double *xPDFErr = new double[N];
      double *upErr = new double[N];
      double *dnErr = new double[N];

      if (!fp->ui->checkBox->isChecked())
        {
          fp->fPDF->initPDF(fp->ui->setmember->value());
          for (int i = 0; i < N; i++) {
            xPDF[i] = fp->fPDF->GetFlvrPDF(x[i],Qf,fl);
            xPDFErr[i] = 0;
          }
        }
      else
        fp->fPDF->GetFlvrPDFCVErr(N,x,Qf,fl,xPDF,xPDFErr,upErr,dnErr);

      for (int i = 0; i < N; i++)
        {
          xPDF[i] /= scales[abs(fl)];
          g->SetPoint(i,x[i],xPDF[i]);

          if (fp->ui->stddev->isChecked())
            g->SetPointError(i,0,xPDFErr[i]);
        }

      delete[] xPDF;
      delete[] xPDFErr;
      delete[] upErr;
      delete[] dnErr;

      if (scales[abs(fl)] != 1)
        leg->AddEntry(g,TString(name[fl+6].toStdString() + "/" + Form("%d",scales[abs(fl)])),"l");
      else
        leg->AddEntry(g,name[fl+6].toStdString().c_str(),"l");
      legindex++;
      mg->Add(g,"le3");
    }

  delete[] x;

  mg->SetTitle(fp->ui->title->text().toStdString().c_str());
  mg->Draw("AL");

  mg->GetXaxis()->SetTitle(fp->ui->xtitle->text().toStdString().c_str());
  mg->GetXaxis()->CenterTitle(true);

  mg->GetYaxis()->SetTitle(fp->ui->ytitle->text().toStdString().c_str());
  mg->GetYaxis()->CenterTitle(true);

  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetLabelSize(0.05);

  mg->GetXaxis()->SetLimits(xmin,xmax);
  if (!fp->ui->automaticrange->isChecked())
    {
      mg->GetYaxis()->SetRangeUser(ymin,ymax);
    }

  leg->AddEntry("",TString("Q = " + QString::number(fp->ui->Qf->text().toDouble(),'g',3).toStdString() + " GeV"),"");
  if (legindex < 13)
    for (int i = legindex; i < 13; i++)
      leg->AddEntry(""," ","");
  leg->Draw("same");

  TLatex l; //l.SetTextAlign(12);
  l.SetTextSize(0.02);
  l.SetTextAngle(90);
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(kBlack);
  l.DrawLatex(0.95,0.15,TString("Generated by APFEL" + TString(APFEL::GetVersion()) + ": V.Bertone, S.Carrazza, J.Rojo (arXiv:1310.1394)"));

  fC->SaveAs(fFileName.toStdString().c_str());
}

void plotthread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}

void plotthread::stop()
{
  fIsTerminated = true;
}
