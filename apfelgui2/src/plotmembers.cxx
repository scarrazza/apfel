#include "plotmembers.h"
#include "ui_plotmembers.h"
#include "pdfdialog.h"
#include <algorithm>
#include <vector>
#include <ctime>
#include <cmath>
#include "common.h"
#include "utils.h"
using namespace std;

#include <QGraphicsScene>
#include <QGraphicsItemGroup>
#include <QGraphicsSvgItem>
#include <QFile>
#include <QtSvg/QSvgWidget>
#include <QDebug>
#include <QFileDialog>
#include <QDesktopWidget>

#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"

PlotMembers::PlotMembers(QWidget *parent, PDFDialog *pdf) :
  QWidget(parent),
  ui(new Ui::PlotMembers),
  thread(NULL),
  fPDF(pdf),
  fPlotName("memberplot.svg"),
  fIsRunning(false)
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  long int t = static_cast<long int> (time(NULL));
  QString str;
  str.append(QString("%1").arg(t));

  fPlotName = "memberplot_" + str + ".svg";
  thread  = new memberthread(this,fPlotName);

  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));

  ui->graphicsView->scale(1.2,1.2);

  if (fPDF->isLHAPDF()) ui->Qi->setEnabled(false);
  if (fPDF->numberPDF() == 1)
    {
      ui->checkBox->setEnabled(false);
      ui->mean->setEnabled(false);
      ui->mean->setChecked(false);
      ui->stddev->setEnabled(false);
      ui->stddev->setChecked(false);
      ui->cl->setEnabled(false);
      ui->cl->setChecked(false);
    }

  if (fPDF->GetErrorType() != ER_MC) {
    ui->cl->setEnabled(false);
    ui->cl->setChecked(false);
  }

  ui->xtitle->setText("x");
  ui->ytitle->setText("");
  ui->title->setText(name[ui->PDFflavor->currentIndex()]+ ", " + fPDF->PDFname()+ " members");
}

PlotMembers::~PlotMembers()
{
  delete ui;
  if (thread)
    {
      thread->terminate();
      delete thread;
    }
}

void PlotMembers::on_playButton_clicked()
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

void PlotMembers::on_checkBox_toggled(bool checked)
{
  if (checked == true)
    {
      ui->setmember->setEnabled(false);
      ui->mean->setEnabled(true);
      ui->stddev->setEnabled(true);
      ui->cl->setEnabled(true);
    }
  else
    {
      ui->setmember->setEnabled(true);
      ui->mean->setEnabled(false);
      ui->mean->setChecked(false);
      ui->stddev->setEnabled(false);
      ui->stddev->setChecked(false);
      ui->cl->setEnabled(false);
      ui->cl->setChecked(false);
    }
}

void PlotMembers::on_automaticrange_toggled(bool checked)
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

void PlotMembers::ThreadFinished()
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

void PlotMembers::ThreadProgress(int i)
{
  ui->progressBar->setValue(i);
}

void PlotMembers::on_saveButton_clicked()
{
  QString selectedFilter;
  QString path = QFileDialog::getSaveFileName(this,
                                              tr("Save as"),"",
                                              tr(".eps;;.ps;;.pdf;;.png;;.root;;.C"),&selectedFilter);
  if (path != 0) thread->SaveCanvas(path + selectedFilter);
}

void PlotMembers::on_PDFflavor_currentIndexChanged(int index)
{
  ui->title->setText(name[index]+ ", " + fPDF->PDFname()+ " members");
}

memberthread::memberthread(QObject *parent, QString filename):
  QThread(parent),
  fp((PlotMembers*)parent),
  fFileName(filename),
  fIsTerminated(false)
{
}

memberthread::~memberthread()
{
}

void memberthread::run()
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
  double xmin = lharanges[fp->ui->PDFflavor->currentIndex()][0];
  double xmax = lharanges[fp->ui->PDFflavor->currentIndex()][1];
  double ymin = lharanges[fp->ui->PDFflavor->currentIndex()][2];
  double ymax = lharanges[fp->ui->PDFflavor->currentIndex()][3];

  if (!fp->ui->logx->isChecked())
    {
      ymin = lharanges[fp->ui->PDFflavor->currentIndex()][4];
      ymax = lharanges[fp->ui->PDFflavor->currentIndex()][5];
    }
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
  const int Nrep = fp->fPDF->numberPDF();
  const int f = fp->ui->PDFflavor->currentIndex()-6;

  int memi = 0;
  int memf = Nrep;
  if (!fp->ui->checkBox->isChecked())
    {
      memf = 1;
      if (fp->ui->setmember->value() > Nrep || fp->ui->setmember->value() < 0) {
          finished();
          return;
        }
    }

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(fp->ui->title->text().toStdString().c_str());

  TLegend *leg = new TLegend(0.603448,0.673729,0.981322,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  int legindex = 0;

  double *x = new double[N];
  for (int i = 0; i < N; i++)
    {
      if (fp->ui->logx->isChecked()) x[i] = exp(log(xmin)+i*(log(xmax)-log(xmin))/(N-1));
      else x[i] = xmin+i*(xmax-xmin)/(N-1);
    }

  //////////////////////////////
  /// REPLICAS
  //////////////////////////////
  for (int r = memi; r <= memf; r++)
    {
      if(fIsTerminated){ fIsTerminated = false; return; }

      emit progress(r*100/memf);

      int rep = r;
      if (fp->fPDF->GetErrorType() == ER_NONE) rep = fp->fPDF->GetReplica();
      if (!fp->ui->checkBox->isChecked()) rep = fp->ui->setmember->value();
      fp->fPDF->initPDF(rep);

      TGraph *g = new TGraph(N);
      g->SetLineWidth(2);
      g->SetLineColor(colors[fp->ui->colorreplica->currentIndex()]);

      for (int i = 0; i < N; i++)
        g->SetPoint(i, x[i], fp->fPDF->GetFlvrPDF(x[i],Qf,f));

      mg->Add(g,"l");
      if (r == 1) { leg->AddEntry(g,"Members","l"); legindex++;}
    }

  //////////////////////////////
  /// Central values / Errors
  //////////////////////////////
  if (fp->ui->mean->isChecked() || fp->ui->stddev->isChecked() || fp->ui->cl->isChecked())
    {
      if(fIsTerminated){ return; }

      double *xPDF = new double[N];
      double *xPDFErr = new double[N];
      double *upErr = new double[N];
      double *dnErr = new double[N];
      fp->fPDF->GetFlvrPDFCVErr(N,x,Qf,f,xPDF,xPDFErr,upErr,dnErr);

      if(fIsTerminated){ fIsTerminated = false; return; }

      if (fp->ui->mean->isChecked())
        {
          TGraph *gcv = new TGraph(N);
          gcv->SetLineWidth(2);
          gcv->SetLineStyle(2);
          gcv->SetLineColor(colors[fp->ui->coloravg->currentIndex()]);

          for (int i = 0; i < N; i++)
            gcv->SetPoint(i, x[i], xPDF[i]);

          mg->Add(gcv,"l");
          leg->AddEntry(gcv,"Central value","l"); legindex++;
        }

      if(fIsTerminated){ fIsTerminated = false;  return; }

      if (fp->ui->stddev->isChecked())
        {
          TGraph *gstd = new TGraph(N);
          gstd->SetLineWidth(2);
          gstd->SetLineColor(colors[fp->ui->colorstddev->currentIndex()]);

          TGraph *gstd2 = new TGraph(N);
          gstd2->SetLineWidth(2);
          gstd2->SetLineColor(colors[fp->ui->colorstddev->currentIndex()]);

          for (int i = 0; i < N; i++)
            {
              gstd->SetPoint(i, x[i], xPDF[i]+xPDFErr[i]);
              gstd2->SetPoint(i, x[i], xPDF[i]-xPDFErr[i]);
            }

          mg->Add(gstd,"l");
          mg->Add(gstd2,"l");
          leg->AddEntry(gstd,"Std. deviation","l"); legindex++;
        }

      if(fIsTerminated){ fIsTerminated = false; return; }

      if (fp->ui->cl->isChecked())
        {
          TGraph *up = new TGraph(N);
          up->SetLineWidth(2);
          up->SetLineColor(colors[fp->ui->color68cl->currentIndex()]);

          TGraph *dn = new TGraph(N);
          dn->SetLineWidth(2);
          dn->SetLineColor(colors[fp->ui->color68cl->currentIndex()]);

          mg->Add(up,"l");
          mg->Add(dn,"l");
          leg->AddEntry(up,"68% c.l.","l"); legindex++;

          for (int i = 0; i < N; i++) {
            up->SetPoint(i, x[i], upErr[i]);
            dn->SetPoint(i, x[i], dnErr[i]);
          }
        }

      delete[] xPDF;
      delete[] xPDFErr;
      delete[] upErr;
      delete[] dnErr;
    }

  delete[] x;

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
  mg->GetYaxis()->SetRangeUser(ymin,ymax);

  leg->AddEntry("",TString("Q = " + QString::number(fp->ui->Qf->text().toDouble(),'g',3).toStdString() + " GeV"),"");
  if (legindex < 4)
    for (int i = legindex; i < 4; i++)
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

void memberthread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}

void memberthread::stop()
{
  fIsTerminated = true;
}
