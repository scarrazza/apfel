#include "inc/plotmembers.h"
#include "ui_plotmembers.h"
#include "pdfdialog.h"
#include <cmath>
#include <algorithm>
#include <vector>
using namespace std;

#include <QGraphicsScene>
#include <QGraphicsItemGroup>
#include <QGraphicsSvgItem>
#include <QFile>
#include <QtSvg/QSvgWidget>
#include <QDebug>
#include <QFileDialog>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"

#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"

double ComputeAVG(int n, int ix, double **x)
{
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    sum += x[i][ix];
  return sum / n;
}

double ComputeStdDev(int n, int ix, double **x)
{
  double sum = 0.0;
  double avg = ComputeAVG(n, ix, x);
  for (int i = 0; i < n; i++)
    sum += (x[i][ix]-avg)*(x[i][ix]-avg);

  sum /= n-1;

  return sqrt(sum);
}

PlotMembers::PlotMembers(QWidget *parent, PDFDialog *pdf) :
  QWidget(parent),
  ui(new Ui::PlotMembers),
  thread(NULL),
  fPDF(pdf)
{
  ui->setupUi(this);

  thread  = new memberthread(this);

  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));

  ui->graphicsView->scale(1.2,1.2);

  if (fPDF->isLHAPDF()) ui->Qi->setEnabled(false);
}

PlotMembers::~PlotMembers()
{
  delete ui;
  delete thread;
}

void PlotMembers::on_playButton_clicked()
{
  ui->playButton->setEnabled(false);
  QApplication::processEvents();
  thread->start();
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
  ui->playButton->setEnabled(true);
  ui->saveButton->setEnabled(true);

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);
  QGraphicsSvgItem * item = new QGraphicsSvgItem("memberplot.svg");
  scene->addItem(item);
  ui->graphicsView->setScene(scene);
  ui->graphicsView->show();

  ui->progressBar->setValue(0);

  QFile::remove("memberplot.svg") ;
}

void PlotMembers::ThreadProgress(int i)
{
  ui->progressBar->setValue(i);
}

void PlotMembers::on_saveButton_clicked()
{
  QString path;
  path = QFileDialog::getSaveFileName(this,tr("Save as"),"",tr("*.eps (*.eps);;All files (*.*)"));

  if(path != 0) thread->SaveCanvas(path + ".eps");

}

memberthread::memberthread(QObject *parent):
  QThread(parent),
  fp((PlotMembers*)parent)
{
}

memberthread::~memberthread()
{
}

double lharanges[][6] = {
  {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // tbar
  {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // bbar
  {1e-5, 1.0, -1.0, 1.0, -0.1, 0.1}, // cbar
  {1e-5, 1.0, -1.5, 2.5, -0.1, 0.1}, // sbar
  {1e-5, 1.0, -0.1, 1.3, -0.1, 0.4}, // ubar
  {1e-5, 1.0, -0.1, 1.3, -0.1, 0.4}, // dbar
  {1e-5, 1.0, -2.5, 7.0, -0.1, 0.6}, // g
  {1e-5, 1.0, -0.1, 1.3, -0.1, 0.6}, // d
  {1e-5, 1.0, -0.1, 1.3, -0.1, 1.0}, // u
  {1e-5, 1.0, -1.5, 2.5, -0.1, 0.1}, // s
  {1e-3, 1.0, -0.1, 1.0, -0.1, 0.1}, // c
  {1e-3, 1.0, -0.1, 1.0, -0.1, 0.5}, // b
  {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // t
  {1e-5, 1.0, -3.0, 3.0, -1.0, 1.0}, // photon
                        };

QString name[] = {
  "x#bar{t}(x,Q)",
  "x#bar{b}(x,Q)",
  "x#bar{c}(x,Q)",
  "x#bar{s}(x,Q)",
  "x#bar{u}(x,Q)",
  "x#bar{d}(x,Q)",
  "xg(x,Q)",
  "xd(x,Q)",
  "xu(x,Q)",
  "xs(x,Q)",
  "xc(x,Q)",
  "xb(x,Q)",
  "xt(x,Q)",
  "x#gamma(x,Q)"
};

void memberthread::run()
{
  fC = new TCanvas();
  fC->SetTickx();
  fC->SetTicky();

  if (fp->ui->log->isChecked())
    fC->SetLogx();

  // Initialize PDFs
  fp->fPDF->InitPDFset();

  const int N = fp->ui->xpoints->value();
  double xmin = lharanges[fp->ui->PDFflavor->currentIndex()][0];
  double xmax = lharanges[fp->ui->PDFflavor->currentIndex()][1];
  double ymin = lharanges[fp->ui->PDFflavor->currentIndex()][2];
  double ymax = lharanges[fp->ui->PDFflavor->currentIndex()][3];

  if (fp->ui->lin->isChecked())
    {
      ymin = lharanges[fp->ui->PDFflavor->currentIndex()][4];
      ymax = lharanges[fp->ui->PDFflavor->currentIndex()][5];
    }

  if (!fp->ui->automaticrange->isChecked())
    {
      xmin = fp->ui->xmin->text().toDouble();
      xmax = fp->ui->xmax->text().toDouble();
      ymin = fp->ui->ymin->text().toDouble();
      ymax = fp->ui->ymax->text().toDouble();
    }

  const int Nrep = fp->fPDF->numberPDF();
  const double Qi = fp->ui->Qi->text().toDouble();
  const double Qf = fp->ui->Qf->text().toDouble();
  const int f = fp->ui->PDFflavor->currentIndex()-6;
  double **xPDF = new double*[Nrep];
  for (int i = 0; i < Nrep; i++)
    {
      xPDF[i] = new double[N];
      for (int j = 0; j < N; j++)
        xPDF[i][j] = 0.0;
    }


  int memi = 1;
  int memf = Nrep;
  if (!fp->ui->checkBox->isChecked())
    {
      memf = 1;
      if (fp->ui->setmember->value() > Nrep || fp->ui->setmember->value() < 0) {
          finished();
          return;
        }
    }

  TLegend *leg = new TLegend(0.603448,0.673729,0.981322,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  int legindex = 0;

  TMultiGraph *mg = new TMultiGraph();

  for (int r = memi; r <= memf; r++)
    {
      emit progress((r-1)*100/memf);

      if (memf == 1)
        fp->fPDF->initPDF(fp->ui->setmember->value());
      else
        fp->fPDF->initPDF(r);

      TGraph *g = new TGraph(N);
      g->SetLineWidth(2);
      g->SetLineColor(kGreen);

      for (int ix = 0; ix < N; ix++)
        {
          double x = 0;          
          if (fp->ui->log->isChecked())
            x = exp(log(xmin)+ix*(log(xmax)-log(xmin)/N));
          else
            x = xmin+ix*(xmax-xmin)/N;

          xPDF[r-1][ix] = fp->fPDF->GetFlvrPDF(x,Qi,Qf,f);
          g->SetPoint(ix, x, xPDF[r-1][ix]);

        }
      mg->Add(g,"l");
      if (r == 1) { leg->AddEntry(g,"Replica","l"); legindex++;}
    }

  if (fp->ui->mean->isChecked())
    {
      TGraph *gcv = new TGraph(N);
      gcv->SetLineWidth(2);
      gcv->SetLineStyle(2);
      gcv->SetLineColor(kRed);

      for (int i = 0; i < N; i++)
        {
          double x = 0;
          if (fp->ui->log->isChecked())
            x= exp(log(xmin)+i*(log(xmax)-log(xmin)/N));
          else
            x= xmin+i*(xmax-xmin)/N;

          gcv->SetPoint(i, x, ComputeAVG(Nrep,i,xPDF));
        }
      mg->Add(gcv,"l");
      leg->AddEntry(gcv,"Mean value","l");
      legindex++;
    }

  if (fp->ui->stddev->isChecked())
    {
      TGraph *gstd = new TGraph(N);
      gstd->SetLineWidth(2);
      gstd->SetLineColor(kBlue);

      TGraph *gstd2 = new TGraph(N);
      gstd2->SetLineWidth(2);
      gstd2->SetLineColor(kBlue);

      for (int i = 0; i < N; i++)
        {
          double x = 0;
          if (fp->ui->log->isChecked())
            x= exp(log(xmin)+i*(log(xmax)-log(xmin)/N));
          else
            x= xmin+i*(xmax-xmin)/N;

          double avg = ComputeAVG(Nrep,i,xPDF);
          double stddev = ComputeStdDev(Nrep,i,xPDF);

          gstd->SetPoint(i, x, avg+stddev);
          gstd2->SetPoint(i, x, avg-stddev);
        }
      mg->Add(gstd,"l");
      mg->Add(gstd2,"l");
      leg->AddEntry(gstd,"Std. deviation","l");
      legindex++;
    }

  if (fp->ui->cl->isChecked())
    {
      TGraph *up = new TGraph(N);
      up->SetLineWidth(2);
      TGraph *dn = new TGraph(N);
      dn->SetLineWidth(2);
      for (int i = 0; i < N; i++)
        {
          double x = 0;
          if (fp->ui->log->isChecked())
            x= exp(log(xmin)+i*(log(xmax)-log(xmin)/N));
          else
            x= xmin+i*(xmax-xmin)/N;                  

          vector<double> xval;
          for (int r = 0; r < Nrep; r++)
            xval.push_back(xPDF[r][i]);
          sort(xval.begin(),xval.end());

          int esc = Nrep*(1-0.68)/2;

          up->SetPoint(i, x, xval[Nrep-esc-1]);
          dn->SetPoint(i, x, xval[esc]);
        }
      mg->Add(up,"l");
      mg->Add(dn,"l");
      leg->AddEntry(up,"68% c.l.","l");
      legindex++;
    }

  mg->SetTitle(TString(name[fp->ui->PDFflavor->currentIndex()].toStdString() + ", " + fp->fPDF->PDFname().toStdString() + " members"));

  mg->Draw("AL");

  mg->GetXaxis()->SetTitle("x");
  mg->GetXaxis()->CenterTitle(true);

  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetLabelSize(0.05);

  mg->GetXaxis()->SetLimits(xmin,xmax);
  mg->GetYaxis()->SetRangeUser(ymin,ymax);

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

  fC->SaveAs("memberplot.svg");

  for(int i = 0; i < Nrep; i++)
    if (xPDF[i]) delete[] xPDF[i];
  delete[] xPDF;

}

void memberthread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}
