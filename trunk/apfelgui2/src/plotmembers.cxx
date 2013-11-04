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

#include "LHAPDF/LHAPDF.h"
#include "APFEL/APFEL.h"

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

  TMultiGraph *mg = new TMultiGraph();

  TLegend *leg = new TLegend(0.603448,0.673729,0.981322,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  int legindex = 0;

  if (fp->fPDF->isLHAPDF())
    {
      LHAPDF::initPDFSetByName(fp->fPDF->PDFname().toStdString());
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

      const int Nrep = LHAPDF::numberPDF();
      double xPDF[Nrep][N];

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

      for (int i = 0; i < Nrep; i++)
        {
          for (int j = 0; j < N; j++)
            {
              double x = 0;
              if (fp->ui->log->isChecked())
                x= exp(log(xmin)+j*(log(xmax)-log(xmin)/N));
              else
                x= xmin+j*(xmax-xmin)/N;

              qDebug() << LHAPDF::xfx((double)x,(double)(fp->ui->Qf->text()).toDouble(),(int)(fp->ui->PDFflavor->currentIndex()-6));
        }
        }

      /*
      for (int i = memi; i <= memf; i++)
        {
          emit progress((i-1)*100/memf);

          if (memf == 1)
            LHAPDF::initPDF(fp->ui->setmember->value());
          else
            LHAPDF::initPDF(i);

          TGraph *g = new TGraph(N);          
          g->SetLineWidth(2);
          g->SetLineColor(kGreen);

          for (int j = 0; j < N; j++)
            {
              double x = 0;
              if (fp->ui->log->isChecked())
                x= exp(log(xmin)+j*(log(xmax)-log(xmin)/N));
              else
                x= xmin+j*(xmax-xmin)/N;

              if (!LHAPDF::hasPhoton())
                {                  
                  xPDF[i-1][j] = LHAPDF::xfx(x,(fp->ui->Qf->text()).toDouble(),fp->ui->PDFflavor->currentIndex()-6);
                  g->SetPoint(j, x, LHAPDF::xfx(x,(fp->ui->Qf->text()).toDouble(),fp->ui->PDFflavor->currentIndex()-6));
                }
              else
                {
                  xPDF[i-1][j] = LHAPDF::xfxphoton(x,(fp->ui->Qf->text()).toDouble(),fp->ui->PDFflavor->currentIndex()-6);
                  g->SetPoint(j, x, LHAPDF::xfxphoton(x,(fp->ui->Qf->text()).toDouble(),fp->ui->PDFflavor->currentIndex()-6));
                }
            }
          mg->Add(g,"l");
          if (i == 1) { leg->AddEntry(g,"Replica","l"); legindex++;}
        }

      if (fp->ui->mean->isChecked())
        {
          TGraph *gcv = new TGraph(N);
          gcv->SetLineWidth(2);
          gcv->SetLineStyle(2);
          gcv->SetLineColor(kRed);

          for (int i = 0; i < N; i++)
            {
              double sum = 0;
              for (int j = 0; j < Nrep; j++)
                sum += xPDF[j][i];

              double x = 0;
              if (fp->ui->log->isChecked())
                x= exp(log(xmin)+i*(log(xmax)-log(xmin)/N));
              else
                x= xmin+i*(xmax-xmin)/N;

              gcv->SetPoint(i, x, sum/Nrep);
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
              double sum = 0;
              for (int j = 0; j < Nrep; j++)
                sum += xPDF[j][i];
              sum/=Nrep;

              double sum2 = 0;
              for (int j = 0; j < Nrep; j++)
                sum2 += pow(xPDF[j][i] - sum,2);
              sum2 /= Nrep-1;
              sum2 = sqrt(sum2);

              double x = 0;
              if (fp->ui->log->isChecked())
                x= exp(log(xmin)+i*(log(xmax)-log(xmin)/N));
              else
                x= xmin+i*(xmax-xmin)/N;

              gstd->SetPoint(i, x, sum+sum2);
              gstd2->SetPoint(i, x, sum-sum2);
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

    }
  else
    {
      // Initialize apfel
      APFEL::SetQLimits(0e0,1e4);
      APFEL::SetPerturbativeOrder(fp->fPDF->ptord());

      if (fp->fPDF->scheme() == 0)
        APFEL::SetVFNS();
      else
        APFEL::SetFFNS(fp->fPDF->nf());

      APFEL::SetTheory(fp->fPDF->theory().toStdString());
/*
      APFEL::SetAlphaQCDRef(ui->lineEdit_13->text().toDouble(),ui->lineEdit_14->text().toDouble());

      APFEL::SetAlphaQEDRef(ui->lineEdit_15->text().toDouble(),ui->lineEdit_16->text().toDouble());

      if (ui->comboBox_4->currentIndex() == 0)
        APFEL::SetPoleMasses(ui->lineEdit_10->text().toDouble(),
                             ui->lineEdit_11->text().toDouble(),
                             ui->lineEdit_12->text().toDouble());
      else
        APFEL::SetMSbarMasses(ui->lineEdit_10->text().toDouble(),
                              ui->lineEdit_11->text().toDouble(),
                              ui->lineEdit_12->text().toDouble());

      APFEL::SetMaxFlavourPDFs(ui->spinBox_2->value());
      APFEL::SetMaxFlavourAlpha(ui->spinBox_3->value());

      APFEL::SetRenFacRatio(ui->doubleSpinBox_10->value());

      APFEL::SetPDFSet(ui->lineEdit->text().toStdString());
      APFEL::SetReplica(ui->spinBox_4->value());

      APFEL::InitializeAPFEL();

      APFEL::EvolveAPFEL(fp->ui->Qi, fp->ui->Qf);
*/
    }

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
  l.DrawLatex(0.95,0.15,TString("Generated by APFEL: V.Bertone, S.Carrazza, J.Rojo (arXiv:1310.1394)"));

  fC->SaveAs("memberplot.svg");
}

void memberthread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}
