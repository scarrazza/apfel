#include "plotlumi.h"
#include "ui_plotlumi.h"
#include "pdfdialog.h"
#include "common.h"
#include "utils.h"
#include <cmath>

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

PlotLumi::PlotLumi(QWidget *parent,std::vector<PDFDialog*> pdf) :
  QWidget(parent),
  ui(new Ui::PlotLumi),
  thread(NULL),
  fPDF(pdf),
  fPlotName("lumiplot.png")
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  long int t = static_cast<long int> (time(NULL));
  QString str;
  str.append(QString("%1").arg(t));

  fPlotName = "luminplot_" + str + ".png";
  thread  = new lumithread(this,fPlotName);

  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));

  //ui->graphicsView->scale(1.2,1.2);

  ui->Qi->setEnabled(false);
  for (int i = 0; i < (int) fPDF.size(); i++)
    if (!fPDF[i]->isLHAPDF())
      ui->Qi->setEnabled(true);

  ui->xtitle->setText("x");
  ui->ytitle->setText("Ratio");
  ui->title->setText("Gluon - Gluon Luminosity");

  for (int i = 0; i < (int) fPDF.size(); i++)
    ui->referenceSet->addItem(fPDF[i]->PDFname());
}

PlotLumi::~PlotLumi()
{
  delete ui;
  delete thread;
}

void PlotLumi::on_referenceSet_currentIndexChanged(int index)
{
  PDFDialog *tmp = fPDF[index];
  PDFDialog *pre = fPDF[0];

  fPDF[0] = tmp;
  fPDF[index] = pre;
}

void PlotLumi::on_automaticrange_toggled(bool checked)
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

void PlotLumi::on_playButton_clicked()
{
  ui->playButton->setEnabled(false);
  QApplication::processEvents();
  thread->start();
}

void PlotLumi::on_saveButton_clicked()
{
  QString selectedFilter;
  QString path = QFileDialog::getSaveFileName(this,
                                              tr("Save as"),"",
                                              tr(".eps;;.ps;;.pdf;;.png;;.root"),&selectedFilter);
  if (path != 0) thread->SaveCanvas(path + selectedFilter);
}


void PlotLumi::ThreadFinished()
{
  ui->playButton->setEnabled(true);
  ui->saveButton->setEnabled(true);

  // plot to canvas
  QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);
  //QGraphicsSvgItem * item = new QGraphicsSvgItem(fPlotName);
  QGraphicsPixmapItem *item = new QGraphicsPixmapItem(fPlotName);
  scene->addItem(item);
  ui->graphicsView->setScene(scene);
  ui->graphicsView->show();

  ui->progressBar->setValue(0);

  QFile::remove(fPlotName) ;
}

void PlotLumi::ThreadProgress(int i)
{
  ui->progressBar->setValue(i);
}

void PlotLumi::on_PDFflavor_currentIndexChanged(int index)
{
  QString titleY [] = {
  "Gluon - Gluon Luminosity",
  "Quark - Antiquark Luminosity",
  "Quark - Gluon Luminosity",
  "Charm - Anticharm Luminosity",
  "Bottom - Antibottom Luminosity",
  "Gluon - Charm Luminosity",
  "Bottom - Gluon Luminosity",
  "Quark - Quark Luminosity",
  "Photon - Photon Luminosity" };

  ui->ytitle->setText(titleY[index]);
}

lumithread::lumithread(QObject *parent, QString filename):
  QThread(parent),
  fp((PlotLumi*)parent),
  fFileName(filename)
{
}

lumithread::~lumithread()
{
}

void lumithread::run()
{
  fC = new TCanvas();
  fC->SetTickx();
  fC->SetTicky();

  if (fp->ui->log->isChecked())
    fC->SetLogx();

  // Initialize PDFs

  const int N = fp->ui->xpoints->value();
  double xmin = 10;
  double xmax = 6e3;
  double ymin = 0.8;
  double ymax = 1.3;
  double eps = fp->ui->integration->text().toDouble();
  double S = pow(fp->ui->cmenergy->text().toDouble(), 2.0);

  vector<string> lumis;
  lumis.push_back("GG");
  lumis.push_back("QQ");
  lumis.push_back("QG");
  lumis.push_back("GC");
  lumis.push_back("BG");
  lumis.push_back("CC");
  lumis.push_back("BB");
  lumis.push_back("Q2");
  lumis.push_back("PP");
  lumis.push_back("PG");

  if (!fp->ui->automaticrange->isChecked())
    {
      xmin = fp->ui->xmin->text().toDouble();
      xmax = fp->ui->xmax->text().toDouble();
      ymin = fp->ui->ymin->text().toDouble();
      ymax = fp->ui->ymax->text().toDouble();
    }

  const double Qi = fp->ui->Qi->text().toDouble();

  TLegend *leg = new TLegend(0.12931,0.673729,0.507184,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);   

  TMultiGraph *mg = new TMultiGraph();

  double *refx = new double[N];
  for (int set = 0; set < (int) fp->fPDF.size(); set++)
    {
      fp->fPDF[set]->InitPDFset(Qi,xmax);
      const int Nrep = fp->fPDF[set]->numberPDF();

      int memi = 1;
      int memf = Nrep;

      TGraphErrors *g = new TGraphErrors(N);
      g->SetLineWidth(2);
      g->SetLineStyle(2);
      g->SetLineColor(fillcolor[set]+2);
      g->SetFillColor(fillcolor[set]);
      g->SetFillStyle(fillStyle[set]);

      TGraph *gcv = new TGraph(N);
      gcv->SetLineWidth(2);
      gcv->SetLineStyle(2);
      gcv->SetLineColor(fillcolor[set]+2);

      TGraph *gup = new TGraph(N);
      gup->SetLineWidth(2);
      gup->SetLineStyle(1);
      gup->SetLineColor(fillcolor[set]+2);

      TGraph *gdn = new TGraph(N);
      gdn->SetLineWidth(2);
      gdn->SetLineStyle(1);
      gdn->SetLineColor(fillcolor[set]+2);

      double *xmH = new double[N];
      for (int imH = 1; imH <= N; imH++)
        {
          emit progress((imH-1)*100/N);
          double mH = xmin * pow(xmax/xmin, double(imH-1)/(N-1));

          double *ggflux = new double[memf];

          for (int r = memi; r <= memf; r++)
            {
              fp->fPDF[set]->initPDF(r);

              ggflux[r-memi] = fp->fPDF[set]->getLum(mH, S, lumis[fp->ui->PDFflavor->currentIndex()],eps);
              xmH[imH-1]   = mH; //sqrt(mH*mH/S);

            }

	  double cv = ggflux[0];
	  if (fp->ui->stddev->isChecked()) cv = ComputeAVG(memf,ggflux);
	  if (set == 0) refx[imH-1] = cv;
	  cv /= refx[imH-1];

	  double cverr = 0;
	  if (fp->ui->stddev->isChecked() && memf > 1) cverr = ComputeStdDev(memf,ggflux);
	  cverr /= refx[imH-1];

	  g->SetPoint(imH-1, xmH[imH-1], cv);
	  g->SetPointError(imH-1,0,cverr);

	  gcv->SetPoint(imH-1,xmH[imH-1], cv);
	  if (fp->ui->stddev->isChecked())
	    {
	      gup->SetPoint(imH-1,xmH[imH-1], cv+cverr);
	      gdn->SetPoint(imH-1,xmH[imH-1], cv-cverr);
	    }

          delete[] ggflux;
        }

      delete[] xmH;

      mg->Add(g,"le3");
      if (set == 0)
        mg->Add(gcv,"l");

      if (fp->ui->stddev->isChecked())
        {
          mg->Add(gup,"l");
          mg->Add(gdn,"l");
        }

      leg->AddEntry(g,TString(fp->fPDF[set]->PDFname().toStdString()),"fl");

    }

  delete[] refx;

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
  mg->GetYaxis()->SetRangeUser(ymin,ymax);

  leg->AddEntry("",TString("#sqrt{S} = " + QString::number(fp->ui->cmenergy->text().toDouble(),'g',3).toStdString() + " GeV"),"");
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

void lumithread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}
