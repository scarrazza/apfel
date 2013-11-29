#include "plotcomparison.h"
#include "ui_plotcomparison.h"
#include "pdfdialog.h"
#include "common.h"
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

PlotComparison::PlotComparison(QWidget *parent,std::vector<PDFDialog*> pdf) :
  QWidget(parent),
  ui(new Ui::PlotComparison),
  thread(NULL),
  fPDF(pdf),
  fPlotName("compareplot.png")
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  long int t = static_cast<long int> (time(NULL));
  QString str;
  str.append(QString("%1").arg(t));

  fPlotName = "comparisonplot_" + str + ".png";
  thread  = new comparisonthread(this,fPlotName);

  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));

  //ui->graphicsView->scale(1.2,1.2);

  ui->Qi->setEnabled(false);
  for (int i = 0; i < (int) fPDF.size(); i++)
    if (!fPDF[i]->isLHAPDF())
      ui->Qi->setEnabled(true);

  ui->xtitle->setText("x");
  ui->ytitle->setText("");
  ui->title->setText(name[ui->PDFflavor->currentIndex()]+ ", comparison plot");

  for (int i = 0; i < (int) fPDF.size(); i++)
    ui->referenceSet->addItem(fPDF[i]->PDFname());
}

PlotComparison::~PlotComparison()
{
  delete ui;
  delete thread;
}

void PlotComparison::on_playButton_clicked()
{
  ui->playButton->setEnabled(false);
  QApplication::processEvents();
  thread->start();
}

void PlotComparison::on_saveButton_clicked()
{
  QString selectedFilter;
  QString path = QFileDialog::getSaveFileName(this,
                                              tr("Save as"),"",
                                              tr(".eps;;.ps;;.pdf;;.png;;.root"),&selectedFilter);
  if (path != 0) thread->SaveCanvas(path + selectedFilter);
}

void PlotComparison::ThreadFinished()
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

void PlotComparison::ThreadProgress(int i)
{
  ui->progressBar->setValue(i);
}

void PlotComparison::on_PDFflavor_currentIndexChanged(int index)
{
  ui->title->setText(name[index]+ ", comparison plot");
}

void PlotComparison::on_referenceSet_currentIndexChanged(int index)
{
  PDFDialog *tmp = fPDF[index];
  PDFDialog *pre = fPDF[0];

  fPDF[0] = tmp;
  fPDF[index] = pre;
}

void PlotComparison::on_ratio_toggled(bool checked)
{
  if (checked == true)
    ui->ytitle->setText("Ratio");
  else
    ui->ytitle->setText("");
}

void PlotComparison::on_automaticrange_toggled(bool checked)
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

comparisonthread::comparisonthread(QObject *parent, QString filename):
  QThread(parent),
  fp((PlotComparison*)parent),
  fFileName(filename)
{
}

comparisonthread::~comparisonthread()
{
}

void comparisonthread::run()
{
  fC = new TCanvas();
  fC->SetTickx();
  fC->SetTicky();

  if (fp->ui->log->isChecked())
    fC->SetLogx();

  // Initialize PDFs

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

  if (fp->ui->ratio->isChecked())
    {
      ymin = 0.5;
      ymax = 1.6;
    }

  if (!fp->ui->automaticrange->isChecked())
    {
      xmin = fp->ui->xmin->text().toDouble();
      xmax = fp->ui->xmax->text().toDouble();
      ymin = fp->ui->ymin->text().toDouble();
      ymax = fp->ui->ymax->text().toDouble();
    }

  const double Qi = fp->ui->Qi->text().toDouble();
  const double Qf = fp->ui->Qf->text().toDouble();

  TLegend *leg = NULL;
  if (fp->ui->log->isChecked())
    leg = new TLegend(0.479885,0.673729,0.859195,0.883475);
  else
    leg = new TLegend(0.12931,0.673729,0.507184,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TMultiGraph *mg = new TMultiGraph();

  double *refx = new double[N];
  for (int set = 0; set < (int) fp->fPDF.size(); set++)
    {
      fp->fPDF[set]->InitPDFset(Qi,Qf);
      const int Nrep = fp->fPDF[set]->numberPDF();
      const int f = fp->ui->PDFflavor->currentIndex()-6;

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

      double **xPDF = new double*[Nrep];
      for (int r = memi; r <= memf; r++)
        {
          emit progress((r-1)*100/memf);

          xPDF[r-memi] = new double[N];
          fp->fPDF[set]->initPDF(r);

          for (int ix = 0; ix < N; ix++)
            {
              double x = 0;
              if (fp->ui->log->isChecked())
                x = exp(log(xmin)+ix*(log(xmax)-log(xmin)/N));
              else
                x = xmin+ix*(xmax-xmin)/N;

              xPDF[r-memi][ix] = fp->fPDF[set]->GetFlvrPDFCV(x,Qf,f);
            }
        }

      for (int ix = 0; ix < N; ix++)
        {
          double x = 0;
          if (fp->ui->log->isChecked())
            x = exp(log(xmin)+ix*(log(xmax)-log(xmin)/N));
          else
            x = xmin+ix*(xmax-xmin)/N;

          double xf = xPDF[0][ix];
          if (memf > 1) xf = ComputeAVG(Nrep, ix, xPDF);

          if (fp->ui->ratio->isChecked()) {
              if (set == 0) refx[ix] = xf;
              xf /= refx[ix];
            }

          g->SetPoint(ix,x, xf);
          gcv->SetPoint(ix, x, xf);

          if (fp->ui->stddev->isChecked() && memf > 1) {
              double xferr = ComputeStdDev(Nrep, ix, xPDF);
              if (fp->ui->ratio->isChecked()) xferr /= refx[ix];
              g->SetPointError(ix,0, xferr);
              gup->SetPoint(ix,x,xf+xferr);
              gdn->SetPoint(ix,x,xf-xferr);
            }
          else
            g->SetPointError(ix,0, 0);
        }

      mg->Add(g,"le3");
      if (set == 0)
        mg->Add(gcv,"l");

      if (fp->ui->stddev->isChecked())
        {
          mg->Add(gup,"l");
          mg->Add(gdn,"l");
        }

      leg->AddEntry(g,TString(fp->fPDF[set]->PDFname().toStdString()),"fl");

      for (int i = 0; i < memf; i++)
        if (xPDF[i]) delete[] xPDF[i];
      delete[] xPDF;

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

  leg->AddEntry("",TString("Q = " + QString::number(fp->ui->Qf->text().toDouble(),'g',3).toStdString() + " GeV"),"");
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

void comparisonthread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}
