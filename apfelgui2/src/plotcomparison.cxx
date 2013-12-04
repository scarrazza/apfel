#include "plotcomparison.h"
#include "ui_plotcomparison.h"
#include "pdfdialog.h"
#include "common.h"
#include <cmath>
#include <vector>

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
using namespace std;

PlotComparison::PlotComparison(QWidget *parent,std::vector<PDFDialog*> pdf) :
  QWidget(parent),
  ui(new Ui::PlotComparison),
  thread(NULL),
  fPDF(pdf),
  fPlotName("compareplot.png"),
  fRef(0),
  fIsRunning(false)
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
  if (thread)
    {
      thread->terminate();
      delete thread;
    }
}

void PlotComparison::on_playButton_clicked()
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

void PlotComparison::on_saveButton_clicked()
{
  QString selectedFilter;
  QString path = QFileDialog::getSaveFileName(this,
                                              tr("Save as"),"",
                                              tr(".eps;;.ps;;.pdf;;.png;;.root;;.C"),&selectedFilter);
  if (path != 0) thread->SaveCanvas(path + selectedFilter);
}

void PlotComparison::ThreadFinished()
{
  ui->playButton->setIcon(QIcon(":/images/StepForwardNormalBlue.png"));
  ui->playButton->setText("Compute");
  fIsRunning = false;

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
  fRef = index;
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
  fFileName(filename),
  fIsTerminated(false)
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

  TLegend *leg =  new TLegend(0.479885,0.673729,0.859195,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TMultiGraph *mg = new TMultiGraph();

  double *refx = new double[N];

  double *x = new double[N];
  for (int i = 0; i < N; i++)
    {
      if (fp->ui->logx->isChecked()) x[i] = exp(log(xmin)+i*(log(xmax)-log(xmin))/(N-1));
      else x[i] = xmin+i*(xmax-xmin)/(N-1);
    }

  vector<int> indRef;
  indRef.push_back(fp->fRef);
  for (int i = 0; i < (int) fp->fPDF.size(); i++)
    if (i != fp->fRef) indRef.push_back(i);

  for (int set = 0; set < (int) fp->fPDF.size(); set++)
    {
      emit progress(set*100/fp->fPDF.size());
      fp->fPDF[indRef[set]]->InitPDFset(Qi,Qf);
      const int f = fp->ui->PDFflavor->currentIndex()-6;

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

      double *xPDF = new double[N];
      double *xPDFErr = new double[N];
      double *upErr = new double[N];
      double *dnErr = new double[N];

      fp->fPDF[indRef[set]]->GetFlvrPDFCVErr(N,x,Qf,f,xPDF,xPDFErr,upErr,dnErr);

      for (int ix = 0; ix < N; ix++)
        {
          if(fIsTerminated){ fIsTerminated = false; return; }

          if (fp->ui->ratio->isChecked()) {
              if (set == 0) refx[ix] = xPDF[ix];
              xPDF[ix] /= refx[ix];
            }

          g->SetPoint(ix,x[ix], xPDF[ix]);
          gcv->SetPoint(ix, x[ix], xPDF[ix]);

          if (fp->ui->stddev->isChecked() && fp->fPDF[indRef[set]]->numberPDF() > 1) {
              if (fp->ui->ratio->isChecked()) xPDFErr[ix] /= refx[ix];
              g->SetPointError(ix,0, xPDFErr[ix]);
              gup->SetPoint(ix,x[ix],xPDF[ix]+xPDFErr[ix]);
              gdn->SetPoint(ix,x[ix],xPDF[ix]-xPDFErr[ix]);
            }
          else
            g->SetPointError(ix,0, 0);
        }

      delete[] xPDF;
      delete[] xPDFErr;
      delete[] upErr;
      delete[] dnErr;

      mg->Add(g,"le3");
      if (set == 0)
        mg->Add(gcv,"l");

      if (fp->ui->stddev->isChecked())
        {
          mg->Add(gup,"l");
          mg->Add(gdn,"l");
        }

      leg->AddEntry(g,TString(fp->fPDF[indRef[set]]->PDFname().toStdString()),"fl");
    }

  delete[] refx;
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

void comparisonthread::stop()
{
  fIsTerminated = true;
}
