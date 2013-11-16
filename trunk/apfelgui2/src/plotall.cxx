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
  fPlotName("allplot.svg")
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
  delete thread;
}

void PlotAll::on_playButton_clicked()
{
  ui->playButton->setEnabled(false);
  QApplication::processEvents();
  thread->start();
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
  QString path;
  QFileDialog d(this,tr("Save as"),"",tr(".eps;;.ps;;.pdf;;.png;;.root"));

  if (d.exec())
    {
      path = d.selectedFiles()[0];
      if(path != 0) thread->SaveCanvas(QString(path + d.selectedNameFilter()));
    }
}

plotthread::plotthread(QObject *parent, QString filename):
  QThread(parent),
  fp((PlotAll*)parent),
  fFileName(filename)
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

  if (fp->ui->log->isChecked())
    fC->SetLogx();

  // Initialize PDFs

  const int N = fp->ui->xpoints->value();
  double xmin = 1e-3;
  double xmax = 1;
  double ymin = 0;
  double ymax = 1;;
  int legindex = 0;

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

  TLegend *leg = NULL;
  if (fp->ui->log->isChecked())
    leg = new TLegend(0.130747,0.491525,0.300287,0.883475);
  else
    leg = new TLegend(0.716954,0.491525,0.886494,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TMultiGraph *mg = new TMultiGraph();

  int nf = fp->ui->maxFlavors->value();
  for (int fl = -nf; fl <= nf; fl++)
    {
      emit progress( (fl+nf)*100./(2*nf+1));

      TGraphErrors *g = new TGraphErrors(N);
      g->SetLineWidth(2);
      g->SetLineColor(colors2[fl+6]);
      g->SetFillColor(colors2[fl+6]);

      double **xPDF = new double*[memf];
      for (int r = memi; r <= memf; r++)
        {
          xPDF[r-memi] = new double[N];
          if (memf == 1)
            fp->fPDF->initPDF(fp->ui->setmember->value());
          else
            fp->fPDF->initPDF(r);

          for (int ix = 0; ix < N; ix++)
            {
              double x = 0;
              if (fp->ui->log->isChecked())
                x = exp(log(xmin)+ix*(log(xmax)-log(xmin)/N));
              else
                x = xmin+ix*(xmax-xmin)/N;

              xPDF[r-memi][ix] = fp->fPDF->GetFlvrPDF(x,Qf,fl);
              if (fl == 0)
                xPDF[r-memi][ix] /= 10;
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

          if (fp->ui->checkBox->isChecked())
            xf = ComputeAVG(Nrep, ix, xPDF);
          g->SetPoint(ix,x, xf);

          if (fp->ui->stddev->isChecked())
            g->SetPointError(ix,0, ComputeStdDev(Nrep, ix, xPDF));
          else
            g->SetPointError(ix,0, 0);
        }

      if (fl != 0)
        leg->AddEntry(g,name[fl+6].toStdString().c_str(),"l");
      else
        leg->AddEntry(g,"xg(x,Q)/10","l");
      legindex++;

      mg->Add(g,"le3");

      for (int i = 0; i < memf; i++)
        if (xPDF[i]) delete[] xPDF[i];
      delete[] xPDF;
    }

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
