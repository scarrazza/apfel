#include <iostream>
#include "plotdis.h"
#include "ui_plotdis.h"
#include "pdfdialog.h"
#include "common.h"
#include <cmath>

#include <QGraphicsScene>
#include <QGraphicsItemGroup>
//#include <QGraphicsSvgItem>
#include <QFile>
//#include <QtSvg/QSvgWidget>
#include <QDebug>
#include <QFileDialog>
#include <QDesktopWidget>

#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
using namespace std;

PlotDIS::PlotDIS(QWidget *parent, PDFDialog *pdf) :
  QWidget(parent),
  ui(new Ui::PlotDIS),
  thread(NULL),
  fPDF(pdf),
  fPlotName("displot.png"),
  fIsRunning(false)
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  long int t = static_cast<long int> (time(NULL));
  QString str;
  str.append(QString("%1").arg(t));

  fPlotName = "displot_" + str + ".png";
  thread  = new disthread(this,fPlotName);

  connect(thread, SIGNAL(finished()), this, SLOT(ThreadFinished()));
  connect(thread, SIGNAL(progress(int)), this, SLOT(ThreadProgress(int)));

  if (fPDF->isLHAPDF()) ui->Qi->setEnabled(false);
  if (fPDF->numberPDF() == 1)
    {
      ui->stddev->setEnabled(false);
      ui->stddev->setChecked(false);
    }

  ui->xtitle->setText("x");
  ui->ytitle->setText("");
  ui->title->setText( "F_{2}^{c}(x), " + fPDF->PDFname());
}

PlotDIS::~PlotDIS()
{
  delete ui;
  if (thread)
    {
      thread->terminate();
      delete thread;
    }
}

void PlotDIS::on_process_currentIndexChanged(int index)
{
  if (index == 0)
    {
      ui->projectile->removeItem(2);
      ui->projectile->removeItem(2);
      ui->target->removeItem(3);
    }
  else if (index == 1)
    {
      ui->projectile->removeItem(2);
      ui->projectile->removeItem(2);
      ui->target->removeItem(3);
    }
  else if (index == 2)
    {      
      ui->projectile->addItem("Neutrino");
      ui->projectile->addItem("Antineutrino");
      ui->target->addItem("Iron");
    }
}

void PlotDIS::on_automaticrange_toggled(bool checked)
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

void PlotDIS::on_playButton_clicked()
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

void PlotDIS::on_saveButton_clicked()
{
  QString selectedFilter;
  QString path = QFileDialog::getSaveFileName(this,
                                              tr("Save as"),"",
                                              tr(".eps;;.ps;;.pdf;;.png;;.root;;.C"),&selectedFilter);
  if (path != 0) thread->SaveCanvas(path + selectedFilter);

}

void PlotDIS::ThreadFinished()
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

void PlotDIS::ThreadProgress(int i)
{
  ui->progressBar->setValue(i);
}

void PlotDIS::on_output_currentIndexChanged(int index)
{
  if (index == 0)
    ui->title->setText( "F_{2}^{l}(x), " + fPDF->PDFname());
  else if (index == 1)
    ui->title->setText( "F_{2}^{c}(x), " + fPDF->PDFname());
  else if (index == 2)
    ui->title->setText( "F_{2}^{b}(x), " + fPDF->PDFname());
  else if (index == 3)
    ui->title->setText( "F_{2}^{t}(x), " + fPDF->PDFname());
  else if (index == 4)
    ui->title->setText( "F_{2}^{p}(x), " + fPDF->PDFname());
  else if (index == 5)
    ui->title->setText( "F_{L}^{l}(x), " + fPDF->PDFname());
  else if (index == 6)
    ui->title->setText( "F_{L}^{c}(x), " + fPDF->PDFname());
  else if (index == 7)
    ui->title->setText( "F_{L}^{b}(x), " + fPDF->PDFname());
  else if (index == 8)
    ui->title->setText( "F_{L}^{t}(x), " + fPDF->PDFname());
  else if (index == 9)
    ui->title->setText( "F_{L}^{p}(x), " + fPDF->PDFname());
  else if (index == 10)
    ui->title->setText( "F_{3}^{l}(x), " + fPDF->PDFname());
  else if (index == 11)
    ui->title->setText( "F_{3}^{c}(x), " + fPDF->PDFname());
  else if (index == 12)
    ui->title->setText( "F_{3}^{b}(x), " + fPDF->PDFname());
  else if (index == 13)
    ui->title->setText( "F_{3}^{t}(x), " + fPDF->PDFname());
  else if (index == 14)
    ui->title->setText( "F_{3}^{p}(x), " + fPDF->PDFname());
  else if (index == 15)
    ui->title->setText( "#sigma^{l}(x), " + fPDF->PDFname());
  else if (index == 16)
    ui->title->setText( "#sigma^{c}(x), " + fPDF->PDFname());
  else if (index == 17)
    ui->title->setText( "#sigma^{b}(x), " + fPDF->PDFname());
  else if (index == 18)
    ui->title->setText( "#sigma^{t}(x), " + fPDF->PDFname());
  else if (index == 19)
    ui->title->setText( "#sigma^{p}(x), " + fPDF->PDFname());

}

disthread::disthread(QObject *parent, QString filename):
  QThread(parent),
  fp((PlotDIS*)parent),
  fFileName(filename),
  fIsTerminated(false)
{
}

disthread::~disthread()
{
}

void disthread::run()
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
  double xmin = lharanges[6][0];
  double xmax = lharanges[6][1];
  double ymin = 0, ymax = 0;

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
  fp->fPDF->InitPDFset2(Qi,Qf);
  const double y  = fp->ui->y->text().toDouble();
  const int pto = fp->ui->ptord->currentIndex();

  string proc;
  if (fp->ui->process->currentIndex() == 0)
    proc = "EM";
  else if (fp->ui->process->currentIndex() == 1)
    proc = "NC";
  else
    proc = "CC";

  string scheme;
  if (fp->ui->scheme->currentIndex() != 3)
    scheme = fp->ui->scheme->currentText().toStdString();

  string target;
  if (fp->ui->target->currentIndex() == 0)
    target = "PROTON";
  else if (fp->ui->target->currentIndex() == 1)
    target = "NEUTRON";
  else if (fp->ui->target->currentIndex() == 2)
    target = "ISOSCALAR";
  else
    target = "IRON";

  string project;
  if (fp->ui->projectile->currentIndex() == 0)
    project = "ELECTRON";
  else if (fp->ui->projectile->currentIndex() == 1)
    project = "POSITRON";
  else if (fp->ui->projectile->currentIndex() == 2)
    project = "NEUTRINO";
  else
    project = "ANTINEUTRINO";

  TLegend *leg =  new TLegend(0.603448,0.673729,0.981322,0.883475);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  int legindex = 0;

  TMultiGraph *mg = new TMultiGraph();

  double *x = new double[N];
  for (int i = 0; i < N; i++)
    {
      if (fp->ui->logx->isChecked()) x[i] = exp(log(xmin)+i*(log(xmax)-log(xmin))/N);
      else x[i] = xmin+i*(xmax-xmin)/N;
    }

  int nset = 1;
  vector<string> vnfs;
  if (fp->ui->scheme->currentIndex() == 3)
    {
      nset = 3;
      vnfs.push_back("FONLL");
      vnfs.push_back("FFNS");
      vnfs.push_back("ZMVN");
    }
  else
    vnfs.push_back(scheme);

  for (int set = 0; set < nset; set++)
    {
      scheme = vnfs[set];

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

      for (int i=0; i < N; i++)
        {

          if(fIsTerminated){ fIsTerminated = false; delete[] x; return; }

          emit progress(i*100/N);

          double res = 0, err = 0;
          double F2[5],F3[5],FL[5],sigma[5];
          double F2err[5],F3err[5],FLerr[5],sigmaerr[5];

          fp->fPDF->DIS(x[i],Qi,Qf,y,proc,scheme,pto,target,project,
                        F2,F3,FL,sigma,F2err,F3err,FLerr,sigmaerr);

          if (fp->ui->output->currentIndex() == 0) {
              res = F2[0]; err = F2err[0]; }
          if (fp->ui->output->currentIndex() == 1) {
            res = F2[1]; err = F2err[1];}
          if (fp->ui->output->currentIndex() == 2) {
            res = F2[2]; err = F2err[2];}
          if (fp->ui->output->currentIndex() == 3) {
            res = F2[3]; err = F2err[3];}
          if (fp->ui->output->currentIndex() == 4) {
            res = F2[4]; err = F2err[4];}

          if (fp->ui->output->currentIndex() == 5) {
              res = F3[0]; err = F3err[0];}
          if (fp->ui->output->currentIndex() == 6) {
            res = F3[1]; err = F3err[1];}
          if (fp->ui->output->currentIndex() == 7) {
            res = F3[2]; err = F3err[2];}
          if (fp->ui->output->currentIndex() == 8) {
            res = F3[3]; err = F3err[3];}
          if (fp->ui->output->currentIndex() == 9) {
            res = F3[4]; err = F3err[4];}

          if (fp->ui->output->currentIndex() == 10) {
            res = FL[0]; err = FLerr[0];}
          if (fp->ui->output->currentIndex() == 11) {
            res = FL[1]; err = FLerr[1];}
          if (fp->ui->output->currentIndex() == 12) {
            res = FL[2]; err = FLerr[2];}
          if (fp->ui->output->currentIndex() == 13) {
            res = FL[3]; err = FLerr[3];}
          if (fp->ui->output->currentIndex() == 14) {
            res = FL[4]; err = FLerr[4];}

          if (fp->ui->output->currentIndex() == 15) {
            res = sigma[0]; err = sigmaerr[0];}
          if (fp->ui->output->currentIndex() == 16) {
            res = sigma[1]; err = sigmaerr[1];}
          if (fp->ui->output->currentIndex() == 17) {
            res = sigma[2]; err = sigmaerr[2];}
          if (fp->ui->output->currentIndex() == 18) {
            res = sigma[3]; err = sigmaerr[3];}
          if (fp->ui->output->currentIndex() == 19) {
            res = sigma[4]; err = sigmaerr[4];}

          g->SetPoint(i,x[i],res);
          gcv->SetPoint(i, x[i], res);

          if (fp->ui->stddev->isChecked() && fp->fPDF->numberPDF() > 1) {
              g->SetPointError(i, 0, err);
              gup->SetPoint(i, x[i], res+err);
              gdn->SetPoint(i, x[i], res-err);
            }
          else
            g->SetPointError(i, 0, 0);
        }

      mg->Add(g,"le3");

      if (set == 0)
        mg->Add(gcv,"l");

      if (fp->ui->stddev->isChecked())
        {
          mg->Add(gup,"l");
          mg->Add(gdn,"l");
        }

      leg->AddEntry(g,TString(scheme.c_str()) + ", " + TString(fp->ui->ptord->currentText().toStdString()),"fl"); legindex++;
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
    mg->GetYaxis()->SetRangeUser(ymin,ymax);

  leg->AddEntry("",TString("Q = " + QString::number(fp->ui->Qf->text().toDouble(),'g',3).toStdString() + " GeV"),"");
  leg->AddEntry("",TString("Target: " + fp->ui->target->currentText().toStdString()),"");
  leg->AddEntry("",TString("Projec: " + fp->ui->projectile->currentText().toStdString()),"");
  if (legindex < 2)
    for (int i = legindex; i < 2; i++)
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

void disthread::SaveCanvas(QString filename)
{
  fC->SaveAs(filename.toStdString().c_str());
}

void disthread::stop()
{
  fIsTerminated = true;
}
