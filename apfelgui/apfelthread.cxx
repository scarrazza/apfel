#include <QMessageBox>
#include <QAction>
#include <QString>
#include <QGraphicsPixmapItem>
#include <QDesktopWidget>
#include <QApplication>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QProgressDialog>
#include <QDebug>

#include "apfelthread.h"
#include "apfelmainwindow.h"
#include "ui_apfelmainwindow.h"

#include "APFEL/APFEL.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TImage.h"
#include "TLatex.h"
#include <cmath>

apfelthread::apfelthread(QObject *parent, Ui::MainWindow *fui,QProgressDialog *fd, int fmod):
  QThread(parent),
  ui(fui),
  d(fd),
  C(0),
  mod(fmod)
{
}

apfelthread::~apfelthread()
{

}

void apfelthread::run()
{

  // Initialize apfel
  APFEL::SetQLimits(0e0,1e4);
  APFEL::SetPerturbativeOrder(ui->comboBox_2->currentIndex());

  if (mod == 1)
    {
      APFEL::SetNumberOfGrids(3);
      APFEL::SetGridParameters(1,150,3,1e-9);
      APFEL::SetGridParameters(2,60,5,1e-1);
      APFEL::SetGridParameters(3,20,5,8e-1);
    }

  if (ui->comboBox_3->currentIndex() == 0)
    APFEL::SetVFNS();
  else
    APFEL::SetFFNS(ui->spinBox->value());

  APFEL::SetTheory(ui->comboBox->currentText().toStdString());

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

  d->setValue(3);

  APFEL::InitializeAPFEL();

  if (mod == 0)
    {

      d->setValue(6);

      APFEL::EvolveAPFEL(ui->lineEdit_2->text().toDouble(), ui->lineEdit_9->text().toDouble());

      ui->lineEdit_3->setText(QString::number(APFEL::AlphaQCD(ui->lineEdit_2->text().toDouble()),'f',10));
      ui->lineEdit_4->setText(QString::number(APFEL::AlphaQED(ui->lineEdit_9->text().toDouble()),'f',10));

      // sum rules
      double momsr = 0.0;
      for (int i = -6; i < 7; i++)
         momsr += APFEL::NPDF(i,2);

      momsr += APFEL::Ngamma(2);
      double uvsr = APFEL::NPDF(2,1) - APFEL::NPDF(-2,1);
      double dvsr = APFEL::NPDF(1,1) - APFEL::NPDF(-1,1);
      double svsr = APFEL::NPDF(3,1) - APFEL::NPDF(-3,1);

      ui->lineEdit_5->setText(QString::number(momsr,'e',10));
      ui->lineEdit_6->setText(QString::number(uvsr,'e',10));
      ui->lineEdit_7->setText(QString::number(dvsr,'e',10));
      ui->lineEdit_8->setText(QString::number(svsr,'e',10));

      // Plot
      C = new TCanvas("c","",625,466);
      C->SetLogx();
      C->SetLogy();
      C->SetTickx();
      C->SetTicky();

      const int N = 100;
      TGraph *g = new TGraph(N);
      TGraph *u = new TGraph(N);
      TGraph *d = new TGraph(N);
      TGraph *s = new TGraph(N);
      TGraph *c = new TGraph(N);
      TGraph *b = new TGraph(N);
      TGraph *p = new TGraph(N);

      TGraph *ub = new TGraph(N);
      TGraph *db = new TGraph(N);
      TGraph *sb = new TGraph(N);
      TGraph *cb = new TGraph(N);
      TGraph *bb = new TGraph(N);

      double logMin = log(1.001e-5);
      double logMax = log(1.0);
      double Delta = (logMax-logMin)/N;

      for (int i = 0; i < N; i++)
        {
          double x = exp(logMin + i*Delta);
          g->SetPoint(i, x, APFEL::xPDF(0,x));
          u->SetPoint(i, x, APFEL::xPDF(2,x));
          d->SetPoint(i, x, APFEL::xPDF(1,x));
          s->SetPoint(i, x, APFEL::xPDF(3,x));
          c->SetPoint(i, x, APFEL::xPDF(4,x));
          b->SetPoint(i, x, APFEL::xPDF(5,x));
          p->SetPoint(i, x, APFEL::xgamma(x));
          ub->SetPoint(i, x, APFEL::xPDF(-2,x));
          db->SetPoint(i, x, APFEL::xPDF(-1,x));
          sb->SetPoint(i, x, APFEL::xPDF(-3,x));
          cb->SetPoint(i, x, APFEL::xPDF(-4,x));
          bb->SetPoint(i, x, APFEL::xPDF(-5,x));
        }

      g->SetLineWidth(2);
      u->SetLineWidth(2);
      d->SetLineWidth(2);
      s->SetLineWidth(2);
      c->SetLineWidth(2);
      b->SetLineWidth(2);
      p->SetLineWidth(2);
      ub->SetLineWidth(2);
      db->SetLineWidth(2);
      sb->SetLineWidth(2);
      cb->SetLineWidth(2);
      bb->SetLineWidth(2);

      ub->SetLineStyle(2);
      db->SetLineStyle(2);
      sb->SetLineStyle(2);
      cb->SetLineStyle(2);
      bb->SetLineStyle(2);

      g->SetLineColor(kRed);
      u->SetLineColor(kBlue);
      d->SetLineColor(kGreen);
      s->SetLineColor(kCyan);
      c->SetLineColor(kMagenta);
      b->SetLineColor(kViolet);
      p->SetLineColor(kBlack);
      ub->SetLineColor(kBlue+2);
      db->SetLineColor(kGreen+2);
      sb->SetLineColor(kCyan+2);
      cb->SetLineColor(kMagenta+2);
      bb->SetLineColor(kViolet+2);

      TMultiGraph *mg = new TMultiGraph();
      mg->SetTitle( "APFEL: " + TString(ui->comboBox->currentText().toStdString() + " evolution @ " + ui->comboBox_2->currentText().toStdString()
                            + ", " + ui->comboBox_3->currentText().toStdString() + ", Q = " + QString::number(ui->lineEdit_9->text().toDouble(),'g',1).toStdString() + " GeV") );

      mg->Add(p,"l");
      mg->Add(g,"l");
      mg->Add(u,"l");
      mg->Add(d,"l");
      mg->Add(s,"l");
      mg->Add(c,"l");
      mg->Add(b,"l");
      mg->Add(ub,"l");
      mg->Add(db,"l");
      mg->Add(sb,"l");
      mg->Add(cb,"l");
      mg->Add(bb,"l");

      mg->SetMinimum(1e-5);

      mg->Draw("a");
      mg->GetXaxis()->SetTitle("x");
      mg->GetXaxis()->CenterTitle(true);
      mg->GetYaxis()->SetTitle("xPDF(x,Q)");
      mg->GetYaxis()->SetTitleOffset(1.2);
      mg->GetYaxis()->CenterTitle(true);
      mg->GetXaxis()->SetLimits(1e-5,1.0);

      TLegend *leg = new TLegend(0.130872,0.132114,0.614094,0.53252);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      leg->AddEntry(g,"xg(x,Q)","l");
      leg->AddEntry(u,"xu(x,Q)","l");
      leg->AddEntry(d,"xd(x,Q)","l");
      leg->AddEntry(s,"xs(x,Q)","l");
      leg->AddEntry(c,"xc(x,Q)","l");
      leg->AddEntry(b,"xb(x,Q)","l");
      leg->AddEntry(p,"x#gamma(x,Q)","l");

      leg->AddEntry(ub,"x#bar{u}(x,Q)","l");
      leg->AddEntry(db,"x#bar{d}(x,Q)","l");
      leg->AddEntry(sb,"x#bar{s}(x,Q)","l");
      leg->AddEntry(cb,"x#bar{c}(x,Q)","l");
      leg->AddEntry(bb,"x#bar{b}(x,Q)","l");

      leg->Draw("same");

      TLatex l; //l.SetTextAlign(12);
      l.SetTextSize(0.02);
      l.SetTextAngle(90);
      l.SetNDC();
      l.SetTextFont(72);
      l.SetTextColor(kBlack);
      l.DrawLatex(0.95,0.20,"Generated by APFEL: V.Bertone, S.Carrazza, J.Rojo (arXiv:1310.1394)");

      C->SaveAs("apfelplot.svg");
    }
  else
    {
      // luminosity
      double S = pow(ui->lineEdit_17->text().toDouble(),2.0);
      double Q0 = ui->lineEdit_2->text().toDouble();
      double Qmin = ui->lineEdit_18->text().toDouble();
      double Qmax = ui->lineEdit_19->text().toDouble();
      int NMX  = ui->spinBox_5->value();

      TGraph *g = new TGraph(NMX-1);
      TGraph *p = new TGraph(NMX-1);
      TGraph *q = new TGraph(NMX-1);
      TGraph *qg = new TGraph(NMX-1);
      TGraph *qq = new TGraph(NMX-1);

      for (int i=1; i < NMX; i++)
        {
          double Q = Qmin * pow(Qmax/Qmin,(double)(i-1)/(double)(NMX-1));
          APFEL::EvolveAPFEL(Q0,Q);

          d->setValue(d->value()+1);
          g->SetPoint(i-1, Q, APFEL::LUMI(0,0,S));

          d->setValue(d->value()+1);
          p->SetPoint(i-1,Q, APFEL::LUMI(7,7,S));

          d->setValue(d->value()+1);
          double result = 0;
          for (int j = 1; j <= 6; j++)
            result += APFEL::LUMI(j,-j,S);
          q->SetPoint(i-1,Q,result);

          d->setValue(d->value()+1);
          double xsinglet = 0.0;
          for (int j = 1; j <= 6; j++)
            xsinglet += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
          qg->SetPoint(i-1,Q,xsinglet);

          d->setValue(d->value()+1);
          double result2 = 0;
          for (int j = 1; j <= 6; j++)
            for (int z = j; z <= 6; z++)
               result2 += APFEL::LUMI(j,z,S);
          qq->SetPoint(i-1,Q,result2);
        }

      g->SetLineWidth(2);
      p->SetLineWidth(2);
      q->SetLineWidth(2);
      qg->SetLineWidth(2);
      qq->SetLineWidth(2);

      g->SetLineColor(kRed);
      p->SetLineColor(kBlack);
      q->SetLineColor(kBlue);
      qg->SetLineColor(kGreen);
      qq->SetLineColor(kViolet);

      C = new TCanvas("c2","",625,466);
      C->SetLogx();
      C->SetLogy();
      C->SetTickx();
      C->SetTicky();

      TMultiGraph *mgr = new TMultiGraph();
      mgr->SetTitle( "APFEL: " + TString(ui->comboBox->currentText().toStdString() + " luminosities @ " + ui->comboBox_2->currentText().toStdString()
                            + ", " + ui->comboBox_3->currentText().toStdString() + ", #sqrt{S} = " + QString::number(ui->lineEdit_17->text().toDouble(),'g',1).toStdString() + " GeV") );

      mgr->Add(g,"l");
      mgr->Add(p,"l");
      mgr->Add(q,"l");
      mgr->Add(qg,"l");
      mgr->Add(qq,"l");

      mgr->SetMinimum(1e-7);

      mgr->Draw("a");
      mgr->GetXaxis()->SetTitle("M_{x}");
      mgr->GetXaxis()->CenterTitle(true);
      mgr->GetYaxis()->SetTitle("Luminosity");
      mgr->GetYaxis()->SetTitleOffset(1.2);
      mgr->GetYaxis()->CenterTitle(true);
      mgr->GetXaxis()->SetLimits(Qmin,Qmax);

      TLegend *leg = new TLegend(0.586151,0.630137,0.954911,0.876712);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);

      leg->AddEntry(g,"gluon-gluon","l");
      leg->AddEntry(p,"#gamma-#gamma","l");
      leg->AddEntry(q,"quark-antiquark","l");
      leg->AddEntry(qg,"quark-gluon","l");
      leg->AddEntry(qq,"quark-quark","l");

      leg->Draw("same");

      TLatex l; //l.SetTextAlign(12);
      l.SetTextSize(0.02);
      l.SetTextAngle(90);
      l.SetNDC();
      l.SetTextFont(72);
      l.SetTextColor(kBlack);
      l.DrawLatex(0.95,0.20,"Generated by APFEL: V.Bertone, S.Carrazza, J.Rojo (arXiv:1310.1394)");

      C->SaveAs("apfelplot2.svg");

      ui->lineEdit_3->setText(QString::number(APFEL::AlphaQCD(ui->lineEdit_2->text().toDouble()),'f',10));
      ui->lineEdit_4->setText(QString::number(APFEL::AlphaQED(ui->lineEdit_9->text().toDouble()),'f',10));

      // sum rules
      double momsr = 0.0;
      for (int i = -6; i < 7; i++)
         momsr += APFEL::NPDF(i,2);

      momsr += APFEL::Ngamma(2);
      double uvsr = APFEL::NPDF(2,1) - APFEL::NPDF(-2,1);
      double dvsr = APFEL::NPDF(1,1) - APFEL::NPDF(-1,1);
      double svsr = APFEL::NPDF(3,1) - APFEL::NPDF(-3,1);

      ui->lineEdit_5->setText(QString::number(momsr,'e',10));
      ui->lineEdit_6->setText(QString::number(uvsr,'e',10));
      ui->lineEdit_7->setText(QString::number(dvsr,'e',10));
      ui->lineEdit_8->setText(QString::number(svsr,'e',10));
    }


}

void apfelthread::savecanvas(string file)
{
  C->SaveAs(file.c_str());
}
