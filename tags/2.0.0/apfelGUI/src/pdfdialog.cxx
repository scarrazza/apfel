#include "pdfdialog.h"
#include "ui_pdfdialog.h"
#include "apfelmainwindow.h"
#include "lumiintegral.h"

#include <QDir>
#include <QProcess>
#include <QDebug>
#include <cmath>

#include "LHAPDF/LHAPDF.h"
#include "APFEL/APFEL.h"

PDFDialog::PDFDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::PDFDialog),
  fisAccept(false),
  fQi(sqrt(2.0)),
  fQf(sqrt(2.0))
{
  ui->setupUi(this);
  ui->groupBox->setEnabled(true);

  // Loading list of LHgrids
  QDir path( ((APFELMainWindow*) parent)->getLHAPDFpath());
  QStringList files = path.entryList(QDir::Files);
  ui->comboPDFset->addItems(files);

  exec();
}

PDFDialog::~PDFDialog()
{
  delete ui;
}

void PDFDialog::InitPDFset(double Q0, double Q)
{
  if (isLHAPDF())
    {
      LHAPDF::initPDFSet(PDFname().toStdString());
    }
  else
    {
      fQi = Q0;
      fQf = Q;
      // Initialize apfel
      APFEL::SetQLimits(0e0,1e4);
      
      APFEL::SetNumberOfGrids(3);
      APFEL::SetGridParameters(1,100,3,1e-7);
      APFEL::SetGridParameters(2,50,5,0.1);
      APFEL::SetGridParameters(3,20,5,8e-1);
      
      APFEL::SetPerturbativeOrder(ptord());

      if (scheme() == 0)
        APFEL::SetVFNS();
      else
        APFEL::SetFFNS(nf());
      
      APFEL::SetTheory(theory().toStdString());
      APFEL::SetAlphaQCDRef(alphas(),Qalphas());
      APFEL::SetAlphaQEDRef(alpha(),Qalpha());
      
      if (ui->comboBox_hqscheme->currentIndex() == 0)
        APFEL::SetPoleMasses(ui->lineEdit_mc->text().toDouble(),
                             ui->lineEdit_mb->text().toDouble(),
                             ui->lineEdit_mt->text().toDouble());
      else
        APFEL::SetMSbarMasses(ui->lineEdit_mc->text().toDouble(),
                              ui->lineEdit_mb->text().toDouble(),
                              ui->lineEdit_mt->text().toDouble());

      APFEL::SetMaxFlavourPDFs(ui->spinBox_maxpdf->value());
      APFEL::SetMaxFlavourAlpha(ui->spinBox_maxa->value());

      APFEL::SetRenFacRatio(ui->doubleSpinBox_murmuf->value());
      
      QString apset = "ToyLH";
      if (ui->comboPDFset->currentIndex() != 0)
        apset = PDFname();

      APFEL::SetPDFSet(apset.toStdString());

      APFEL::InitializeAPFEL();
    }
}

void PDFDialog::InitPDFset2(double Q0, double Q)
{
  if (!isLHAPDF())
    {
      fQi = Q0;
      fQf = Q;
      // Initialize apfel
      APFEL::SetQLimits(0e0,1e4);

      APFEL::SetNumberOfGrids(3);
      APFEL::SetGridParameters(1,100,3,1e-7);
      APFEL::SetGridParameters(2,50,5,0.1);
      APFEL::SetGridParameters(3,20,5,8e-1);

      APFEL::SetPerturbativeOrder(ptord());

      if (scheme() == 0)
        APFEL::SetVFNS();
      else
        APFEL::SetFFNS(nf());

      APFEL::SetTheory(theory().toStdString());
      APFEL::SetAlphaQCDRef(alphas(),Qalphas());
      APFEL::SetAlphaQEDRef(alpha(),Qalpha());

      if (ui->comboBox_hqscheme->currentIndex() == 0)
        APFEL::SetPoleMasses(ui->lineEdit_mc->text().toDouble(),
                             ui->lineEdit_mb->text().toDouble(),
                             ui->lineEdit_mt->text().toDouble());
      else
        APFEL::SetMSbarMasses(ui->lineEdit_mc->text().toDouble(),
                              ui->lineEdit_mb->text().toDouble(),
                              ui->lineEdit_mt->text().toDouble());

      APFEL::SetMaxFlavourPDFs(ui->spinBox_maxpdf->value());
      APFEL::SetMaxFlavourAlpha(ui->spinBox_maxa->value());

      APFEL::SetRenFacRatio(ui->doubleSpinBox_murmuf->value());

      QString apset = "ToyLH";
      if (ui->comboPDFset->currentIndex() != 0)
        apset = PDFname();

      APFEL::SetPDFSet(apset.toStdString());
    }
}

void PDFDialog::on_comboBox_theory_currentIndexChanged(int index)
{
  if (index == 0)
    {
      ui->lineEdit_a->setEnabled(false);
      ui->lineEdit_Qa->setEnabled(false);
      ui->lineEdit_as->setEnabled(true);
      ui->lineEdit_Qas->setEnabled(true);
    }
  else
    {
      ui->lineEdit_a->setEnabled(true);
      ui->lineEdit_Qa->setEnabled(true);

      if (index == 1)
        {
          ui->lineEdit_as->setEnabled(false);
          ui->lineEdit_Qas->setEnabled(false);
        }
      else
        {
          ui->lineEdit_as->setEnabled(true);
          ui->lineEdit_Qas->setEnabled(true);
        }
    }
}

void PDFDialog::on_comboBox_scheme_currentIndexChanged(int index)
{
  if (index == 0)
    {
      ui->comboBox_hqscheme->setEnabled(true);
      ui->lineEdit_mc->setEnabled(true);
      ui->lineEdit_mb->setEnabled(true);
      ui->lineEdit_mt->setEnabled(true);
      ui->spinBox_maxa->setEnabled(true);
      ui->spinBox_maxpdf->setEnabled(true);
    }
  else
    {
      ui->comboBox_hqscheme->setEnabled(false);
      ui->lineEdit_mc->setEnabled(false);
      ui->lineEdit_mb->setEnabled(false);
      ui->lineEdit_mt->setEnabled(false);
      ui->spinBox_maxa->setEnabled(false);
      ui->spinBox_maxpdf->setEnabled(false);
    }
}

void PDFDialog::on_dialDLAGP_valueChanged(int value)
{
  if (value == 0)
    ui->groupBox->setEnabled(false);
  else
    ui->groupBox->setEnabled(true);
}

void PDFDialog::on_buttonBox_accepted()
{
    fisAccept = true;
}

void PDFDialog::on_buttonBox_rejected()
{
    fisAccept = false;
}

QString PDFDialog::PDFname()
{
  return ui->comboPDFset->currentText();
}

void PDFDialog::on_comboPDFset_currentIndexChanged(int index)
{
  if (index == 0) {
    ui->dialDLAGP->setValue(1);
    ui->dialDLAGP->setEnabled(false);    
    ui->comboPDFerror->setEnabled(false);
    ui->comboPDFerror->setCurrentIndex(0);
    ui->member->setEnabled(false);
    ui->member->setValue(0);
  } else {
    ui->dialDLAGP->setValue(0);
    ui->dialDLAGP->setEnabled(true);
    ui->comboPDFerror->setEnabled(true);
    if (ui->comboPDFerror->currentIndex() == 0)
      {
        ui->member->setEnabled(true);
        ui->member->setValue(0);
      }
    else
      {
        ui->member->setEnabled(false);
        ui->member->setValue(0);
      }
  }
}

bool PDFDialog::isLHAPDF()
{
  if (ui->dialDLAGP->value() == 0)
    return true;
  else
    return false;
}

int PDFDialog::ptord()
{
  return ui->comboBox_ptord->currentIndex();
}

int PDFDialog::scheme()
{
  return ui->comboBox_scheme->currentIndex();
}

int PDFDialog::nf()
{
  return ui->comboBox_scheme->currentIndex()+2;
}

QString PDFDialog::theory()
{
  return ui->comboBox_theory->currentText();
}

double PDFDialog::alpha()
{
  return ui->lineEdit_a->text().toDouble();
}

double PDFDialog::Qalpha()
{
  return ui->lineEdit_Qa->text().toDouble();
}

double PDFDialog::alphas()
{
  return ui->lineEdit_as->text().toDouble();
}

double PDFDialog::Qalphas()
{
  return ui->lineEdit_Qas->text().toDouble();
}

int PDFDialog::numberPDF()
{
  if (isLHAPDF())
    {
      if (ui->comboPDFerror->currentIndex() != ER_NONE)
        return LHAPDF::numberPDF();
      else
        return 1;
    }
  else
    {
      if (ui->comboPDFerror->currentIndex() != ER_NONE)
        return LHAPDF::numberPDF();
      else
        return 1;
    }
}

void PDFDialog::initPDF(int i)
{
  if (isLHAPDF())
    LHAPDF::initPDF(i);
  else
    Evolve(i,fQi,fQf);
}

void PDFDialog::GetFlvrPDFCVErr(int N, double *x, double Q, int f,
                                double *xPDF, double *xPDFerr,double* upErr,double* dnErr)
{
  int Etype = ui->comboPDFerror->currentIndex();

  switch (Etype) {
    case ER_NONE:
      {
        initPDF(ui->member->value());
        for (int i = 0; i < N; i++) {
          xPDF[i] = GetFlvrPDF(x[i],Q,f);
          xPDFerr[i] = upErr[i] = dnErr[i] = 0;
        }
        break;
      }
    case ER_EIG:
    case ER_EIG90:
    case ER_SYMEIG:
      {
        initPDF(0);
        for (int i = 0; i < N; i++) {
          xPDF[i] = GetFlvrPDF(x[i],Q,f);
          upErr[i] = dnErr[i] = 0;
        }

        double **y = new double*[numberPDF()];
        for (int i = 0; i < numberPDF(); i++)
          {
            initPDF(i+1);
            y[i] = new double[N];
            for (int j = 0; j < N; j++)
              y[i][j] = GetFlvrPDF(x[j],Q,f);
          }

        for (int i = 0; i < N; i++)
          {
            if (Etype != ER_SYMEIG)
              {
                xPDFerr[i] = ComputeEigErr(numberPDF(),i,y);
                if (Etype == ER_EIG90)
                  xPDFerr[i] /= 1.64485;
              }
            else
              xPDFerr[i] = ComputeSymEigErr(numberPDF(),i,xPDF[i],y);
          }

        for (int i = 0; i < numberPDF(); i++)
          if (y[i]) delete[] y[i];
        delete[] y;

        break;
      }
    case ER_MC:
      {
        vector<vector<double> > yval;
        double **y = new double*[numberPDF()];

        for (int i = 0; i < numberPDF(); i++)
          {
            initPDF(i+1);
            y[i] = new double[N];
            for (int j = 0; j < N; j++)
              y[i][j] = GetFlvrPDF(x[j],Q,f);
          }

        for (int i = 0; i < N; i++)
          {
            vector<double> w;
            for (int j = 0; j < numberPDF(); j++)
              w.push_back(y[j][i]);
            yval.push_back(w);
          }

        int esc = numberPDF()*(1-0.68)/2;
        for (int i = 0; i < N; i++) {
          xPDF[i] = ComputeAVG(numberPDF(),i,y);
          xPDFerr[i] = ComputeStdDev(numberPDF(),i,y);
          sort(yval[i].begin(),yval[i].end());
          upErr[i] = yval[i][numberPDF()-esc-1];
          dnErr[i] = yval[i][esc];
        }

        for (int i = 0; i < numberPDF(); i++)
          if (y[i]) delete[] y[i];
        delete[] y;

        yval.clear();

        break;
      }
    }
}

double PDFDialog::GetFlvrPDFCV(double x, double Q, int f)
{
  double avg = 0;
  int Etype = ui->comboPDFerror->currentIndex();

  if (isLHAPDF())
    {
      switch (Etype) {
        case ER_NONE:
          {
            initPDF(ui->member->value());
            if (LHAPDF::hasPhoton() == true)
              avg = LHAPDF::xfxphoton(x,Q,f);
            else
              avg = LHAPDF::xfx(x,Q,f);
            break;
          }
        case ER_EIG:
        case ER_EIG90:
        case ER_SYMEIG:
          {
            initPDF(0);
            if (LHAPDF::hasPhoton() == true)
              avg = LHAPDF::xfxphoton(x,Q,f);
            else
              avg = LHAPDF::xfx(x,Q,f);
            break;
          }
        case ER_MC:
          {
            double *y = new double[numberPDF()];
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (LHAPDF::hasPhoton() == true)
                  y[i] = LHAPDF::xfxphoton(x,Q,f);
                else
                  y[i] = LHAPDF::xfx(x,Q,f);
              }
            avg = ComputeAVG(numberPDF(),y);
            delete[] y;
            break;
          }
        }
    }
  else
    {
      switch (Etype) {
        case ER_NONE:
          {
            initPDF(ui->member->value());
            if (f != 7)
              avg = APFEL::xPDF(f,x);
            else
              avg = APFEL::xgamma(x);
            break;
          }
        case ER_EIG:
        case ER_EIG90:
        case ER_SYMEIG:
          {
            initPDF(0);
            if (f != 7)
              avg = APFEL::xPDF(f,x);
            else
              avg = APFEL::xgamma(x);
            break;
          }
        case ER_MC:
          {
            double *y = new double[numberPDF()];
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (f != 7)
                  y[i] = APFEL::xPDF(f,x);
                else
                  y[i] = APFEL::xgamma(x);
              }
            avg = ComputeAVG(numberPDF(),y);
            delete[] y;
            break;
          }
        }
    }

  return avg;
}

double PDFDialog::GetFlvrError(double x, double Q, int f,double &uperr, double &dnerr)
{
  double err = 0;
  int Etype = ui->comboPDFerror->currentIndex();

  if (isLHAPDF())
    {
      switch (Etype) {
        case ER_NONE:
          {
            break;
          }
        case ER_MC:
          {
            double *y = new double[numberPDF()];
            vector<double> yval;
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (LHAPDF::hasPhoton() == true)
                  y[i] = LHAPDF::xfxphoton(x,Q,f);
                else
                  y[i] = LHAPDF::xfx(x,Q,f);
                yval.push_back(y[i]);
              }
            err = ComputeStdDev(numberPDF(),y);

            delete[] y;

            sort(yval.begin(), yval.end());
            int esc = numberPDF()*(1-0.68)/2;

            uperr = yval[numberPDF()-esc-1];
            dnerr = yval[esc];

            break;
          }
        case ER_EIG:
        case ER_EIG90:
          {
            double *y = new double[numberPDF()];
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (LHAPDF::hasPhoton() == true)
                  y[i] = LHAPDF::xfxphoton(x,Q,f);
                else
                  y[i] = LHAPDF::xfx(x,Q,f);
              }
            err = ComputeEigErr(numberPDF(),y);
            if (Etype == ER_EIG90) err /= 1.64485;

            delete[] y;

            break;
          }
        case ER_SYMEIG:
          {
            double *y = new double[numberPDF()];
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (LHAPDF::hasPhoton() == true)
                  y[i] = LHAPDF::xfxphoton(x,Q,f);
                else
                  y[i] = LHAPDF::xfx(x,Q,f);
              }

            err = ComputeSymEigErr(numberPDF(),GetFlvrPDFCV(x,Q,f),y);

            delete[] y;
            break;
          }
        }
    }
  else
    {
      switch (Etype) {
        case ER_NONE:
          {
            break;
          }
        case ER_MC:
          {
            double *y = new double[numberPDF()];
            vector<double> yval;
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (f != 7)
                  y[i] = APFEL::xPDF(f,x);
                else
                  y[i] = APFEL::xgamma(x);
                yval.push_back(y[i]);
              }
            err = ComputeStdDev(numberPDF(),y);

            delete[] y;

            sort(yval.begin(), yval.end());
            int esc = numberPDF()*(1-0.68)/2;

            uperr = yval[numberPDF()-esc-1];
            dnerr = yval[esc];

            break;
          }
        case ER_EIG:
        case ER_EIG90:
          {
            double *y = new double[numberPDF()];
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (f != 7)
                  y[i] = APFEL::xPDF(f,x);
                else
                  y[i] = APFEL::xgamma(x);
              }
            err = ComputeEigErr(numberPDF(),y);
            if (Etype == ER_EIG90) err /= 1.64485;

            delete[] y;

            break;
          }
        case ER_SYMEIG:
          {
            double *y = new double[numberPDF()];
            for (int i = 0; i < numberPDF(); i++)
              {
                initPDF(i+1);
                if (f != 7)
                  y[i] = APFEL::xPDF(f,x);
                else
                  y[i] = APFEL::xgamma(x);
              }

            err = ComputeSymEigErr(numberPDF(),GetFlvrPDFCV(x,Q,f),y);

            delete[] y;
            break;
          }
        }
    }

  return err;
}

double PDFDialog::GetFlvrPDF(double x, double Q, int f)
{
  double res = 0;  

  if (isLHAPDF())
    {
      if (LHAPDF::hasPhoton() == true)
        res = LHAPDF::xfxphoton(x,Q,f);
      else
        res = LHAPDF::xfx(x,Q,f);
    }
  else
    {
      if (f != 7)
        res = APFEL::xPDF(f,x);
      else
        res = APFEL::xgamma(x);
    }

  return res;
}

void PDFDialog::Evolve(int i,double Q0,double Q)
{
  APFEL::SetReplica(i);
  APFEL::EvolveAPFEL(Q0,Q);
}

double PDFDialog::getLum(double i, double S, std::string lumi, double eps, double& err)
{
  double res = 0;
  err = 0;
  int Etype = ui->comboPDFerror->currentIndex();

  if (isLHAPDF())
    {
      switch (Etype) {
        case ER_NONE:
          {
            initPDF(ui->member->value());
            LumiIntegral *l = new LumiIntegral(eps);
            res = l->getLum(i,S,lumi);
            delete l;
            break;
          }
        case ER_EIG:
        case ER_EIG90:
          {
            initPDF(0);
            LumiIntegral *l = new LumiIntegral(eps);
            res = l->getLum(i,S,lumi);
            delete l;

            double *ggflux = new double[numberPDF()];
            for (int r = 0; r < numberPDF(); r++)
              {
                initPDF(r+1);
                LumiIntegral *w = new LumiIntegral(eps);
                ggflux[r] = w->getLum(i,S,lumi);
                delete w;
              }

            err = ComputeEigErr(numberPDF(),ggflux);
            if (Etype == ER_EIG90) err /= 1.64485;

            delete[] ggflux;
            break;
          }
        case ER_SYMEIG:
          {
            initPDF(0);
            LumiIntegral *l = new LumiIntegral(eps);
            res = l->getLum(i,S,lumi);
            delete l;

            double *ggflux = new double[numberPDF()];
            for (int r = 0; r < numberPDF(); r++)
              {
                initPDF(r+1);
                LumiIntegral *w = new LumiIntegral(eps);
                ggflux[r] = w->getLum(i,S,lumi);
                delete w;
              }

            err = ComputeSymEigErr(numberPDF(),res,ggflux);

            delete[] ggflux;
            break;
          }
        case ER_MC:
          {
            double *ggflux = new double[numberPDF()];
            for (int r = 0; r < numberPDF(); r++)
              {
                initPDF(r+1);
                LumiIntegral *l = new LumiIntegral(eps);
                ggflux[r] = l->getLum(i,S,lumi);
                delete l;
              }
            res = ComputeAVG(numberPDF(),ggflux);
            err = ComputeStdDev(numberPDF(),ggflux);
            delete[] ggflux;
            break;
          }
        }
    }
  else
    {
      switch (Etype) {
        case ER_NONE:
          {
            APFEL::SetReplica(ui->member->value());
            APFEL::EvolveAPFEL(fQi,i);

            if (lumi == "GG")
              res = APFEL::LUMI(0,0,S);
            else if (lumi == "PP")
              res = APFEL::LUMI(7,7,S);
            else if (lumi == "QG")
              {
                for (int j = 1; j <= 6; j++)
                  res += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
              }
            else if (lumi == "QQ")
              {
                for (int j = 1; j <= 6; j++)
                  res += APFEL::LUMI(j,-j,S);
              }
            else if (lumi == "Q2")
              {
                for (int j = 1; j <= 6; j++)
                  for (int z = j; z <= 6; z++)
                     res += APFEL::LUMI(j,z,S);
              }
            else if (lumi == "BB")
              res = APFEL::LUMI(5,-5,S);
            else if (lumi == "CC")
              res = APFEL::LUMI(4,-4,S);
            else if (lumi == "BG")
              res = APFEL::LUMI(5,0,S);
            else if (lumi == "GC")
              res = APFEL::LUMI(4,0,S);
            else if (lumi == "PG")
              res = APFEL::LUMI(7,0,S);

            break;
          }
        case ER_EIG:
        case ER_EIG90:
          {
            APFEL::SetReplica(0);
            APFEL::EvolveAPFEL(fQi,i);

            if (lumi == "GG")
              res = APFEL::LUMI(0,0,S);
            else if (lumi == "PP")
              res = APFEL::LUMI(7,7,S);
            else if (lumi == "QG")
              {
                for (int j = 1; j <= 6; j++)
                  res += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
              }
            else if (lumi == "QQ")
              {
                for (int j = 1; j <= 6; j++)
                  res += APFEL::LUMI(j,-j,S);
              }
            else if (lumi == "Q2")
              {
                for (int j = 1; j <= 6; j++)
                  for (int z = j; z <= 6; z++)
                     res += APFEL::LUMI(j,z,S);
              }
            else if (lumi == "BB")
              res = APFEL::LUMI(5,-5,S);
            else if (lumi == "CC")
              res = APFEL::LUMI(4,-4,S);
            else if (lumi == "BG")
              res = APFEL::LUMI(5,0,S);
            else if (lumi == "GC")
              res = APFEL::LUMI(4,0,S);
            else if (lumi == "PG")
              res = APFEL::LUMI(7,0,S);


            double *ggflux = new double[numberPDF()];
            for (int r = 0; r < numberPDF(); r++)
              {
                APFEL::SetReplica(r+1);
                APFEL::EvolveAPFEL(fQi,i);

                if (lumi == "GG")
                  ggflux[r] = APFEL::LUMI(0,0,S);
                else if (lumi == "PP")
                  ggflux[r] = APFEL::LUMI(7,7,S);
                else if (lumi == "QG")
                  {
                    for (int j = 1; j <= 6; j++)
                      ggflux[r] += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
                  }
                else if (lumi == "QQ")
                  {
                    for (int j = 1; j <= 6; j++)
                      ggflux[r] += APFEL::LUMI(j,-j,S);
                  }
                else if (lumi == "Q2")
                  {
                    for (int j = 1; j <= 6; j++)
                      for (int z = j; z <= 6; z++)
                         ggflux[r] += APFEL::LUMI(j,z,S);
                  }
                else if (lumi == "BB")
                  ggflux[r] = APFEL::LUMI(5,-5,S);
                else if (lumi == "CC")
                  ggflux[r] = APFEL::LUMI(4,-4,S);
                else if (lumi == "BG")
                  ggflux[r] = APFEL::LUMI(5,0,S);
                else if (lumi == "GC")
                  ggflux[r] = APFEL::LUMI(4,0,S);
                else if (lumi == "PG")
                  ggflux[r] = APFEL::LUMI(7,0,S);
              }

            err = ComputeEigErr(numberPDF(),ggflux);
            if (Etype == ER_EIG90) err /= 1.64485;

            delete[] ggflux;
            break;
          }
        case ER_SYMEIG:
          {
            APFEL::SetReplica(0);
            APFEL::EvolveAPFEL(fQi,i);

            if (lumi == "GG")
              res = APFEL::LUMI(0,0,S);
            else if (lumi == "PP")
              res = APFEL::LUMI(7,7,S);
            else if (lumi == "QG")
              {
                for (int j = 1; j <= 6; j++)
                  res += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
              }
            else if (lumi == "QQ")
              {
                for (int j = 1; j <= 6; j++)
                  res += APFEL::LUMI(j,-j,S);
              }
            else if (lumi == "Q2")
              {
                for (int j = 1; j <= 6; j++)
                  for (int z = j; z <= 6; z++)
                     res += APFEL::LUMI(j,z,S);
              }
            else if (lumi == "BB")
              res = APFEL::LUMI(5,-5,S);
            else if (lumi == "CC")
              res = APFEL::LUMI(4,-4,S);
            else if (lumi == "BG")
              res = APFEL::LUMI(5,0,S);
            else if (lumi == "GC")
              res = APFEL::LUMI(4,0,S);
            else if (lumi == "PG")
              res = APFEL::LUMI(7,0,S);


            double *ggflux = new double[numberPDF()];
            for (int r = 0; r < numberPDF(); r++)
              {
                APFEL::SetReplica(r+1);
                APFEL::EvolveAPFEL(fQi,i);

                if (lumi == "GG")
                  ggflux[r] = APFEL::LUMI(0,0,S);
                else if (lumi == "PP")
                  ggflux[r] = APFEL::LUMI(7,7,S);
                else if (lumi == "QG")
                  {
                    for (int j = 1; j <= 6; j++)
                      ggflux[r] += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
                  }
                else if (lumi == "QQ")
                  {
                    for (int j = 1; j <= 6; j++)
                      ggflux[r] += APFEL::LUMI(j,-j,S);
                  }
                else if (lumi == "Q2")
                  {
                    for (int j = 1; j <= 6; j++)
                      for (int z = j; z <= 6; z++)
                         ggflux[r] += APFEL::LUMI(j,z,S);
                  }
                else if (lumi == "BB")
                  ggflux[r] = APFEL::LUMI(5,-5,S);
                else if (lumi == "CC")
                  ggflux[r] = APFEL::LUMI(4,-4,S);
                else if (lumi == "BG")
                  ggflux[r] = APFEL::LUMI(5,0,S);
                else if (lumi == "GC")
                  ggflux[r] = APFEL::LUMI(4,0,S);
                else if (lumi == "PG")
                  ggflux[r] = APFEL::LUMI(7,0,S);
              }

            err = ComputeSymEigErr(numberPDF(),res,ggflux);

            delete[] ggflux;
            break;
          }
        case ER_MC:
          {
            double *ggflux = new double[numberPDF()];
            for (int r = 0; r < numberPDF(); r++)
              {
                APFEL::SetReplica(r+1);
                APFEL::EvolveAPFEL(fQi,i);

                if (lumi == "GG")
                  ggflux[r] = APFEL::LUMI(0,0,S);
                else if (lumi == "PP")
                  ggflux[r] = APFEL::LUMI(7,7,S);
                else if (lumi == "QG")
                  {
                    for (int j = 1; j <= 6; j++)
                      ggflux[r] += APFEL::LUMI(j,0,S) + APFEL::LUMI(-j,0,S);
                  }
                else if (lumi == "QQ")
                  {
                    for (int j = 1; j <= 6; j++)
                      ggflux[r] += APFEL::LUMI(j,-j,S);
                  }
                else if (lumi == "Q2")
                  {
                    for (int j = 1; j <= 6; j++)
                      for (int z = j; z <= 6; z++)
                         ggflux[r] += APFEL::LUMI(j,z,S);
                  }
                else if (lumi == "BB")
                  ggflux[r] = APFEL::LUMI(5,-5,S);
                else if (lumi == "CC")
                  ggflux[r] = APFEL::LUMI(4,-4,S);
                else if (lumi == "BG")
                  ggflux[r] = APFEL::LUMI(5,0,S);
                else if (lumi == "GC")
                  ggflux[r] = APFEL::LUMI(4,0,S);
                else if (lumi == "PG")
                  ggflux[r] = APFEL::LUMI(7,0,S);
              }

            res = ComputeAVG(numberPDF(),ggflux);
            err = ComputeStdDev(numberPDF(),ggflux);

            delete[] ggflux;
            break;
          }
        }
    }

  return res;
}

int PDFDialog::GetErrorType()
{
  return ui->comboPDFerror->currentIndex();
}


void PDFDialog::on_comboPDFerror_currentIndexChanged(int index)
{
  if (index == 0)
    {
      ui->member->setEnabled(true);
      ui->member->setValue(0);
    }
  else
    {
      ui->member->setEnabled(false);
      ui->member->setValue(0);
    }
}

int PDFDialog::GetReplica()
{
  return ui->member->value();
}

void PDFDialog::DIS(double x,double qi,double qf,double y,
		    const std::string& proc,const std::string& scheme,
		    int pto, const std::string& target, const std::string& proj,
		    double *F2, double *F3, double *FL, double *sigma,
		    double *F2err, double *F3err, double *FLerr, double *sigmaerr)
{
  int Etype = ui->comboPDFerror->currentIndex();

  for (int i = 0; i < 5; i++){
      F2[i] = F3[i] = FL[i] = sigma[i] = 0.0;
      F2err[i] = F3err[i] = FLerr[i] = sigmaerr[i] = 0.0;
    }

  switch (Etype) {
    case ER_NONE:
      {
        int irep = ui->member->value();

        QString pdfset = "APFEL";
        if (ui->comboPDFset->currentIndex() != 0)
          pdfset = PDFname();

        APFEL::DIS_xsec(x,qi,qf,y,proc,scheme,pto,pdfset.toStdString(),irep,
                        target,proj,F2,F3,FL,sigma);
        break;
      }
    case ER_EIG:
    case ER_EIG90:
      {
        QString pdfset = "APFEL";
        if (ui->comboPDFset->currentIndex() != 0)
          pdfset = PDFname();

        APFEL::DIS_xsec(x,qi,qf,y,proc,scheme,pto,pdfset.toStdString(),0,
                        target,proj,F2,F3,FL,sigma);

        double **F2a = new double*[numberPDF()];
        double **F3a = new double*[numberPDF()];
        double **FLa = new double*[numberPDF()];
        double **sigmaa = new double*[numberPDF()];

        for (int r = 0; r < numberPDF(); r++)
          {
            F2a[r] = new double[5];
            F3a[r] = new double[5];
            FLa[r] = new double[5];
            sigmaa[r] = new double[5];
            APFEL::DIS_xsec(x,qi,qf,y,proc,scheme,pto,pdfset.toStdString(),r+1,
                            target,proj,F2a[r],F3a[r],FLa[r],sigmaa[r]);
          }

        for (int i = 0; i < 5; i++)
          {
            F2err[i] = ComputeEigErr(numberPDF(),i,F2a);
            F3err[i] = ComputeEigErr(numberPDF(),i,F3a);
            FLerr[i] = ComputeEigErr(numberPDF(),i,FLa);
            sigmaerr[i] = ComputeEigErr(numberPDF(),i,sigmaa);

            if (Etype == ER_EIG90)
              {
                F2err[i] /= 1.64485;
                F3err[i] /= 1.64485;
                FLerr[i] /= 1.64485;
                sigmaerr[i] /= 1.64485;
              }
          }

        for (int i = 0; i < numberPDF(); i++){
            if (F2a[i]) delete[] F2a[i];
            if (F3a[i]) delete[] F3a[i];
            if (FLa[i]) delete[] FLa[i];
            if (sigmaa[i]) delete[] sigmaa[i];
          }
        delete[] F2a;
        delete[] F3a;
        delete[] FLa;
        delete[] sigmaa;

        break;
      }
    case ER_SYMEIG:
      {
        QString pdfset = "APFEL";
        if (ui->comboPDFset->currentIndex() != 0)
          pdfset = PDFname();

        APFEL::DIS_xsec(x,qi,qf,y,proc,scheme,pto,pdfset.toStdString(),0,
                        target,proj,F2,F3,FL,sigma);

        double **F2a = new double*[numberPDF()];
        double **F3a = new double*[numberPDF()];
        double **FLa = new double*[numberPDF()];
        double **sigmaa = new double*[numberPDF()];

        for (int r = 0; r < numberPDF(); r++)
          {
            F2a[r] = new double[5];
            F3a[r] = new double[5];
            FLa[r] = new double[5];
            sigmaa[r] = new double[5];
            APFEL::DIS_xsec(x,qi,qf,y,proc,scheme,pto,pdfset.toStdString(),r+1,
                            target,proj,F2a[r],F3a[r],FLa[r],sigmaa[r]);
          }

        for (int i = 0; i < 5; i++)
          {
            F2err[i] = ComputeSymEigErr(numberPDF(),i,F2[i],F2a);
            F3err[i] = ComputeSymEigErr(numberPDF(),i,F3[i],F3a);
            FLerr[i] = ComputeSymEigErr(numberPDF(),i,FL[i],FLa);
            sigmaerr[i] = ComputeSymEigErr(numberPDF(),i,sigma[i],sigmaa);
          }

        for (int i = 0; i < numberPDF(); i++){
            if (F2a[i]) delete[] F2a[i];
            if (F3a[i]) delete[] F3a[i];
            if (FLa[i]) delete[] FLa[i];
            if (sigmaa[i]) delete[] sigmaa[i];
          }
        delete[] F2a;
        delete[] F3a;
        delete[] FLa;
        delete[] sigmaa;

        break;
      }
    case ER_MC:
      {
        QString pdfset = "APFEL";
        if (ui->comboPDFset->currentIndex() != 0)
          pdfset = PDFname();

        double **F2a = new double*[numberPDF()];
        double **F3a = new double*[numberPDF()];
        double **FLa = new double*[numberPDF()];
        double **sigmaa = new double*[numberPDF()];

        for (int r = 0; r < numberPDF(); r++)
          {
            F2a[r] = new double[5];
            F3a[r] = new double[5];
            FLa[r] = new double[5];
            sigmaa[r] = new double[5];
            APFEL::DIS_xsec(x,qi,qf,y,proc,scheme,pto,pdfset.toStdString(),r+1,
                            target,proj,F2a[r],F3a[r],FLa[r],sigmaa[r]);
          }

        for (int i = 0; i < 5; i++)
          {
            F2[i] = ComputeAVG(numberPDF(),i,F2a);
            F3[i] = ComputeAVG(numberPDF(),i,F3a);
            FL[i] = ComputeAVG(numberPDF(),i,FLa);
            sigma[i] = ComputeAVG(numberPDF(),i,sigmaa);

            F2err[i] = ComputeStdDev(numberPDF(),i,F2a);
            F3err[i] = ComputeStdDev(numberPDF(),i,F3a);
            FLerr[i] = ComputeStdDev(numberPDF(),i,FLa);
            sigmaerr[i] = ComputeStdDev(numberPDF(),i,sigmaa);
          }

        for (int i = 0; i < numberPDF(); i++){
            if (F2a[i]) delete[] F2a[i];
            if (F3a[i]) delete[] F3a[i];
            if (FLa[i]) delete[] FLa[i];
            if (sigmaa[i]) delete[] sigmaa[i];
          }
        delete[] F2a;
        delete[] F3a;
        delete[] FLa;
        delete[] sigmaa;

        break;
      }

    }
}
