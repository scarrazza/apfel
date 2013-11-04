#include "pdfdialog.h"
#include "ui_pdfdialog.h"
#include "apfelmainwindow.h"

#include <QDir>
#include <QProcess>
#include <QDebug>
#include <cmath>

#include "LHAPDF/LHAPDF.h"
#include "APFEL/APFEL.h"

PDFDialog::PDFDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::PDFDialog),
  fisAccept(false)
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

void PDFDialog::InitPDFset()
{
  if (isLHAPDF())
    {
      LHAPDF::initPDFSetByName(PDFname().toStdString());
    }
  else
    {
      // Initialize apfel
      APFEL::SetQLimits(0e0,1e4);
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

      qDebug() << apset;
      APFEL::SetPDFSet(apset.toStdString());

      //APFEL::SetReplica(0);

      APFEL::InitializeAPFEL();
      APFEL::EvolveAPFEL(sqrt(2.0),sqrt(2.0));
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
    ui->comboPDFerror->setEnabled(false);
  } else {
    ui->dialDLAGP->setValue(0);
    ui->comboPDFerror->setEnabled(true);
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
      return LHAPDF::numberPDF();
    }
  else
    {
      if (ui->comboPDFset->currentIndex() != 0)
        return LHAPDF::numberPDF();
      else
        return 1;
    }
}
