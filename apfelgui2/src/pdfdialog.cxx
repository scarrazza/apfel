#include "pdfdialog.h"
#include "ui_pdfdialog.h"
#include "apfelmainwindow.h"

#include <QDir>
#include <QProcess>
#include <QDebug>

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
