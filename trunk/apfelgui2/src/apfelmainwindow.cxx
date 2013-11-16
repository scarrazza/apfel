#include "apfelmainwindow.h"
#include "ui_apfelmainwindow.h"
#include "pdfdialog.h"
#include "plotmembers.h"
#include "plotall.h"
#include "plotcomparison.h"
#include "plotlumi.h"

#include <QInputDialog>
#include <QDir>
#include <QString>
#include <QMessageBox>
#include <QSettings>
#include <QDebug>
#include <QListView>
#include <QListWidgetItem>
#include <QDesktopWidget>

#include "APFEL/APFEL.h"

APFELMainWindow::APFELMainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::APFELMainWindow),
  fPDFs(),
  fLHAPDFpath("")
{
  ui->setupUi(this);

  QRect frect = frameGeometry();
  frect.moveCenter(QDesktopWidget().availableGeometry().center());
  move(frect.topLeft());

  connect(ui->tablePDFs, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(showDialog(QListWidgetItem*)));

  // Set Status bar text
  statusBar()->showMessage(QString("APFEL ") +
                           APFEL::GetVersion() +
                           QString(": V. Bertone, S. Carrazza and J. Rojo (arXiv:1310.1394)"));

  QFile input("apfelconfig.ini");
  if (!input.exists())
    {
      bool ok;
      fLHAPDFpath = QInputDialog::getText(this, tr("Welcome to APFEL GUI setup"),
                                          tr("Please set the LHAPDF sets directory where the *.LHgrid files are stored:"),
                                          QLineEdit::Normal,"/usr/local/lhapdf/share/lhapdf/PDFsets", &ok);
      if (ok && !fLHAPDFpath.isEmpty() && QDir(fLHAPDFpath).exists() == true)
        {
          // write to config file
          QSettings settings("apfelconfig.ini",QSettings::IniFormat);
          settings.beginGroup("LHAPDF");
          settings.setValue("LHGRIDPATH",fLHAPDFpath);
          settings.endGroup();
        }
      else
        {
          QMessageBox::warning(this,"Warning","Please set a valid path name!");
          exit(-1);
        }
    }
  else
    {
      // Read path from file
      QSettings settings("apfelconfig.ini",QSettings::IniFormat);
      settings.beginGroup("LHAPDF");
      fLHAPDFpath = settings.value("LHGRIDPATH").toString();
      settings.endGroup();
    }
}

APFELMainWindow::~APFELMainWindow()
{
  delete ui;

  for (size_t i = 0; i < fPDFs.size(); i++)
    if (fPDFs[i]) delete fPDFs[i];
  fPDFs.clear();
}

void APFELMainWindow::on_buttonAddPDF_clicked()
{
  PDFDialog *newPDF = new PDFDialog(this);

  if (newPDF->isAccept() == true)
    {
      fPDFs.push_back(newPDF);

      // add new pdf to table
      QListWidgetItem *item = new QListWidgetItem;
      item->setData( Qt::DisplayRole, newPDF->PDFname());
      item->setData( Qt::CheckStateRole, Qt::Checked );
      ui->tablePDFs->addItem(item);
    }
}

void APFELMainWindow::on_actionPreferences_triggered()
{
  bool ok;
  fLHAPDFpath = QInputDialog::getText(this, tr("Welcome to APFEL GUI setup"),
                                      tr("Please set the LHAPDF sets directory where the *.LHgrid files are stored:"),
                                      QLineEdit::Normal,fLHAPDFpath, &ok);
  if (ok)
    {
      if (!fLHAPDFpath.isEmpty() && QDir(fLHAPDFpath).exists() == true)
        {
          // write to config file
          QSettings settings("apfelconfig.ini",QSettings::IniFormat);
          settings.beginGroup("LHAPDF");
          settings.setValue("LHGRIDPATH",fLHAPDFpath);
          settings.endGroup();
        }
      else
        {
          QMessageBox::warning(this,"Warning","Please set a valid path name!");
        }
    }
}

void APFELMainWindow::on_buttonDelPDF_clicked()
{
  int id = ui->tablePDFs->currentRow();
  if (id >= 0)
    {
      delete ui->tablePDFs->takeItem(id);
      fPDFs.erase(fPDFs.begin()+id);
    }
}

void APFELMainWindow::showDialog(QListWidgetItem *)
{
  fPDFs[ui->tablePDFs->currentRow()]->exec();
  ui->tablePDFs->item(ui->tablePDFs->currentRow())->setData( Qt::DisplayRole, fPDFs[ui->tablePDFs->currentRow()]->PDFname());
}

void APFELMainWindow::on_pushButton_clicked()
{
  if (ui->tablePDFs->count() != 0)
    {
      int i1 = 0, i2 = 0;
      for (int i = 0; i < ui->tablePDFs->count(); i++)
        {
          i1 += ui->tablePDFs->item(i)->checkState()/2;
          if (i1 == 1) i2 = i;
        }

      if (i1 == 1)
        {
          PlotMembers *A = new PlotMembers(0,fPDFs[i2]);
          A->show();
        }
      else
        QMessageBox::information(this,"Information","Please select only one PDF set before continue.");
    }
  else
    QMessageBox::information(this,"Information","Please add a PDF set before continue.");
}

void APFELMainWindow::closeEvent(QCloseEvent*)
{
    qApp->quit();
}

void APFELMainWindow::on_pushButton_2_clicked()
{
    qApp->quit();
}

void APFELMainWindow::on_pushButton_5_clicked()
{
  if (ui->tablePDFs->count() != 0)
    {
      int i1 = 0, i2 = 0;
      for (int i = 0; i < ui->tablePDFs->count(); i++)
        {
          i1 += ui->tablePDFs->item(i)->checkState()/2;
          if (i1 == 1) i2 = i;
        }

      if (i1 == 1)
        {
          PlotAll *A = new PlotAll(0,fPDFs[i2]);
          A->show();
        }
      else
        QMessageBox::information(this,"Information","Please select only one PDF set before continue.");
    }
  else
    QMessageBox::information(this,"Information","Please add a PDF set before continue.");
}

void APFELMainWindow::on_pushButton_6_clicked()
{
  if (ui->tablePDFs->count() != 0)
    {
      std::vector<PDFDialog*> pdfs;
      for (int i = 0; i < ui->tablePDFs->count(); i++)
        {
          if (ui->tablePDFs->item(i)->checkState()/2 == 1)
            pdfs.push_back(fPDFs[i]);
        }

      if (pdfs.size() > 0)
        {
          if (pdfs.size() > 40)
            {
              QMessageBox::information(this,"Information","Too many PDFs selected, remove sets before continue.");
            }
          else
            {
              PlotComparison *A = new PlotComparison(0,pdfs);
              A->show();
            }
        }
      else
        QMessageBox::information(this,"Information","Please select at least one PDF set before continue.");
    }
  else
    QMessageBox::information(this,"Information","Please add a PDF set before continue.");
}

void APFELMainWindow::on_pushButton_7_clicked()
{
  if (ui->tablePDFs->count() != 0)
    {
      std::vector<PDFDialog*> pdfs;
      for (int i = 0; i < ui->tablePDFs->count(); i++)
        {
          if (ui->tablePDFs->item(i)->checkState()/2 == 1)
            pdfs.push_back(fPDFs[i]);
        }

      if (pdfs.size() > 0)
        {
          if (pdfs.size() > 40)
            {
              QMessageBox::information(this,"Information","Too many PDFs selected, remove sets before continue.");
            }
          else
            {
              PlotLumi *A = new PlotLumi(0,pdfs);
              A->show();
            }
        }
      else
        QMessageBox::information(this,"Information","Please select at least one PDF set before continue.");
    }
  else
    QMessageBox::information(this,"Information","Please add a PDF set before continue.");
}
