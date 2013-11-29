#ifndef PDFDIALOG_H
#define PDFDIALOG_H

#include <QDialog>

class QString;
class PDFset;

namespace Ui {
  class PDFDialog;
}


class PDFDialog : public QDialog
{
  Q_OBJECT
  
public:
  explicit PDFDialog(QWidget *parent = 0);
  ~PDFDialog();
  void InitPDFset(double,double);
  void    initPDF(int);
  QString PDFname();  
  QString theory();
  bool    isAccept() { return fisAccept; }
  bool    isLHAPDF();  
  int     ptord();
  int     scheme();
  int     nf();
  int     numberPDF();
  double  alpha();
  double  Qalpha();
  double  alphas();
  double  Qalphas();
  double  GetFlvrPDFCV(double,double,int);
  double  GetFlvrPDF(double,double,int);
  double  getLum(double,double,std::string,double);
  double  GetFlvrError(double x, double Q, int f,double& uperr, double& dnerr);
  int     GetErrorType();

private slots:
  void on_comboBox_theory_currentIndexChanged(int index);
  void on_comboBox_scheme_currentIndexChanged(int index);
  void on_dialDLAGP_valueChanged(int value);
  void on_buttonBox_accepted();
  void on_buttonBox_rejected();
  void on_comboPDFset_currentIndexChanged(int index);

private:
  Ui::PDFDialog *ui;
  bool fisAccept;
  double fQi;
  double fQf;
  void Evolve(int,double,double);
};

#endif // PDFDIALOG_H
