#ifndef APFELMAINWINDOW_H
#define APFELMAINWINDOW_H

#include <QMainWindow>
#include <vector>

class PDFDialog;
class QString;
class QListWidgetItem;

namespace Ui {
  class APFELMainWindow;
}

class APFELMainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  explicit APFELMainWindow(QWidget *parent = 0);
  ~APFELMainWindow();
  QString getLHAPDFpath() { return fLHAPDFpath; }
  
private slots:
  void on_buttonAddPDF_clicked();

  void on_actionPreferences_triggered();

  void on_buttonDelPDF_clicked();

  void showDialog(QListWidgetItem*);

  void on_pushButton_clicked();

  void on_pushButton_2_clicked();

  void on_pushButton_5_clicked();

protected:
  void  closeEvent(QCloseEvent*);

private:
  Ui::APFELMainWindow *ui;
  std::vector<PDFDialog*> fPDFs;
  QString fLHAPDFpath;
};

#endif // APFELMAINWINDOW_H
