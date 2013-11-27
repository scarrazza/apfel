#include "apfelmainwindow.h"
#include <QApplication>
#include <QLocale>

int main(int argc, char *argv[])
{
  QLocale::setDefault(QLocale::c());
  QApplication a(argc, argv);
  APFELMainWindow w;
  w.show();
  
  return a.exec();
}
