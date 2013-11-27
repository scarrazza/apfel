#include "apfelmainwindow.h"
#include <QApplication>
#include <QLocale>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  APFELMainWindow w;
  w.show();
  
  return a.exec();
}
