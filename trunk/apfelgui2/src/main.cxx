#include "apfelmainwindow.h"
#include <QApplication>
#include <QLocale>

int main(int argc, char *argv[])
{
  QLocale loc = QLocale::system(); // current locale
  loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
  QLocale::setDefault(loc); // set as default
 
  QApplication a(argc, argv);
  APFELMainWindow w;
  w.show();
  
  return a.exec();
}
