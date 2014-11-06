#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "stdio.h"
#include "stdlib.h"

using namespace std;

// lhapdf routines
#include "LHAPDF/LHAPDF.h"
extern "C" void   evolvepdf_(const double& , const double& , double* ); 
extern "C" double alphaspdf_(const double& Q);

int main(int argc, char* argv[]) {

  // Check that the input is correct
  if (argc!=3) {
    cout << "Invalid Parameters:" << endl;
    cout << "Syntax: ./TestConvolution <APPLgrid file 1> <APPLgrid file 2>" << endl;
    exit(1);
  }

  const string pdfname = "MSTW2008nlo68cl.LHgrid";
  int iset = 0; // Central set
  LHAPDF::initPDFSet(pdfname,iset);
  LHAPDF::initPDF(iset);

  appl::grid *g1  = new appl::grid(argv[1]);
  const appl::igrid* ig1;
  appl::grid *g2 = new appl::grid(argv[2]);
  const appl::igrid* ig2;

  for(int nloop=0; nloop<2; nloop++) {
    if(nloop == 0) cout << "Computing LO" << endl;
    if(nloop == 1) cout << "Computing NLO" << endl;
    cout << " " << endl;

    vector<double> dist1 = g1->vconvolute(evolvepdf_, alphaspdf_, -nloop);
    vector<double> dist2 = g2->vconvolute(evolvepdf_, alphaspdf_, -nloop);
    //cout << g1->subProcesses() << "   " << g2->subProcesses() << endl;

    if(dist1.size() != dist2.size()) {
      cout << "Error: Mismatch in size" << endl;
      cout << "size 1 = " << dist1.size() << endl;
      cout << "size 2 = " << dist2.size() << endl;
      //      exit(-10);
    }

    cout << "Original    Evolved     RelDiff[%]" << endl;
    for (unsigned int i=0; i<dist1.size(); i++) {
      ig1 = g1->weightgrid(nloop,i);
      ig2 = g2->weightgrid(nloop,i);
      //cout << ig1->fQ2(ig1->gettau(0)) << "   " << ig2->fQ2(ig2->gettau(0)) << endl;

      double obs1 = dist1.at(i);
      double obs2 = dist2.at(i);
      double reldiff = 100 * ( obs2 - obs1 ) / obs1;
      //double ref1 = g1->getReference()->GetBinContent(i+1);
      //double ref2 = g2->getReference()->GetBinContent(i+1);
      cout << obs1 << "  " << obs2 << "  " << reldiff << endl;
    }
    cout << " " << endl;
    /*
    cout << "*************************************" << endl;
    double sum1 = 0; 
    for(int iproc=0; iproc<g1->subProcesses(); iproc++) {
      vector<double> dist1 = g1->vconvolute_subproc(iproc, evolvepdf_, alphaspdf_,  -nloop);

      for (unsigned int i=0; i<1; i++) {
	double obs1 = dist1.at(i);
	sum1 += obs1;
	cout << iproc << "  " << obs1 << endl;
      }
    }
    cout << "sum1 = " << sum1 << endl;
    double sum2 = 0;
    for(int iproc=0; iproc<g2->subProcesses(); iproc++) {
      vector<double> dist2 = g2->vconvolute_subproc(iproc, evolvepdf_, alphaspdf_, -nloop);

      for (unsigned int i=0; i<1; i++) {
	double obs2 = dist2.at(i);
	sum2 += obs2;
	cout << iproc << "  " << obs2 << endl;
      }
    }
    cout << "sum2 = " << sum2 << endl;
    */
  }
  return 0;
}


