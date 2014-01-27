// APPLgridEvolutionFilter.cc
// Add to a pre-existing APPLgrid root file the PDF evolution provided by APFEL
//

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

// APFEL
#include "APFEL/APFEL.h"

// APPLgrid
#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"
#include "appl_grid/SparseMatrix3d.h"

#include "APPLgridInterface.h"

using namespace std;

int main(int argc, char* argv[]) {

  if (argc!=2) {
    cout << "Invalid Parameters:" << endl;
    cout << "Syntax: ./appl_kinvar <APPLgrid file>" << endl;
    exit(1);
  }
  
  // Init applgrid (forget about FastNLO for the moment)
  const igrid* igrid;
  const appl::grid *g = NULL;

  g = new appl::grid(argv[1]);

  cout << "  " << endl;
  cout << "APPLgridEvolutionFilter() information:" << endl;
  cout << "   Transform: " << g->getTransform() << endl;
  cout << "   GenPDF:    " << g->getGenpdf() << endl;
  cout << "   SubProc:   " << g->subProcesses() << endl;
  cout << "   Nbins:     " << g->Nobs() << endl;
  cout << "   Nloops:    " << g->nloops() << endl;
  if(g->isSymmetric()) cout << "   The process is symmetric" << endl;
  else                 cout << "   The process is NOT symmetric" << endl;
  cout << "  " << endl;

  // x and Q2 grid
  double Q0 = sqrt(2); // Initial scale

  for (int ibin=0; ibin<g->Nobs(); ibin++) {
    cout << "APPLgridEvolutionFilter() information:" << endl;
    cout << "   Observable counter = " << ibin << endl;
    cout << "  " << endl;

    for (int o=0; o<=g->nloops(); o++) {
      cout << "APPLgridEvolutionFilter() information:" << endl;
      cout << "   Perturbative order = " << o << endl;
      cout << "  " << endl;

      int pto = o;
      igrid = g->weightgrid(pto,ibin);

      if (igrid==NULL) cerr << "APPLgridEvolutionFilter() Error: grid not found"<<endl;

      cout << "APPLgridEvolutionFilter() information:" << endl;
      cout << "   Number of Q2 grid points = " << igrid->Ntau() << endl;
      cout << "   Number of x1 grid points = " << igrid->Ny1() << endl;
      cout << "   Number of x2 grid points = " << igrid->Ny2() << endl;
      cout << "  " << endl;

      for (int itau=0; itau<igrid->Ntau(); itau++) {

	double Q = sqrt(igrid->fQ2(igrid->gettau(itau)));
	cout << "Q = " << Q << " GeV" << endl;

	double delta,deltap;
	double a,ap;
	double eps = 1e-8;

	// Grid 1
	int    n1 = igrid->Ny1()-1;
	double *xext1 = new double[igrid->Ny1()];
	double *M1 = NULL;
	int ix1 = n1 + 1;

	for (int iy1=0; iy1<igrid->Ny1(); iy1++) {

	  ix1--; // inverse order
	  xext1[ix1] = igrid->fx(igrid->gety1(iy1));

	  // Parameter of the grid (see eq. (1) of arXiv:0911.2985)
	  if(iy1 > 0) {
	    deltap = delta;
	    ap     = a;

	    delta = igrid->gety1(iy1) - igrid->gety1(iy1-1);
	    a     = ( delta - log( xext1[ix1+1] / xext1[ix1] ) ) / ( xext1[ix1+1] - xext1[ix1] );

	    // Check that "delta" and "a" are constant all over the grid
	    if(iy1 > 1) {
	      if( abs(delta - deltap) > eps ) {
		cout << "APPLgridEvolutionFilter() Error: The step is not constant " << endl;
		exit(-10);
	      }

	      if( abs(a - ap) > eps ) {
		cout << "APPLgridEvolutionFilter() Error: The parameter 'a' is not constant " << endl;
		exit(-10);
	      }

	    }
	  }
	}

	// If the Upper limit of the grid if different from 1, 
	// compute all the grid points up to one using "delta"
	// and "a" coomputed in the loop above.
	int n1b = n1;
	if( abs(xext1[n1]-1) > eps ) {
	  double yp = igrid->gety1(0);
	  cout << "APPLgridEvolutionFilter() Warning: Extending grid 1" << endl;
	  for(ix1=n1; xext1[n1b]<1; ix1++) {
	    yp -= delta;
	    n1b++;
	    xext1[n1b] = igrid->fx(yp);
	  }
	  // Force the last point to be equal to one
	  xext1[n1b] = 1;
	}

	//for(ix1=0; ix1<=n1b; ix1++) cout << ix1 << "  " << xext1[ix1] << endl;

	// Call APFEL and compute the evolution operator M1
	APFEL::ExternalEvolutionOperator(Q0,Q,n1b,xext1,M1);

	// Grid 2, if the x-space grids are not equal for the PDFs
	// INFO: Usually APPLgrid says that they are not equal while they are.
	if(!g->isSymmetric()) {
	  int    n2 = igrid->Ny2()-1;
	  double *xext2 = new double[igrid->Ny2()];
	  double *M2 = NULL;
	  int ix2 = n2 + 1;

	  for (int iy2=0; iy2<igrid->Ny2(); iy2++) {

	    ix2--; // inverse order
	    xext2[ix2] = igrid->fx(igrid->gety2(iy2));

	    // Parameter of the grid (see eq. (1) of arXiv:0911.2985)
	    if(iy2 > 0) {
	      deltap = delta;
	      ap     = a;

	      delta = igrid->gety2(iy2) - igrid->gety2(iy2-1);
	      a     = ( delta - log( xext2[ix2+1] / xext2[ix2] ) ) / ( xext2[ix2+1] - xext2[ix2] );

	      // Check that "delta" and "a" are constant all over the grid
	      if(iy2 > 1) {
		if( abs(delta - deltap) > eps ) {
		  cout << "APPLgridEvolutionFilter() Error: The step is not constant " << delta << "  " << deltap << endl;
		  exit(-10);
		}

		if( abs(a - ap) > eps ) {
		  cout << "APPLgridEvolutionFilter() Error: The parameter 'a' is not constant " << endl;
		  exit(-10);
		}

	      }
	    }
	  }

	  // If the Upper limit of the grid if different from 1, 
	  // compute all the grid points up to one using "delta"
	  // and "a" coomputed in the loop above.
	  int n2b = n2;
	  if( abs(xext2[n2]-1) > eps ) {
	    double yp = igrid->gety2(0);
	    cout << "APPLgridEvolutionFilter() Warning: Extending grid 2" << endl;
	    for(ix2=n2; xext2[n2b]<1; ix2++) {
	      yp -= delta;
	      n2b++;
	      xext2[n2b] = igrid->fx(yp);
	    }
	    // Force the last point to be equal to one
	    xext2[n2b] = 1;
	  }

	  //for(ix2=0; ix2<=n2b; ix2++) cout << ix2 << "  " << xext2[ix2] << endl;

	  // Call APFEL and compute the evolution operator M2
	  APFEL::ExternalEvolutionOperator(Q0,Q,n2b,xext2,M2);
	}
      }
      cout << "  " << endl;
    }
  }
  delete g;
  return 0;
}
