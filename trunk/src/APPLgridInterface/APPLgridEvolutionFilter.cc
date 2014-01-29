// APPLgridEvolutionFilter.cc
// Add to a pre-existing APPLgrid root file the PDF evolution provided by APFEL
//

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

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

  // Initial scale
  double Q0 = sqrt(2);
  
  // Init applgrid (forget about FastNLO for the moment)
  const igrid* igrid;
  const appl::grid *g = NULL;

  g = new appl::grid(argv[1]);

  string pdflbl = g->getGenpdf();
  int nsubproc  = g->subProcesses();
  int nbin      = g->Nobs();
  int no        = g->nloops();
  bool sym      = g->isSymmetric();

  // Call PDF subprocess generator
  appl_pdf *genpdf = appl_pdf::getpdf(pdflbl);

  cout << "  " << endl;
  cout << "APPLgridEvolutionFilter() Info:" << endl;
  //cout << "   Transform: " << g->getTransform() << endl;
  //cout << "   GenPDF:    " << pdflbl << endl;
  cout << "   SubProc:   " << nsubproc << endl;
  cout << "   Nbins:     " << nbin << endl;
  cout << "   Nloops:    " << no << endl;
  if(sym) cout << "   The process is symmetric" << endl;
  else                 cout << "   The process is NOT symmetric" << endl;
  cout << "  " << endl;

  for (int ibin=0; ibin<nbin; ibin++) {
    cout << "APPLgridEvolutionFilter() Info: Observable counter = " << ibin << endl;
    cout << "  " << endl;

    for (int o=0; o<=no; o++) {
      cout << "APPLgridEvolutionFilter() Info: Perturbative order = " << o << endl;
      cout << "  " << endl;

      int pto = o;
      igrid = g->weightgrid(pto,ibin);

      if (igrid==NULL) cerr << "APPLgridEvolutionFilter() Error: grid not found"<<endl;

      int ntau = igrid->Ntau();
      int nx1  = igrid->Ny1();
      int nx2  = igrid->Ny2();

      cout << "APPLgridEvolutionFilter() Info: Number of Q2 grid points = " << ntau << endl;
      //cout << "   Number of Q2 grid points = " << ntau << endl;
      //cout << "   Number of x1 grid points = " << nx1 << endl;
      //cout << "   Number of x2 grid points = " << nx2 << endl;
      cout << "  " << endl;

      for (int itau=0; itau<ntau; itau++) {
	double Q = sqrt(igrid->fQ2(igrid->gettau(itau)));
	//Q0 = Q;
	cout << "APPLgridEvolutionFilter() Info: Evolution between " << Q0 << " GeV and " << Q << " GeV" << endl;

	double delta = 0,deltap = 0;
	double a,ap;
	double eps = 1e-8;

	// Grid 1
	int n1b;
	double *M1 = NULL;

	int iy1 = nx1;
	vector<double> xext1v;

	for (int ix1=0; ix1<nx1; ix1++) {
	  iy1--; // inverse order
	  xext1v.push_back (igrid->fx(igrid->gety1(iy1)));
	  //cout << ix1 << "  " << xext1v[ix1] << endl;

	  // Parameter of the grid (see eq. (1) of arXiv:0911.2985)
	  if(ix1 > 0) {
	    deltap = delta;
	    ap     = a;

	    delta = igrid->gety1(iy1) - igrid->gety1(iy1+1);
	    a     = ( delta - log( xext1v[ix1-1] / xext1v[ix1] ) ) / ( xext1v[ix1-1] - xext1v[ix1] );

	    // Check that "delta" and "a" are constant all over the grid
	    if(ix1 > 1) {
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

	if((xext1v.size()-nx1) != 0) {
	  cout << "APPLgridEvolutionFilter() Error: Mismatch in size for xext1v, expected "<< nx1 << ", found " << xext1v.size() << endl;
	  exit(-10);
	}

	// If the Upper limit of the grid if different from 1, 
	// compute all the grid points up to one using "delta"
	// and "a" coomputed in the loop above.
	if( abs(xext1v[nx1-1]-1) > eps ) {
	  double yp = igrid->gety1(0);
	  cout << "APPLgridEvolutionFilter() Info: Extending grid 1 ..." << endl;
	  for(int ix1=0; xext1v[xext1v.size()-1]<1; ix1++) {
	    yp += delta;
	    xext1v.push_back (igrid->fx(yp));
	  }
	  // Force the last point to be equal to one
	  xext1v[xext1v.size()-1] = 1;
	}
	n1b = xext1v.size()-1;

	// Call APFEL and compute the evolution operator M1
	cout << "APPLgridEvolutionFilter() Info: Computing evolution operator on grid 1 ..." << endl;
	M1 = new double[14*14*(n1b+1)*(n1b+1)];
	APFEL::ExternalEvolutionOperator(Q0,Q,n1b,&xext1v[0],M1);

	// Grid 2
	// if the x-space grids are not equal for the PDFs ...
	// INFO: Usually APPLgrid says that they are not equal while they are.
	int n2b;
	double *M2 = NULL;
	if(!sym) {
	  int iy2 = nx2;
	  vector<double> xext2v;

	  for (int ix2=0; ix2<nx2; ix2++) {
	    iy2--; // inverse order
	    xext2v.push_back (igrid->fx(igrid->gety2(iy2)));
	    //cout << ix2 << "  " << xext2v[ix2] << endl;

	    // Parameter of the grid (see eq. (1) of arXiv:0911.2985)
	    if(ix2 > 0) {
	      deltap = delta;
	      ap     = a;

	      delta = igrid->gety2(iy2) - igrid->gety2(iy2+1);
	      a     = ( delta - log( xext2v[ix2-1] / xext2v[ix2] ) ) / ( xext2v[ix2-1] - xext2v[ix2] );

	      // Check that "delta" and "a" are constant all over the grid
	      if(ix2 > 1) {
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

	  if((xext2v.size()-nx2) != 0) {
	    cout << "APPLgridEvolutionFilter() Error: Mismatch in size for xext2v, expected "<< nx2 << ", found " << xext2v.size() << endl;
	    exit(-10);
	  }

	  // If the Upper limit of the grid if different from 1, 
	  // compute all the grid points up to one using "delta"
	  // and "a" coomputed in the loop above.
	  if( abs(xext2v[nx2-1]-1) > eps ) {
	    double yp = igrid->gety2(0);
	    cout << "APPLgridEvolutionFilter() Info: Extending grid 2 ..." << endl;
	    for(int ix2=0; xext2v[xext2v.size()-1]<1; ix2++) {
	      yp += delta;
	      xext2v.push_back (igrid->fx(yp));
	    }
	    // Force the last point to be equal to one
	    xext2v[xext2v.size()-1] = 1;
	  }
	  n2b = xext2v.size()-1;

	  // Call APFEL and compute the evolution operator M2
	  cout << "APPLgridEvolutionFilter() Info: Computing evolution operator on grid 2 ..." << endl;
	  M2 = new double[14*14*(n2b+1)*(n2b+1)];
	  APFEL::ExternalEvolutionOperator(Q0,Q,n2b,&xext2v[0],M2);
        }
	// ... otherwise copy "n1b" and "M1" into "n2b" and "M2"
	else {
	  n2b = n1b;
	  M2 = new double[14*14*(n2b+1)*(n2b+1)];
	  copy(M1, M1+14*14*(n1b+1)*(n1b+1), M2);
	}

	// Extract weights from the APPLgrid root file (in the inverse order)
	double ***W = new double**[nx1];
	for (int rho=0; rho<nx1; rho++) {
	  W[rho] = new double*[nx2];
	  for (int sigma=0; sigma<nx2; sigma++) {
	    W[rho][sigma] = new double[nsubproc];
	    for (int ip=0; ip<nsubproc; ip++) {
	      W[rho][sigma][ip] = (*(const SparseMatrix3d*) igrid->weightgrid(ip))(itau,nx1-1-rho,nx2-1-sigma);
	    }
	  }
	}

	// Put the first evolution operator in a suitable matrix form
	double ****M1a = new double***[n1b+1];
	for(int alpha=0; alpha<n1b+1; alpha++) {
	  M1a[alpha] = new double**[n1b+1];
	  for(int rho=0; rho<n1b+1; rho++) {
	    M1a[alpha][rho] = new double*[13];
	    for(int i=0; i<13; i++) {
	      M1a[alpha][rho][i] = new double[13];
	      for(int k=0; k<13; k++) {
		int m1 = ( k + 1 ) + 14 * ( ( i + 1 ) + 14 * ( rho   + ( n1b + 1 ) * alpha ) );
		M1a[alpha][rho][i][k] = M1[m1];
	      }
	    }
	  }
	}
	delete[] M1;

	// Put the second evolution operator in a suitable matrix form
	double ****M2a = new double***[n2b+1];
	for(int beta=0; beta<n2b+1; beta++) {
	  M2a[beta] = new double**[n2b+1];
	  for(int sigma=0; sigma<n2b+1; sigma++) {
	    M2a[beta][sigma] = new double*[13];
	    for(int j=0; j<13; j++) {
	      M2a[beta][sigma][j] = new double[13];
	      for(int l=0; l<13; l++) {
		int m2 = ( l + 1 ) + 14 * ( ( j + 1 ) + 14 * ( sigma   + ( n2b + 1 ) * beta ) );
		M2a[beta][sigma][j][l] = M2[m2];
	      }
	    }
	  }
	}
	delete[] M2;

	// Combine APPLgrid weights "W" with the evolution operators "M1a" and "M2a" and produce "Wt"
	double *H = new double[nsubproc];
	double ****Wt = new double***[n1b+1];
	for(int alpha=0; alpha<n1b+1; alpha++) {
	  double perc = (double) 100 * alpha / n1b;
	  cout << "APPLgridEvolutionFilter() Info: Combining APPLgrid weights with evolution operators " << setprecision(3) << setw(4) << perc << " %" << endl;
	  if(alpha != n1b) cout << "\x1b[A";
	  Wt[alpha] = new double**[n2b+1];
	  for(int beta=0; beta<n2b+1; beta++) {
	    Wt[alpha][beta] = new double*[13];
	    for(int i=0; i<13; i++) {
	      Wt[alpha][beta][i] = new double[13];
	      for(int j=0; j<13; j++) {
		Wt[alpha][beta][i][j] = 0;




		for(int rho=0; rho<nx1; rho++) {
		  for(int sigma=0; sigma<nx2; sigma++) {



		    double *f1 = M1a[alpha][rho][i];
		    double *f2 = M2a[beta][sigma][j];

		    genpdf->evaluate(f1,f2,H);

		    for (int ip=0; ip<nsubproc; ip++) {
		      Wt[alpha][beta][i][j] += W[rho][sigma][ip] * H[ip];
		    }



		  }
		}



	      }
	    }
	  }
	}
	delete[] W;
	delete[] H;
	delete[] M1a;
	delete[] M2a;




	cout << "  " << endl;
      }
    }
  }
  delete g;

  cout << "APPLgridEvolutionFilter(): Evolution included in " << argv[1] << endl;
  cout << "  " << endl;
  return 0;
}
