// APPLgridEvolutionFilter.cc
// Add to a pre-existing APPLgrid root file the PDF evolution provided by APFEL
//

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

// APFEL
#include "APFEL/APFEL.h"

// APPLgrid
#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"
#include "appl_grid/lumi_pdf.h"
#include "appl_grid/SparseMatrix3d.h"

#include "APPLgridInterface.h"

using namespace std;

int main(int argc, char* argv[]) {

  // Check that the input is correct
  if (argc!=2) {
    cout << "Invalid Parameters:" << endl;
    cout << "Syntax: ./appl_kinvar <APPLgrid file>" << endl;
    exit(1);
  }

  // Initial scale
  double Q0 = sqrt(2);
  double Q20 = Q0 * Q0;

  // Flavour map
  bool FlavMap[13][13]={{0},{0}};
  double Channels[13][13]={{0},{0}};
  string flavour[13]={"t~","b~","c~","s~","u~","d~","g ","d ","u ","s ","c ","b ","t "};

  // Input applgrid (forget about FastNLO for the moment)
  const igrid* igrid;
  const appl::grid *g = NULL;
  g = new appl::grid(argv[1]);

  // Parameters of the input grid
  string pdflbl = g->getGenpdf();
  int nsubproc  = g->subProcesses();
  int nbin      = g->Nobs();
  int no        = g->nloops();
  bool sym      = g->isSymmetric();

  // Call PDF subprocess generator
  appl_pdf *genpdf = appl_pdf::getpdf(pdflbl);

  // Display parameter
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

  // Declare the 6-dim array of weights and other stuff
  double ******Wt = new double*****[nbin];
  double **xmin = new double*[nbin];
  int **n1b = new int*[nbin];
  int **n2b = new int*[nbin];

  // Loop over the bins
  for(int ibin=0; ibin<nbin; ibin++) {
  //for(int ibin=0; ibin<1; ibin++) {
    cout << "APPLgridEvolutionFilter() Info: bin " << ibin+1 << " of " << nbin << endl;
    cout << "  " << endl;

    Wt[ibin]   = new double****[no];
    xmin[ibin] = new double[no];
    n1b[ibin]   = new int[no];
    n2b[ibin]   = new int[no];

    // Loop over the perturbative orders
    for(int o=0; o<=no; o++) {
    //for(int o=1; o<=1; o++) {
      if(o == 0) cout << "APPLgridEvolutionFilter() Info: Perturbative order: LO " << endl;
      if(o == 1) cout << "APPLgridEvolutionFilter() Info: Perturbative order: NLO " << endl;
      if(o == 2) cout << "APPLgridEvolutionFilter() Info: Perturbative order: NNLO " << endl;
      cout << "  " << endl;

      // Access internal grid for "pto" and "ibin"
      int pto = o;
      igrid = g->weightgrid(pto,ibin);
      if (igrid==NULL) cerr << "APPLgridEvolutionFilter() Error: grid not found"<<endl;

      // Parameter of the current internal grid
      int ntau = igrid->Ntau();
      int nx1  = igrid->Ny1();
      int nx2  = igrid->Ny2();
      //cout << "   Number of Q2 grid points = " << ntau << endl;
      //cout << "   Number of x1 grid points = " << nx1 << endl;
      //cout << "   Number of x2 grid points = " << nx2 << endl;
      //cout << "  " << endl;

      // Loop over the grid in Q2
      for(int itau=0; itau<ntau; itau++) {
      //for(int itau=1; itau<2; itau++) {
	cout << "APPLgridEvolutionFilter() Info: Q2 grid point " << itau+1 << " of " << ntau << endl;
	cout << "  " << endl;

	double Q = sqrt(igrid->fQ2(igrid->gettau(itau)));
	Q0  = Q;
	Q20 = Q * Q;
	cout << "APPLgridEvolutionFilter() Info: Evolution between " << Q0 << " GeV and " << Q << " GeV" << endl;

	double delta = 0,deltap = 0;
	double a,ap;
	double eps = 1e-8;

	// Generate the evolution table for the first grid in x
	double *M1 = NULL;
	vector<double> xext1v;
	int iy1 = nx1;

	for(int ix1=0; ix1<nx1; ix1++) {
	  iy1--; // inverse order
	  xext1v.push_back(igrid->fx(igrid->gety1(iy1)));

	  // Parameter of the grid (see eq. (1) of arXiv:0911.2985)
	  if(ix1 > 0) {
	    deltap = delta;
	    ap     = a;
	    delta = igrid->gety1(iy1) - igrid->gety1(iy1+1);
	    a     = ( delta - log( xext1v[ix1-1] / xext1v[ix1] ) ) / ( xext1v[ix1-1] - xext1v[ix1] );

	    // Check that "delta" and "a" are constant all over the grid
	    if(ix1 > 1) {
	      if(abs(delta - deltap) > eps) {
		cout << "APPLgridEvolutionFilter() Error: The step is not constant " << endl;
		exit(-10);
	      }
	      if(abs(a - ap) > eps) {
		cout << "APPLgridEvolutionFilter() Error: The parameter 'a' is not constant " << endl;
		exit(-10);
	      }
	    }
	  }
	}
	// Check that the vector "xext1v" is not empty
	if((xext1v.size()-nx1) != 0) {
	  cout << "APPLgridEvolutionFilter() Error: Mismatch in size for xext1v, expected "<< nx1 << ", found " << xext1v.size() << endl;
	  exit(-10);
	}

	// If the Upper limit of the grid if different from 1, 
	// compute all the grid points up to one using "delta"
	// and "a" coomputed in the loop above.
	if(abs(xext1v[nx1-1]-1) > eps) {
	  double yp = igrid->gety1(0);
	  cout << "APPLgridEvolutionFilter() Info: Extending grid 1 ..." << endl;
	  for(int ix1=0; xext1v[xext1v.size()-1]<1; ix1++) {
	    yp += delta;
	    xext1v.push_back(igrid->fx(yp));
	  }
	  // Force the last point to be equal to one
	  xext1v[xext1v.size()-1] = 1;
	}
	n1b[ibin][pto] = xext1v.size()-1;

	// Call APFEL and compute the evolution operator M1
	cout << "APPLgridEvolutionFilter() Info: Computing evolution operator on grid 1 ..." << endl;
	M1 = new double[14*14*(n1b[ibin][pto]+1)*(n1b[ibin][pto]+1)];
	APFEL::SetPerturbativeOrder(0);
	//APFEL::SetFFNS(3);
	APFEL::ExternalEvolutionOperator(Q0,Q,n1b[ibin][pto],&xext1v[0],M1);

	// Generate the evolution table for the second grid in x.
	// if the x-space grids are not equal for the PDFs ...
	// INFO: Usually APPLgrid says that they are not equal while they are.
	double *M2 = NULL;
	vector<double> xext2v;
	if(!sym) {
	  int iy2 = nx2;

	  for(int ix2=0; ix2<nx2; ix2++) {
	    iy2--; // inverse order
	    xext2v.push_back(igrid->fx(igrid->gety2(iy2)));
	    //cout << ix2 << "  " << xext2v[ix2] << endl;

	    // Parameter of the grid (see eq. (1) of arXiv:0911.2985)
	    if(ix2 > 0) {
	      deltap = delta;
	      ap     = a;
	      delta = igrid->gety2(iy2) - igrid->gety2(iy2+1);
	      a     = ( delta - log( xext2v[ix2-1] / xext2v[ix2] ) ) / ( xext2v[ix2-1] - xext2v[ix2] );

	      // Check that "delta" and "a" are constant all over the grid
	      if(ix2 > 1) {
		if(abs(delta - deltap) > eps) {
		  cout << "APPLgridEvolutionFilter() Error: The step is not constant " << endl;
		  exit(-10);
		}
		if(abs(a - ap) > eps) {
		  cout << "APPLgridEvolutionFilter() Error: The parameter 'a' is not constant " << endl;
		  exit(-10);
		}
	      }
	    }
	  }
	  // Check that the vector "xext2v" is not empty
	  if((xext2v.size()-nx2) != 0) {
	    cout << "APPLgridEvolutionFilter() Error: Mismatch in size for xext2v, expected "<< nx2 << ", found " << xext2v.size() << endl;
	    exit(-10);
	  }

	  // If the Upper limit of the grid if different from 1, 
	  // compute all the grid points up to one using "delta"
	  // and "a" coomputed in the loop above.
	  if(abs(xext2v[nx2-1]-1) > eps) {
	    double yp = igrid->gety2(0);
	    cout << "APPLgridEvolutionFilter() Info: Extending grid 2 ..." << endl;
	    for(int ix2=0; xext2v[xext2v.size()-1]<1; ix2++) {
	      yp += delta;
	      xext2v.push_back(igrid->fx(yp));
	    }
	    // Force the last point to be equal to one
	    xext2v[xext2v.size()-1] = 1;
	  }
	  n2b[ibin][pto] = xext2v.size()-1;

	  // Call APFEL and compute the evolution operator M2
	  cout << "APPLgridEvolutionFilter() Info: Computing evolution operator on grid 2 ..." << endl;
	  M2 = new double[14*14*(n2b[ibin][pto]+1)*(n2b[ibin][pto]+1)];
	  APFEL::ExternalEvolutionOperator(Q0,Q,n2b[ibin][pto],&xext2v[0],M2);
        }
	// ... otherwise copy "n1b", "xext1v" and "M1" into "n2b", "xext2v" and "M2"
	else {
	  n2b[ibin][pto] = n1b[ibin][pto];
	  xext2v = xext1v;
	  M2 = new double[14*14*(n2b[ibin][pto]+1)*(n2b[ibin][pto]+1)];
	  copy(M1, M1+14*14*(n1b[ibin][pto]+1)*(n1b[ibin][pto]+1), M2);
	}

	// Find the minimum value of x
	if(sym) xmin[ibin][pto] = xext1v[0];
	else    xmin[ibin][pto] = min(xext1v[0],xext2v[0]);

	// Extract weights from the APPLgrid root file (in the inverse order)
	double ***W = new double**[nx1];
	for(int rho=0; rho<nx1; rho++) {
	  W[rho] = new double*[nx2];
	  for(int sigma=0; sigma<nx2; sigma++) {
	    W[rho][sigma] = new double[nsubproc];
	    for(int ip=0; ip<nsubproc; ip++) {
	      W[rho][sigma][ip] = (*(const SparseMatrix3d*) igrid->weightgrid(ip))(itau,nx1-1-rho,nx2-1-sigma);
	    }
	  }
	}

	// Put the first evolution operator in a suitable matrix form
	double ****M1a = new double***[n1b[ibin][pto]+1];
	for(int alpha=0; alpha<n1b[ibin][pto]+1; alpha++) {
	  M1a[alpha] = new double**[n1b[ibin][pto]+1];
	  for(int rho=0; rho<n1b[ibin][pto]+1; rho++) {
	    M1a[alpha][rho] = new double*[13];
	    for(int i=0; i<13; i++) {
	      M1a[alpha][rho][i] = new double[13];
	      for(int k=0; k<13; k++) {
		int m1 = ( k + 1 ) + 14 * ( ( i + 1 ) + 14 * ( rho   + ( n1b[ibin][pto] + 1 ) * alpha ) );
		M1a[alpha][rho][i][k] = M1[m1];
	      }
	    }
	  }
	}
	delete[] M1;

	// Put the second evolution operator in a suitable matrix form
	double ****M2a = new double***[n2b[ibin][pto]+1];
	for(int beta=0; beta<n2b[ibin][pto]+1; beta++) {
	  M2a[beta] = new double**[n2b[ibin][pto]+1];
	  for(int sigma=0; sigma<n2b[ibin][pto]+1; sigma++) {
	    M2a[beta][sigma] = new double*[13];
	    for(int j=0; j<13; j++) {
	      M2a[beta][sigma][j] = new double[13];
	      for(int l=0; l<13; l++) {
		int m2 = ( l + 1 ) + 14 * ( ( j + 1 ) + 14 * ( sigma   + ( n2b[ibin][pto] + 1 ) * beta ) );
		M2a[beta][sigma][j][l] = M2[m2];
	      }
	    }
	    /*
	    //if(beta == sigma) {
	      cout << "beta = " << beta << " sigma = " << sigma << endl;
	      for(int j=0; j<13; j++) {
		for(int l=0; l<13; l++) {
		  cout << M2a[beta][sigma][l][j] << "  ";
		}
		cout << endl;
	      }
	      cout << "****************************************" << endl;
	    //}
	    */
	  }
	}
	delete[] M2;

	// Combine APPLgrid weights "W" with the evolution operators "M1a" and "M2a" and produce "Wt"
	double *H = new double[nsubproc];
	Wt[ibin][pto] = new double***[n1b[ibin][pto]+1];
	for(int alpha=0; alpha<n1b[ibin][pto]+1; alpha++) {
	  double perc = (double) 100 * alpha / n1b[ibin][pto];
	  cout << "APPLgridEvolutionFilter() Info: Combining APPLgrid weights with evolution operators: " << setprecision(3) << setw(4) << perc << " % done" << endl;
	  if(alpha != n1b[ibin][pto]) cout << "\x1b[A";
	  Wt[ibin][pto][alpha] = new double**[n2b[ibin][pto]+1];
	  for(int beta=0; beta<n2b[ibin][pto]+1; beta++) {
	    Wt[ibin][pto][alpha][beta] = new double*[13];
	    for(int i=0; i<13; i++) {
	      Wt[ibin][pto][alpha][beta][i] = new double[13];
	      for(int j=0; j<13; j++) {
		Wt[ibin][pto][alpha][beta][i][j] = 0;
		for(int rho=0; rho<nx1; rho++) {
		  for(int sigma=0; sigma<nx2; sigma++) {
		    double *f1 = M1a[alpha][rho][i];
		    double *f2 = M2a[beta][sigma][j];
		    // Evaluate luminosities
		    genpdf->evaluate(f1,f2,H);
		    // Combine luminosities with weights
		    for(int ip=0; ip<nsubproc; ip++) {
		      Wt[ibin][pto][alpha][beta][i][j] += W[rho][sigma][ip] * H[ip];
		    }
		    if(Wt[ibin][pto][alpha][beta][i][j] != 0) {
		      FlavMap[i][j] = true;
		      Channels[i][j] += Wt[ibin][pto][alpha][beta][i][j];
		    }
		  }
		}
	      }
	    }
	  }
	}
	cout << "  " << endl;
	delete[] W;
	delete[] H;
	delete[] M1a;
	delete[] M2a;
      }
    }
  }
  /*
  // Write Flavour Map
  cout << "Final flavour map:" << endl;
  cout << "  " << endl;
  cout << "   ";
  for(int i=0; i<13; i++) {
    cout << flavour[i] << " ";
  }
  cout << endl;
  for(int i=0; i<13; i++) {
    cout << flavour[i] << " ";
    for(int j=0; j<13; j++) {
      cout << FlavMap[i][j] << "  ";
    }
    cout << endl;
  }
  cout << "  " << endl;
  */
  // Find the independent channels.
  // First of all count how many non-zero elements are there.
  // This sets the maximum number of channels.
  int MaxChannels = 0;
  for(int i=0; i<13; i++) {
    for(int j=0; j<13; j++) {
      if(Channels[i][j] != 0) MaxChannels++;
    }
  }
  // Find channels
  vector<string> PartChann;
  string proc;
  for(int i=0; i<13; i++) {
    for(int j=0; j<13; j++) {
      if(FlavMap[i][j]) {
	FlavMap[i][j] = false;
	proc = flavour[i] + flavour[j];
	for(int k=0; k<13; k++) {
	  for(int l=0; l<13; l++) {
	    if(!(k == i && l == j)) {
	      double eps = 1e-6;
	      double ratio = abs( 1 - Channels[k][l] / Channels[i][j] );
	      if(ratio < eps) {
		proc += " + " + flavour[k] + flavour[l];
		FlavMap[k][l] = false;
	      }
	    }
	  }
	}
	PartChann.push_back(proc);
      }
    }
  }
  // Report how many channels have been found
  int Nchannels = PartChann.size();
  cout << "APPLgridEvolutionFilter() Info: Number of independent channels after evolution: " << Nchannels << endl;
  if(Nchannels > MaxChannels) {
    cout << "APPLgridEvolutionFilter() Error: The number of independent channels exceeds the maximum allowed" << endl;
    exit(-10);
  }
  //for(int ip=0; ip<Nchannels; ip++) cout << PartChann[ip] << endl;
  cout << "  " << endl;

  // Output grid
  // Fill vector of luminosities
  std::vector<int> luminosities;
  luminosities.push_back(Nchannels);

  for(int ip=0 ; ip<Nchannels; ip++) {

    luminosities.push_back(ip);
    int nproc = count(PartChann[ip].begin(),PartChann[ip].end(),'+') + 1;
    luminosities.push_back(nproc);

    // loop over the parton-parton combinations within each luminosity
    for(int iproc=0; iproc<nproc; iproc++) {
      int part1 = FlavourCode(PartChann[ip].substr(7*iproc,2));
      int part2 = FlavourCode(PartChann[ip].substr(7*iproc+2,2));

      luminosities.push_back(part1);
      luminosities.push_back(part2);
    }
  }

  string pdffile = "apfel.config";
  new lumi_pdf(pdffile,luminosities);

  // Open the output grid with some temporay parameters that will be
  // redefined below.
  appl::grid *og = NULL;
  og = new appl::grid(1,   Q20, Q20, 0,
		      30, 1e-5, 1e0, 5,    // Temporary parameters to be redefined below
		      g->Nobs(), g->obsmin(), g->obsmax(), 
		      pdffile,
		      g->leadingOrder(), g->nloops());

  for(int ibin=0; ibin<nbin; ibin++) {
    for(int o=0; o<=no; o++) {
  //for(int ibin=0; ibin<1; ibin++) {
    //for(int o=1; o<=1; o++) {
      int pto = o;

      // Redefine parameters
      int nx = max(n1b[ibin][pto],n2b[ibin][pto]) + 1;
      og->redefine(ibin, pto, 1, Q20, Q20, nx, xmin[ibin][pto], 1e0);

      for(int alpha=0; alpha<n1b[ibin][pto]+1; alpha++) {
	for(int beta=0; beta<n2b[ibin][pto]+1; beta++) {

	  double *weight = new double[Nchannels];
	  for(int ip=0; ip<Nchannels; ip++) {
	    int i = FlavourCode(PartChann[ip].substr(0,2)) + 6;
	    int j = FlavourCode(PartChann[ip].substr(2,2)) + 6;

	    weight[ip] = Wt[ibin][pto][alpha][beta][i][j];
	  }

	  // Fill output grid with weights
	  og->fill_index(alpha, beta, 0, ibin, weight, pto);
	}
      }
    }
  }

  // Finally dump to file 
  //og->Write("pollo.root");

  /*
    NQ2      = igrid->getNQ2();
    Nx1      = igrid->getNx1();
    Nx2      = igrid->getNx2();
    tauorder = igrid->tauorder();
    xintOrd  = igrid->yorder();
    Q2min    = igrid->getQ2min();
    Q2max    = igrid->getQ2max();
    x1min    = igrid->getx1min();
    x2min    = igrid->getx2min();
  */

  delete Wt;
  delete g;
  delete og;

  cout << "APPLgridEvolutionFilter() Info: Evolution included in " << argv[1] << endl;
  cout << "  " << endl;
  return 0;
}

// Function that returns the flavour code
int FlavourCode(string Flavour)
{
  int code;

  if(Flavour == "t~")      code = -6;
  else if(Flavour == "b~") code = -5;
  else if(Flavour == "c~") code = -4;
  else if(Flavour == "s~") code = -3;
  else if(Flavour == "u~") code = -2;
  else if(Flavour == "d~") code = -1;
  else if(Flavour == "g ") code =  0;
  else if(Flavour == "d ") code =  1;
  else if(Flavour == "u ") code =  2;
  else if(Flavour == "s ") code =  3;
  else if(Flavour == "c ") code =  4;
  else if(Flavour == "b ") code =  5;
  else if(Flavour == "t ") code =  6;
  else {
    cout << "APPLgridEvolutionFilter() Error: Unknown flavour" << endl;
    exit(-10);
  }
  return code;
}
