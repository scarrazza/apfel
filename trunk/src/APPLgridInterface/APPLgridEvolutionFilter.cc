// APPLgridEvolutionFilter.cc
// Add to a pre-existing APPLgrid root file the PDF evolution provided by APFEL
//
// valerio.bertone@cern.ch

// Standard libraries
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

// APFEL
#include "APFEL/APFEL.h"
#include "APPLgridInterface.h"

// APPLgrid
#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"
#include "appl_grid/lumi_pdf.h"
#include "appl_grid/SparseMatrix3d.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(int argc, char* argv[]) {

  // Check that the input is correct
  if (argc!=2) {
    cout << "Invalid Parameters:" << endl;
    cout << "Syntax: ./APPLgridEvolutionFilter <APPLgrid file>" << endl;
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

  // Declare the arrays of weights and other stuff
  double *******Wt = new double******[nbin];
  double ******W   = new double*****[nbin];
  double ***x1v   = new double**[nbin];
  double ***x2v   = new double**[nbin];
  double **xmin    = new double*[nbin];
  int **n1b        = new int*[nbin];
  int **n2b        = new int*[nbin];
  int **ntau       = new int*[nbin];

  // Loop over the bins
  //for(int ibin=0; ibin<nbin; ibin++) {
  for(int ibin=0; ibin<1; ibin++) {
    cout << "APPLgridEvolutionFilter() Info: bin " << ibin+1 << " of " << nbin << endl;
    cout << "  " << endl;

    Wt[ibin]   = new double*****[no];
    W[ibin]    = new double****[no];
    x1v[ibin]  = new double*[no];
    x2v[ibin]  = new double*[no];
    xmin[ibin] = new double[no];
    n1b[ibin]  = new int[no];
    n2b[ibin]  = new int[no];
    ntau[ibin] = new int[no];

    // Loop over the perturbative orders
    //for(int pto=0; pto<=no; pto++) {
    for(int pto=1; pto<2; pto++) {
      if(pto == 0) cout << "APPLgridEvolutionFilter() Info: Perturbative order: LO " << endl;
      if(pto == 1) cout << "APPLgridEvolutionFilter() Info: Perturbative order: NLO " << endl;
      if(pto == 2) cout << "APPLgridEvolutionFilter() Info: Perturbative order: NNLO " << endl;
      cout << "  " << endl;

      // Access internal grid for "pto" and "ibin"
      igrid = g->weightgrid(pto,ibin);
      if (igrid==NULL) cerr << "APPLgridEvolutionFilter() Error: grid not found"<<endl;

      // Parameter of the current internal grid
      ntau[ibin][pto] = igrid->Ntau();
      int nx1         = igrid->Ny1();
      int nx2         = igrid->Ny2();
      //cout << "   Number of Q2 grid points = " << ntau[ibin][pto] << endl;
      //cout << "   Number of x1 grid points = " << nx1 << endl;
      //cout << "   Number of x2 grid points = " << nx2 << endl;
      //cout << "  " << endl;

      Wt[ibin][pto]  = new double****[ntau[ibin][pto]];
      W[ibin][pto]   = new double***[ntau[ibin][pto]];

      // Loop over the grid in Q2
      for(int itau=0; itau<ntau[ibin][pto]; itau++) {
      //for(int itau=1; itau<2; itau++) {
	cout << "APPLgridEvolutionFilter() Info: Q2 grid point " << itau+1 << " of " << ntau[ibin][pto] << endl;
	cout << "  " << endl;

	double Q = sqrt(igrid->fQ2(igrid->gettau(itau)));
	//Q0  = sqrt(igrid->fQ2(igrid->gettau(1)));;
	//Q20 = Q0 * Q0;
	cout << "APPLgridEvolutionFilter() Info: Evolution between " << Q0 << " GeV and " << Q << " GeV" << endl;

	double delta = 0,deltap = 0;
	double a,ap;
	double eps = 1e-8;

	// Generate the evolution table for the first grid in x
	vector<double> xext1v;
	int iy1 = nx1;

	for(int ix1=0; ix1<nx1; ix1++) {
	  iy1--; // inverse order
	  xext1v.push_back(igrid->fx(igrid->gety1(iy1)));

	  // Parameters of the grid (see eq. (1) of arXiv:0911.2985)
	  if(ix1 > 0) {
	    deltap = delta;
	    ap     = a;
	    delta  = igrid->gety1(iy1) - igrid->gety1(iy1+1);
	    a      = ( delta - log( xext1v[ix1-1] / xext1v[ix1] ) ) / ( xext1v[ix1-1] - xext1v[ix1] );

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

	// Generate the evolution table for the second grid in x.
	// if the x-space grids are not equal for the PDFs ...
	// INFO: Usually APPLgrid says that they are not equal while they are.
	vector<double> xext2v;
	if(!sym) {
	  int iy2 = nx2;

	  for(int ix2=0; ix2<nx2; ix2++) {
	    iy2--; // inverse order
	    xext2v.push_back(igrid->fx(igrid->gety2(iy2)));

	    // Parameters of the grid (see eq. (1) of arXiv:0911.2985)
	    if(ix2 > 0) {
	      deltap = delta;
	      ap     = a;
	      delta  = igrid->gety2(iy2) - igrid->gety2(iy2+1);
	      a      = ( delta - log( xext2v[ix2-1] / xext2v[ix2] ) ) / ( xext2v[ix2-1] - xext2v[ix2] );

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
	}
	// ... otherwise copy "n1b", "xext1v" and "M1" into "n2b", "xext2v" and "M2"
	else {
	  n2b[ibin][pto] = n1b[ibin][pto];
	  xext2v = xext1v;
	}

	// Extract weights from the APPLgrid root file (in the inverse order)
	cout << "APPLgridEvolutionFilter() Info: Extracting original weights ..." << endl;
	bool NonZeroWeights = false;
	W[ibin][pto][itau] = new double**[nx1];
	for(int alpha=0; alpha<nx1; alpha++) {
	  W[ibin][pto][itau][alpha] = new double*[nx2];
	  for(int beta=0; beta<nx2; beta++) {
	    W[ibin][pto][itau][alpha][beta] = new double[nsubproc];
	    for(int ip=0; ip<nsubproc; ip++) {
	      W[ibin][pto][itau][alpha][beta][ip] = (*(const SparseMatrix3d*) igrid->weightgrid(ip))(itau,nx1-1-alpha,nx2-1-beta);
	      if(W[ibin][pto][itau][alpha][beta][ip] != 0) NonZeroWeights = true;
	    }
	  }
	}

	if(NonZeroWeights) {
	  // Call APFEL and compute the evolution operators
	  APFEL::SetPerturbativeOrder(0);
	  //APFEL::SetMaxFlavourPDFs(4);
	  //APFEL::SetMaxFlavourAlpha(4);
	  APFEL::SetFFNS(3);
	  double *M1 = NULL;
	  M1 = new double[14*14*(n1b[ibin][pto]+1)*(n1b[ibin][pto]+1)];
	  double *M2 = NULL;
	  M2 = new double[14*14*(n2b[ibin][pto]+1)*(n2b[ibin][pto]+1)];
	  cout << "APPLgridEvolutionFilter() Info: Computing evolution operator on grid 1 ..." << endl;
	  APFEL::ExternalEvolutionOperator(Q0,Q,n1b[ibin][pto],&xext1v[0],M1);
	  cout << "APPLgridEvolutionFilter() Info: Computing evolution operator on grid 2 ..." << endl;
	  if(!sym) {
	    APFEL::ExternalEvolutionOperator(Q0,Q,n2b[ibin][pto],&xext2v[0],M2);
	    xmin[ibin][pto] = min(xext1v[0],xext2v[0]);
	  }
	  else {
	    copy(M1,M1+14*14*(n1b[ibin][pto]+1)*(n1b[ibin][pto]+1),M2);
	    xmin[ibin][pto] = xext1v[0];
	  }

	  x1v[ibin][pto] = new double[n1b[ibin][pto]+1];
	  x2v[ibin][pto] = new double[n2b[ibin][pto]+1];
	  copy(xext1v.begin(),xext1v.end(),x1v[ibin][pto]);
	  copy(xext2v.begin(),xext2v.end(),x2v[ibin][pto]);

	  // Put the first evolution operator in a suitable matrix form
	  double ****M1a = new double***[n1b[ibin][pto]+1];
	  for(int alpha=0; alpha<n1b[ibin][pto]+1; alpha++) {
	    M1a[alpha] = new double**[n1b[ibin][pto]+1];
	    for(int rho=0; rho<n1b[ibin][pto]+1; rho++) {
	      M1a[alpha][rho] = new double*[13];
	      for(int k=0; k<13; k++) {
		M1a[alpha][rho][k] = new double[13];
		for(int i=0; i<13; i++) {
		  int m1 = ( i + 1 ) + 14 * ( ( k + 1 ) + 14 * ( alpha + ( n1b[ibin][pto] + 1 ) * rho ) );
		  M1a[alpha][rho][k][i] = M1[m1];
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
	      for(int l=0; l<13; l++) {
		M2a[beta][sigma][l] = new double[13];
		for(int j=0; j<13; j++) {
		  int m2 = ( j + 1 ) + 14 * ( ( l + 1 ) + 14 * ( beta + ( n2b[ibin][pto] + 1 ) * sigma ) );
		  M2a[beta][sigma][l][j] = M2[m2];
		}
	      }
	    }
	  }
	  delete[] M2;
	  /*
	  // Test Evolution operators
	  APFEL::EnableWelcomeMessage(true);
	  APFEL::EnableEvolutionOperator(false);
	  APFEL::SetPDFSet("toyLH_FFN_LO.LHgrid");
	  APFEL::InitializeAPFEL();
	  APFEL::EvolveAPFEL(Q0,Q);

	  LHAPDF::initPDFSet("toyLH_FFN_LO.LHgrid");
	  LHAPDF::initPDF(0);

	  double fOP[n1b[ibin][pto]+1][13];
	  for(int alpha=0; alpha<n1b[ibin][pto]+1; alpha++) {
	    for(int i=3; i<10; i++) {
	      fOP[alpha][i] = 0;
	      for(int rho=0; rho<n1b[ibin][pto]+1; rho++) {
		for(int k=0; k<13; k++) {
                  fOP[alpha][i] += M1a[alpha][rho][k][i] * LHAPDF::xfx(x1v[ibin][pto][rho],Q0,k-6);
		}
	      }
	      //cout << i-6 << "  " << x1v[ibin][pto][alpha] << "  " << LHAPDF::xfx(x1v[ibin][pto][alpha],Q,i-6) << "   " << APFEL::xPDF(i-6,x1v[ibin][pto][alpha]) << "  " << fOP[alpha][i] << endl;
	      double reldiff = 100 * ( LHAPDF::xfx(x1v[ibin][pto][alpha],Q,i-6) - fOP[alpha][i] ) / LHAPDF::xfx(x1v[ibin][pto][alpha],Q,i-6);
	      if(abs(reldiff) > 1) cout << i-6 << "  " << x1v[ibin][pto][alpha] << "  " << reldiff << endl;
	    }
	  }
	  cout << "   " << endl;
	  for(int beta=0; beta<n2b[ibin][pto]+1; beta++) {
	    for(int j=3; j<10; j++) {
	      fOP[beta][j] = 0;
	      for(int sigma=0; sigma<n2b[ibin][pto]+1; sigma++) {
		for(int l=0; l<13; l++) {
                  fOP[beta][j] += M1a[beta][sigma][l][j] * LHAPDF::xfx(x1v[ibin][pto][sigma],Q0,l-6);
		}
	      }
	      //cout << i-6 << "  " << x1v[ibin][pto][beta] << "  " << LHAPDF::xfx(x1v[ibin][pto][beta],Q,j-6) << "   " << APFEL::xPDF(j-6,x1v[ibin][pto][beta]) << "  " << fOP[beta][j] << endl;
	      double reldiff = 100 * ( LHAPDF::xfx(x1v[ibin][pto][beta],Q,j-6) - fOP[beta][j] ) / LHAPDF::xfx(x1v[ibin][pto][beta],Q,j-6);
	      if(abs(reldiff) > 1) cout << j-6 << "  " << x1v[ibin][pto][beta] << "  " << reldiff << endl;
	    }
	  }
	  */

	  // Combine APPLgrid weights "W" with the evolution operators "M1a" and "M2a" and produce "Wt"
	  double *H = new double[nsubproc];
	  Wt[ibin][pto][itau] = new double***[n1b[ibin][pto]+1];
	  for(int rho=0; rho<n1b[ibin][pto]+1; rho++) {
	    double perc = (double) 100 * rho / n1b[ibin][pto];
	    cout << "APPLgridEvolutionFilter() Info: Combining APPLgrid weights with evolution operators: " << setprecision(3) << setw(4) << perc << " % done" << endl;
	    if(rho != n1b[ibin][pto]) cout << "\x1b[A";
	    Wt[ibin][pto][itau][rho] = new double**[n2b[ibin][pto]+1];
	    for(int sigma=0; sigma<n2b[ibin][pto]+1; sigma++) {
	      Wt[ibin][pto][itau][rho][sigma] = new double*[13];
	      for(int k=0; k<13; k++) {
		Wt[ibin][pto][itau][rho][sigma][k] = new double[13];
		for(int l=0; l<13; l++) {
		  Wt[ibin][pto][itau][rho][sigma][k][l] = 0;
		  for(int alpha=0; alpha<nx1; alpha++) {
		    for(int beta=0; beta<nx2; beta++) {
		      double *f1 = M1a[alpha][rho][k];
		      double *f2 = M2a[beta][sigma][l];
		      // Evaluate luminosities
		      genpdf->evaluate(f1,f2,H);
		      // Combine luminosities with weights
		      for(int ip=0; ip<nsubproc; ip++) {
			Wt[ibin][pto][itau][rho][sigma][k][l] += W[ibin][pto][itau][alpha][beta][ip] * H[ip];
		      }
		    }
		  }
		  if(Wt[ibin][pto][itau][rho][sigma][k][l] != 0) {
		    FlavMap[k][l] = true;
		    Channels[k][l] += Wt[ibin][pto][itau][rho][sigma][k][l];
		  }
		}
	      }
	    }
	  }
	  delete[] H;
	  delete[] M1a;
	  delete[] M2a;
	}
	else {
	  cout << "APPLgridEvolutionFilter() Info: All the weights are null: Evolution not needed" << endl;
	  Wt[ibin][pto][itau] = new double***[n1b[ibin][pto]+1];
	  for(int rho=0; rho<n1b[ibin][pto]+1; rho++) {
	    Wt[ibin][pto][itau][rho] = new double**[n2b[ibin][pto]+1];
	    for(int sigma=0; sigma<n2b[ibin][pto]+1; sigma++) {
	      Wt[ibin][pto][itau][rho][sigma] = new double*[13];
	      for(int k=0; k<13; k++) {
		Wt[ibin][pto][itau][rho][sigma][k] = new double[13];
		for(int l=0; l<13; l++) {
		  Wt[ibin][pto][itau][rho][sigma][k][l] = 0;
		}
	      }
	    }
	  }
	}
	cout << "  " << endl;
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
  // Check whether we are dealing with W+, W- or none of them
  int ckmcharge = genpdf->getckmcharge();
  if(ckmcharge == 1 || ckmcharge == -1) luminosities.push_back(ckmcharge);
  else                                  luminosities.push_back(0);

  string pdffile = "apfel.config";
  new lumi_pdf(pdffile,luminosities);

  // Open the output grid with some temporay parameters that will be
  // redefined below.
  const appl::igrid* oigrid;
  appl::grid *og = NULL;
  og = new appl::grid( 1,  Q20, Q20, 0,
  		      30, 1e-5,   1, 5,    // Temporary parameters to be adjusted below
  		      g->Nobs()+1, g->obsmin(), g->obsmax(), 
  		      pdffile,
  		      g->leadingOrder(), g->nloops());

  //for(int ibin=0; ibin<nbin; ibin++) {
  for(int ibin=0; ibin<1; ibin++) {
    //for(int pto=0; pto<=no; pto++) {
    for(int pto=1; pto<2; pto++) {

      // Redefine parameters
      og->redefine(ibin, pto, 1, Q20, Q20, 30, xmin[ibin][pto], 1);

      // Access internal grid
      oigrid = og->weightgrid(pto,ibin);

      for(int rho=0; rho<n1b[ibin][pto]+1; rho++) {
         for(int sigma=0; sigma<n2b[ibin][pto]+1; sigma++) {
	  double *weight = new double[Nchannels];
	  for(int ip=0; ip<Nchannels; ip++) {
	    int k = FlavourCode(PartChann[ip].substr(0,2)) + 6;
	    int l = FlavourCode(PartChann[ip].substr(2,2)) + 6;
	    weight[ip] = 0;
	    for(int itau=0; itau<ntau[ibin][pto]; itau++) {
	    //for(int itau=1; itau<2; itau++) {
	      weight[ip] += Wt[ibin][pto][itau][rho][sigma][k][l];
	    }
	  }
	  // Fill output grid with weights
	  og->fill(x1v[ibin][pto][rho], x2v[ibin][pto][sigma], Q20, og->obs(ibin), weight, pto);
	}
      }
    }
  }

  // Finally dump to file 
  og->Write(string(argv[1]).substr(0,string(argv[1]).size()-5) + "-Evolved.root");

  delete W;
  delete Wt;
  delete g;
  delete genpdf;
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
