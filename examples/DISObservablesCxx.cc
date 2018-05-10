#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 
		1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options

  //
  // Settings
  //
  string proc = "NC";
  //APFEL::SetMassScheme("ZM-VFNS");
  APFEL::SetMassScheme("FONLL-C");
  APFEL::SetProcessDIS(proc);
  APFEL::SetQLimits(1.4,250);
  //APFEL::SetPolarizationDIS(0);
  //APFEL::SetProjectileDIS("electron");
  //APFEL::SetTargetDIS("proton");
  //APFEL::EnableTargetMassCorrections(false);
  //APFEL::EnableDampingFONLL(true);
  //APFEL::SetFastEvolution(true);
  //APFEL::LockGrids(true);
  //APFEL::EnableEvolutionOperator(true);
  //APFEL::SetFFNS(3);
  //APFEL::SetTheory("QUniD");
  //APFEL::SetTheory("QED");
  //APFEL::SetTheory("QUniD");
  //APFEL::EnableLeptonEvolution(true);
  //APFEL::SetTauMass(1e10);
  //APFEL::SetPerturbativeOrder(0);
  //APFEL::SetPDFEvolution("exactalpha");
  //APFEL::SetPDFSet("NNPDF30_nlo_as_0118");
  //APFEL::SetPDFSet("MRST2004qed");
  //APFEL::SetNumberOfGrids(1);
  //APFEL::SetGridParameters(1,30,3,1e-5);
  //APFEL::SetGridParameters(2,30,3,2e-1);
  //APFEL::SetGridParameters(3,30,3,8e-1);
  //APFEL::SetPDFSet("NNPDF30_nnlo_as_0118");
  //APFEL::SetAlphaQCDRef(0.118,91.2);
  //APFEL::SetAlphaEvolution("expanded");
  //APFEL::SetPDFEvolution("expandalpha");
  //APFEL::SetPoleMasses(1.275,4.18,173.03);
  //APFEL::SetMaxFlavourPDFs(5);
  //APFEL::SetMaxFlavourAlpha(5);

  // Initializes integrals on the grids
  APFEL::InitializeAPFEL_DIS();

  double Q02, Q2, eps = 1e-10;
  cout << "Enter initial and final scale in GeV^2" << endl;
  cin >> Q02 >> Q2;

  // Load evolution
  double Q0 = sqrt(Q02) - eps;
  double Q  = sqrt(Q2);

  APFEL::ComputeStructureFunctionsAPFEL(Q0,Q);

  cout << scientific << setprecision(5) << endl;
  // Tabulate PDFs for the LHA x values
  cout << "alpha_QCD(mu2F) = " << APFEL::AlphaQCD(Q) << endl;
  cout << "alpha_QED(mu2F) = " << APFEL::AlphaQED(Q) << endl;
  cout << endl;

  cout << "   x   " 
       << setw(11) << "   F2light   " 
       << setw(11) << "   F2charm   " 
       << setw(11) << "   F2bottom  " 
       << setw(11) << "   F2total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::F2light(xlha[i])  << "  "
      	 << setw(11) << APFEL::F2charm(xlha[i])  << "  "
      	 << setw(11) << APFEL::F2bottom(xlha[i]) << "  "
      	 << setw(11) << APFEL::F2total(xlha[i])	 << endl;
  cout << "      " << endl;

  cout << "   x   " 
       << setw(11) << "   FLlight   " 
       << setw(11) << "   FLcharm   " 
       << setw(11) << "   FLbottom  " 
       << setw(11) << "   FLtotal   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::FLlight(xlha[i])  << "  "
      	 << setw(11) << APFEL::FLcharm(xlha[i])  << "  "
      	 << setw(11) << APFEL::FLbottom(xlha[i]) << "  "
      	 << setw(11) << APFEL::FLtotal(xlha[i])  << endl;
  cout << "      " << endl;

  cout << "   x   " 
       << setw(11) << "   F3light   " 
       << setw(11) << "   F3charm   " 
       << setw(11) << "   F3bottom  " 
       << setw(11) << "   F3total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::F3light(xlha[i])  << "  "
      	 << setw(11) << APFEL::F3charm(xlha[i])  << "  "
      	 << setw(11) << APFEL::F3bottom(xlha[i]) << "  "
      	 << setw(11) << APFEL::F3total(xlha[i])	 << endl;
  cout << "      " << endl;
  //
  // Cache Structure functions
  //
  APFEL::CacheStructureFunctionsAPFEL(Q0);

  cout << "   x   " 
       << setw(11) << "   F2light   " 
       << setw(11) << "   F2charm   " 
       << setw(11) << "   F2bottom  " 
       << setw(11) << "   F2total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","light",xlha[i],Q)  << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","charm",xlha[i],Q)  << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","bottom",xlha[i],Q) << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","total",xlha[i],Q)  << endl;
  cout << "      " << endl;

  cout << "   x   " 
       << setw(11) << "   FLlight   " 
       << setw(11) << "   FLcharm   " 
       << setw(11) << "   FLbottom  " 
       << setw(11) << "   FLtotal   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","light",xlha[i],Q)  << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","charm",xlha[i],Q)  << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","bottom",xlha[i],Q) << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","total",xlha[i],Q)  << endl;
  cout << "      " << endl;

  cout << "   x   " 
       << setw(11) << "   F3light   " 
       << setw(11) << "   F3charm   " 
       << setw(11) << "   F3bottom  " 
       << setw(11) << "   F3total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","light",xlha[i],Q)  << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","charm",xlha[i],Q)  << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","bottom",xlha[i],Q) << "  "
      	 << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","total",xlha[i],Q)  << endl;
  cout << "      " << endl;

  return 0;
}
