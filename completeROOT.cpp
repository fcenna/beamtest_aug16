// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVirtualFFT.h"

// Standard includes
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

int completeROOT()
{
	TFile f_complete ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/75V_100V_130V_150V_170V_180V_190V_200V_210V_220Va_240Va.root"); 

	TTree *t_complete = (TTree*)f_complete.Get("Analysis");

	//TCut cut1001_ch0 = ("ampl[0]>100 && ampl[0]<800 && t_frac30[0]>6 && t_frac30[0]<12");
	TCut Run75V  = ( "nrun == 75");
	TCut Run100V = ( "nrun == 100");
	TCut Run130V = ( "nrun == 130");
	TCut Run150V = ( "nrun == 150");
	TCut Run170V = ( "nrun == 170");
	TCut Run180V = ( "nrun == 180");
	TCut Run190V = ( "nrun == 190");
	TCut Run200V = ( "nrun == 200");
	TCut Run210V = ( "nrun == 210");
	TCut Run220V = ( "nrun == 220");
	TCut Run240V = ( "nrun == 240");

	TCut C0_75V  = ( "ampl[0]>60  && ampl[0]<750 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_75V  = ( "ampl[1]>70  && ampl[1]<750 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

	TCut C0_100V = ( "ampl[0]>100 && ampl[0]<750 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_100V = ( "ampl[1]>100 && ampl[1]<750 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

	TCut C0_130V = ( "ampl[0]>100 && ampl[0]<750 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_130V = ( "ampl[1]>100 && ampl[1]<750 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

	TCut C0_170V = ( "ampl[0]>200 && ampl[0]<950 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_170V = ( "ampl[1]>200 && ampl[1]<950 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

	TCut C0_180V = ( "ampl[0]>300 && ampl[0]<1100 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_180V = ( "ampl[1]>300 && ampl[1]<1100 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

	TCut C0_210V = ( "ampl[0]>400 && ampl[0]<1100 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_210V = ( "ampl[1]>400 && ampl[1]<1100 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

	TCut C0_220V = ( "ampl[0]>260 && ampl[0]<800 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_220V = ( "ampl[1]>260 && ampl[1]<800 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

//	TCanvas *c_CACCA = new TCanvas(" CACCA ",640,480);

	t_complete->Draw("(2*ampl[0])>>h_ampl_C0", Run220V && C0_220V );

//	t_complete->Draw("(ampl[0])>>h_ampl_C0", Run75V && C0_75V );
	h_ampl_C0->Fit("landau"); //"fit_ampl_C0", 
	Double_t MPVampl_C0 = landau->GetParameter(1); //Get most probable value of amplitude
	Double_t eMPVampl_C0 = landau->GetParameter(2); //Get error of MPV as sigma of landau

	t_complete->Draw("(2*area[0])>> h_area_C0 ", Run220V && C0_220V );
	h_area_C0->Fit("landau"); //"fit_area_C0_75V", 
	Double_t  MPVarea_C0 = landau->GetParameter(1);  //Get most probable value of area
	Double_t eMPVarea_C0 = landau->GetParameter(2);  //Get error of MPV as sigma of landau
	

	t_complete->Draw("(t_frac30[1]-t_frac30[0])>>h_deltat30", Run220V && C0_220V && C1_220V);
	h_deltat30->Fit("gaus"); //"fit_ampl_C0", 
	Double_t sigmat = (1/sqrt(2))*gaus->GetParameter(2); //Get time resolution as the stddev of Gaussian distribution

	cout << endl << "********* RESULTS ********"<< endl;
	cout << " MPVampl  [mV]= " << MPVampl_C0  << endl;
	cout << " eMPVampl [mV]= " << eMPVampl_C0 << endl;
	cout << " MPVarea  [mV]= " << MPVarea_C0  << endl;
	cout << " eMPVarea [mV]= " << eMPVarea_C0 << endl;
	cout << " sigma    [ns]= " << sigmat      << endl;
	cout << "**************************"<< endl;

	landau->Delete();
	gaus->Delete();
	//c1->close();

	t_complete->Delete();

	f_complete.Close();

	return 0;
} 

