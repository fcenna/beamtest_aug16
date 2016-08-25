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

	TCut C0_75V = ( "ampl[0]>60 && ampl[0]<750 && t_frac30[0]>=7 && t_frac30[0]<=10 ");
	TCut C1_75V = ( "ampl[1]>70 && ampl[1]<750 && t_frac30[1]>=7 && t_frac30[1]<=10 ");

//	TCanvas *c_CACCA = new TCanvas(" CACCA ",640,480);
	
	t_complete->Draw("(ampl[0])>>h_ampl_C0_75V", Run75V && C0_75V );
	h_ampl_C0_75V->Fit("landau"); //"fit_ampl_C0_75V", 
	Double_t MPVampl_C0_75V = landau->GetParameter(1); //Get most probable value of amplitude
	Double_t eMPVampl_C0_75V = landau->GetParameter(2); //Get error of MPV as sigma of landau

	t_complete->Draw("(area[0])>> h_area_C0_75V ", Run75V && C0_75V );
	h_area_C0_75V->Fit("landau"); //"fit_area_C0_75V", 
	Double_t  MPVarea_C0_75V = landau->GetParameter(1);  //Get most probable value of area
	Double_t eMPVarea_C0_75V = landau->GetParameter(2);  //Get error of MPV as sigma of landau
	

	t_complete->Draw("(t_frac30[1]-t_frac30[0])>>h_deltat30_75V", Run75V && C0_75V && C1_75V);
	h_deltat30_75V->Fit("gaus"); //"fit_ampl_C0_75V", 
	Double_t sigmat_75V = (1/sqrt(2))*gaus->GetParameter(2); //Get time resolution as the stddev of Gaussian distribution

	cout << "MPVampl_C0_75V  [mV]= " << MPVampl_C0_75V  << endl;
	cout << "eMPVampl_C0_75V [mV]= " << eMPVampl_C0_75V << endl;
	cout << "MPVarea_C0_75V  [mV]= " << MPVarea_C0_75V  << endl;
	cout << "eMPVarea_C0_75V [mV]= " << eMPVarea_C0_75V << endl;
	cout << "sigmat_75V      [ns]= " << sigmat_75V      << endl;

	landau->Delete();
	gaus->Delete();

	t_complete->Delete();

	f_complete.Close();

	return 0;
} 

