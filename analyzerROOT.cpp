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

int analyzerROOT()
{
	TFile f_1001 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run1001.root"); // C1_100V_C2_100V_C3_100V_C4_28V_I
	TFile f_1302 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run1302.root"); // C1_130V_C2_130V_C3_130V_C4_28V_II
	TFile f_1502 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run1502.root"); // C1_150V_C2_150V_C3_150V_C4_28V_II
	TFile f_1901 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run1901.root"); // C1_190V_C2_190V_C3_190V_C4_28V_I
	TFile f_2001 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run2001.root"); // C1_200V_C2_200V_C3_200V_C4_28V_I

	TTree *tree1001 = (TTree*)f_1001.Get("Analysis");
	TTree *tree1302 = (TTree*)f_1302.Get("Analysis");
	TTree *tree1502 = (TTree*)f_1502.Get("Analysis");
	TTree *tree1901 = (TTree*)f_1901.Get("Analysis");
	TTree *tree2001 = (TTree*)f_2001.Get("Analysis");

	TCut cut1001_ch0 = ("ampl[0]>100 && ampl[0]<800 && t_frac30[0]>6 && t_frac30[0]<12");
	TCut cut1001_ch1 = ("ampl[1]>100 && ampl[1]<800 && t_frac30[1]>6 && t_frac30[1]<12");

	TCut cut2001_ch0 = ("ampl[0]>200 && ampl[0]<1050 && t_frac30[0]>6 && t_frac30[0]<12");
	TCut cut2001_ch1 = ("ampl[1]>200 && ampl[1]<1050 && t_frac30[1]>6 && t_frac30[1]<12");

	TCanvas *1001_deltatfrac30_distr = new TCanvas("1001_deltatfrac30_distr","UFSD Signal Time Gap (constant fraction 30% @100V)",1200,1200);
	
	tree1001
	tree1001->Draw("(t_frac30[0]-t_frac30[1])>>h1001_deltatfrac30",cut1001_ch0&&cut1001_ch1);
	
	
	tree1502->Draw("t_frac30[0]-t_frac30[1]",cut1001_ch0&&cut1001_ch1);
	tree1901->Draw("t_frac30[0]-t_frac30[1]",cut1901_ch0&&cut1901_ch1);
	tree2001->Draw("t_frac30[0]-t_frac30[1]",cut2001_ch0&&cut2001_ch1);

//	TCanvas *1001_ampl_distribution = new TCanvas("chprofEfficiencyVsZ","UFSD Signal Amplitude Distribution @100V",1200,1200);
//  	tree1001->Draw("ampl[0]",cut1001_ch0);
//  	tree1001->Draw("ampl[1]",cut1001_ch1, "SAME");
//  	tree1001->Draw("ampl[2]",cut1001_ch2, "SAME");

	tree1001.Delete();
	tree1302.Delete();
	tree1502.Delete();
	tree2001.Delete();

	f_1001.Close();
	f_1302.Close();
	f_1502.Close();
	f_2001.Close();

	return 0;
} 

