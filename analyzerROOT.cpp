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
	TFile f_1501 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run1502.root"); // C1_150V_C2_150V_C3_150V_C4_28V_II

	TTree *tree1001 = (TTree*)f_1001.Get("Analysis");
	TTree *tree1501 = (TTree*)f_1501.Get("Analysis");

	TCut cut1001 = ("ampl[0]>100 && ampl[0]<800 && ampl[1]>100 && ampl[1]<800 && ampl[2]>100 && ampl[2]<800");

	TCanvas *1001_ampl_distribution = new TCanvas("chprofEfficiencyVsZ","UFSD Signal Amplitude Distribution @100V",1200,1200);
  	tree1001->Draw("ampl[0]",cut1001);
  	tree1001->Draw("ampl[1]",cut1001, "SAME");
  	tree1001->Draw("ampl[2]",cut1001, "SAME");

	f_1001.Close();
	f_1501.Close();

} 

