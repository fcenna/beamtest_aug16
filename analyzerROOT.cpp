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

// Lorenzo includes
#include "Analysis_40Gs.h"

using namespace std;

int analyzerROOT()
{
	TFile f_1501 ("/afs/cern.ch/work/r/rmulargi/private/BeamTestAug2016/DATA/Run1502.root"); // C1_150V_C2_150V_C3_150V_C4_28V_II
	TTree *tree1501 = (TTree*)f_1501.Get("Analysis");
	f_1501.Close();



	
} 

