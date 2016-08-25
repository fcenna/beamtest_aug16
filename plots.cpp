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

int plots()
{
	// [V]
	Double_t Voltage[]  = {  75.     , 100.        , 130.     , 150.        , 170.                            };

	// [mV]
	Double_t MPVampl[]  = { 118.851  , 167.903     , 233.088  , 288.704     , 383.407                                };
	Double_t eMPVampl[] = {  23.1859 ,  24.8838    ,  34.1929 ,  41.2505    ,  52.2331                               };

	// [ns*mV] ????????
	Double_t MPVarea[]  = { 168.915  , 193.867     , 247.868  , 287.227     , 428.053                                };
	Double_t eMPVarea[] = {  31.8203 ,  28.1663    ,  31.6028 ,  27.8198    ,  33.6734                               };

	// [ns]
	Double_t sigma[]    = {   0.10305,   0.0723633 , 0.054047 ,   0.0474093 ,   0.0408785                       };

}
