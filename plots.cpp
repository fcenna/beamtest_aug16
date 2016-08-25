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
	Double_t Voltage[]  = {  75.     , 100.        , 130.     , 150.        , 170.         , 180.        , 190.         , 200.        , 210.        , 220.                  };

	// [mV]
	Double_t MPVampl[]  = { 118.851  , 167.903     , 233.088  , 288.704     , 383.407      , 408.379     , 487.687      , 541.138     , 631.794     , 744.075               };
	Double_t eMPVampl[] = {  23.1859 ,  24.8838    ,  34.1929 ,  41.2505    ,  52.2331     ,  51.1958    ,  62.1559     ,  66.3681    ,  77.8171    ,  78.6304              };

	// [ns*mV] ????????
	Double_t MPVarea[]  = { 168.915  , 193.867     , 247.868  , 287.227     , 428.053      , 379.652     , 471.48       , 538.631     , 629.418     , 742.799               };
	Double_t eMPVarea[] = {  31.8203 ,  28.1663    ,  31.6028 ,  27.8198    ,  33.6734     ,  40.3216    ,  42.3391     ,  40.5427    ,  45.0933    ,  51.7833                                                                                        };

	// [ns]
	Double_t sigma[]    = {   0.10305,   0.0723633 , 0.054047 ,   0.0474093 ,   0.0408785  ,   0.0408986 ,   0.0377921  ,   0.0356133 ,   0.0331712 , 0.0373301             };
	
	Int_t nPoints = sizeof(Voltage)/sizeof(Double_t);
	Double_t esigma[];

	Int_t i=0;
	cout << nPoints << "	"<< sigma[0] << endl;
	for (i=0; i<nPoints ; i++) {
		cout << i << endl;
//		esigma[i] = 0.01*sigma[i];
	}
		


	TGraphErrors *sigmaVSvoltage = new TGraphErrors( nPoints, Voltage, sigma, NULL, esigma);
	TCanvas *Ccacca = new TCanvas("cacca", "cacca", 1200, 1200);

	sigmaVSvoltage->SetMarkerStyle(21);
	sigmaVSvoltage->SetMarkerSize(1.);
	sigmaVSvoltage->Draw("AP");

}
