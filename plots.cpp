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

	// [mV] MPV and Sigma from Landau fit
	Double_t MPVampl[]  = { 118.851  , 167.903     , 233.088  , 288.704     , 383.407      , 408.379     , 487.687      , 541.138     , 631.794     , 744.075               };
	Double_t eMPVampl[] = {  23.1859 ,  24.8838    ,  34.1929 ,  41.2505    ,  52.2331     ,  51.1958    ,  62.1559     ,  66.3681    ,  77.8171    ,  78.6304              };

	// [ns*mV] MPV and sigma from Landau fit
	Double_t MPVarea[]  = { 168.915  , 193.867     , 247.868  , 287.227     , 428.053      , 379.652     , 471.48       , 538.631     , 629.418     , 742.799               };
	Double_t eMPVarea[] = {  31.8203 ,  28.1663    ,  31.6028 ,  27.8198    ,  33.6734     ,  40.3216    ,  42.3391     ,  40.5427    ,  45.0933    ,  51.7833                                                                                        };

	// [ns] Time resolution comparing 2 LGAD  -- TFrac30
	Double_t sigma[]    = {   103.05,   72.3633 , 54.047 ,   47.4093 ,   40.8785  ,   40.8986 ,   37.7921  ,   35.6133 ,   33.1712 , 37.3301             };
	
	// from UCSC slides GAIN @200V = 14.1284759029
	// gain(@V)= ( gain(@200V)/MPVampl(@200V) ) *MPVampl(@V)
	Double_t GainUCSC_200V = 14.1284759029 ;
	Double_t scaling_factor = GainUCSC_200V/MPVampl[7]; // !!! change MPVampl index if more runs are added
	
	Int_t nPoints = sizeof(Voltage)/sizeof(Double_t);

	Double_t *gain = new Double_t [nPoints]; // dynamic allocation in order to use non static constant index 
	Double_t *egain = new Double_t [nPoints]; // dynamic allocation in order to use non static constant index 
//	Double_t *esigma = new Double_t [nPoints]; // dynamic allocation in order to use non static constant index 

	for (Int_t i=0; i<nPoints ; i++) {
		gain[i] = scaling_factor * MPVampl[i];
		egain[i] = scaling_factor * eMPVampl[i];
//		esigma[i] = 0.01*sigma[i];
	}	
	

/*		
TCanvas *c1 = new TCanvas ("c1","Time Resolution",1200,1200);		
c1->Divide(2,2);
c1->cd(1);
TGraphErrors *sigmaVsVoltage = new TGraphErrors( nPoints, Voltage, sigma, NULL, NULL);
sigmaVsVoltage->SetTitle("Time resolution versus Bias Voltage");
sigmaVsVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
sigmaVsVoltage->GetYaxis()->SetTitle("Time resolution [ps]");
sigmaVsVoltage->SetMarkerStyle(21);
sigmaVsVoltage->SetMarkerSize(1.);
sigmaVsVoltage->Draw("AP");

c1->cd(2);
TGraphErrors *sigmaVsMPV = new TGraphErrors( nPoints, MPVampl, sigma, eMPVampl, NULL);
sigmaVsMPV->SetTitle("Time resolution versus Signal Amplitude");
sigmaVsMPV->GetXaxis()->SetTitle("MPV [mV]");
sigmaVsMPV->GetYaxis()->SetTitle("Time resolution [ps]");
sigmaVsMPV->SetMarkerStyle(21);
sigmaVsMPV->SetMarkerSize(1.);
sigmaVsMPV->Draw("AP");

c1->cd(3);
TGraphErrors *sigmaVsArea = new TGraphErrors( nPoints, MPVarea, sigma, eMPVarea, NULL);
sigmaVsArea->SetTitle("Time resolution versus Signal Area");
sigmaVsArea->GetXaxis()->SetTitle("Area [mV*ns]");
sigmaVsArea->GetYaxis()->SetTitle("Time resolution [ps]");
sigmaVsArea->SetMarkerStyle(21);
sigmaVsArea->SetMarkerSize(1.);
sigmaVsArea->Draw("AP");

c1->cd(4);
TGraphErrors *MPVsArea = new TGraphErrors( nPoints, MPVarea, MPVampl, eMPVarea, eMPVampl);
MPVsArea->SetTitle("Signal Amplitude versus Signal Area");
MPVsArea->GetXaxis()->SetTitle("MPV [mV]");
MPVsArea->GetYaxis()->SetTitle("Area [mV*ns]");
MPVsArea->SetMarkerStyle(21);
MPVsArea->SetMarkerSize(1.);
MPVsArea->Draw("AP");

TCanvas *c1a = new TCanvas ("c1a","Time Resolution",1200,1200);		
c1a->cd(1);
TGraphErrors *sigmaVsVoltage = new TGraphErrors( nPoints, Voltage, sigma, NULL, NULL);
sigmaVsVoltage->SetTitle("Time resolution vs Bias Voltage");
sigmaVsVoltage->GetXaxis()->SetTitleOffset(1.03);
sigmaVsVoltage->GetYaxis()->SetTitleOffset(1.03);
sigmaVsVoltage->GetXaxis()->SetTitleSize(0.045);
sigmaVsVoltage->GetYaxis()->SetTitleSize(0.045);
sigmaVsVoltage->GetXaxis()->SetLabelSize(0.045);
sigmaVsVoltage->GetYaxis()->SetLabelSize(0.045);
sigmaVsVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
sigmaVsVoltage->GetYaxis()->SetTitle("Time resolution [ps]");
sigmaVsVoltage->SetMarkerStyle(21);
sigmaVsVoltage->SetMarkerSize(1.);
sigmaVsVoltage->Draw("AP");
c1a->SaveAs("../../../public/plots/TimeRes_V.pdf");

TCanvas *c1b = new TCanvas ("c1b","Time Resolution",1200,1200);		
c1b->cd(1);
TGraphErrors *sigmaVsMPV = new TGraphErrors( nPoints, MPVampl, sigma, eMPVampl, NULL);
sigmaVsMPV->GetXaxis()->SetTitleOffset(1.03);
sigmaVsMPV->GetYaxis()->SetTitleOffset(1.03);
sigmaVsMPV->GetXaxis()->SetTitleSize(0.045);
sigmaVsMPV->GetYaxis()->SetTitleSize(0.045);
sigmaVsMPV->GetXaxis()->SetLabelSize(0.045);
sigmaVsMPV->GetYaxis()->SetLabelSize(0.045);
sigmaVsMPV->SetTitle("Time resolution vs Signal Amplitude");
sigmaVsMPV->GetXaxis()->SetTitle("Amplitude [mV]");
sigmaVsMPV->GetYaxis()->SetTitle("Time resolution [ps]");
sigmaVsMPV->SetMarkerStyle(21);
sigmaVsMPV->SetMarkerSize(1.);
sigmaVsMPV->Draw("AP");
c1b->SaveAs("../../../public/plots/TimeRes_Ampl.pdf");

TCanvas *c1c = new TCanvas ("c1c","Time Resolution",1200,1200);		
c1c->cd(1);
TGraphErrors *sigmaVsArea = new TGraphErrors( nPoints, MPVarea, sigma, eMPVarea, NULL);
sigmaVsArea->SetTitle("Time resolution vs Signal Area");
sigmaVsArea->GetXaxis()->SetTitleOffset(1.03);
sigmaVsArea->GetYaxis()->SetTitleOffset(1.03);
sigmaVsArea->GetXaxis()->SetTitleSize(0.045);
sigmaVsArea->GetYaxis()->SetTitleSize(0.045);
sigmaVsArea->GetXaxis()->SetLabelSize(0.045);
sigmaVsArea->GetYaxis()->SetLabelSize(0.045);
sigmaVsArea->GetXaxis()->SetTitle("Area [mV*ns]");
sigmaVsArea->GetYaxis()->SetTitle("Time resolution [ps]");
sigmaVsArea->SetMarkerStyle(21);
sigmaVsArea->SetMarkerSize(1.);
sigmaVsArea->Draw("AP");
c1c->SaveAs("../../../public/plots/TimeRes_Area.pdf");

TCanvas *c1d = new TCanvas ("c1d","Time Resolution",1200,1200);		
c1d->cd(1);
TGraphErrors *MPVsArea = new TGraphErrors( nPoints, MPVampl, MPVarea, eMPVampl, eMPVarea);
MPVsArea->SetTitle("Signal Area vs Signal Amplitude");
MPVsArea->GetXaxis()->SetTitleOffset(1.03);
MPVsArea->GetYaxis()->SetTitleOffset(1.03);
MPVsArea->GetXaxis()->SetTitleSize(0.045);
MPVsArea->GetYaxis()->SetTitleSize(0.045);
MPVsArea->GetXaxis()->SetLabelSize(0.045);
MPVsArea->GetYaxis()->SetLabelSize(0.045);
MPVsArea->GetXaxis()->SetTitle("Amplitude [mV]");
MPVsArea->GetYaxis()->SetTitle("Area [mV*ns]");
MPVsArea->SetMarkerStyle(21);
MPVsArea->SetMarkerSize(1.);
MPVsArea->Draw("AP");
c1d->SaveAs("../../../public/plots/Ampl_Area.pdf");
*/


TCanvas *c1e = new TCanvas ("c1e","Amplitude",1920,1080);		
c1e->cd(1);
TGraphErrors *amplMPVvsVoltage = new TGraphErrors( nPoints, Voltage, MPVampl, NULL, eMPVampl);
amplMPVvsVoltage->SetTitle("Signal Amplitude vs Bias Voltage");
amplMPVvsVoltage->GetXaxis()->SetTitleOffset(1.03);
amplMPVvsVoltage->GetYaxis()->SetTitleOffset(1.03);
amplMPVvsVoltage->GetXaxis()->SetTitleSize(0.040);
amplMPVvsVoltage->GetYaxis()->SetTitleSize(0.040);
amplMPVvsVoltage->GetXaxis()->SetLabelSize(0.040);
amplMPVvsVoltage->GetYaxis()->SetLabelSize(0.040);
amplMPVvsVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
amplMPVvsVoltage->GetYaxis()->SetTitle("Signal Amplitude [mV]");
amplMPVvsVoltage->SetMarkerStyle(21);
amplMPVvsVoltage->SetMarkerSize(1.);
amplMPVvsVoltage->Draw("AP");
c1e->SaveAs("../../../public/plots/Amplitude_V.pdf");

TCanvas *c1f = new TCanvas ("c1f","Gain",1920,1080);		
c1f->cd(1);
TGraphErrors *GainvsVoltage = new TGraphErrors( nPoints, Voltage, gain, NULL, egain);
GainvsVoltage->SetTitle("Gain vs Bias Voltage");
GainvsVoltage->GetXaxis()->SetTitleOffset(1.03);
GainvsVoltage->GetYaxis()->SetTitleOffset(1.03);
GainvsVoltage->GetXaxis()->SetTitleSize(0.040);
GainvsVoltage->GetYaxis()->SetTitleSize(0.040);
GainvsVoltage->GetXaxis()->SetLabelSize(0.040);
GainvsVoltage->GetYaxis()->SetLabelSize(0.040);
GainvsVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
GainvsVoltage->GetYaxis()->SetTitle("Gain");
GainvsVoltage->SetMarkerStyle(21);
GainvsVoltage->SetMarkerSize(1.);

GainvsVoltage->Fit("expo");

GainvsVoltage->Draw("AP");
expo->Draw("SAME");
c1f->SaveAs("../../../public/plots/Gain_V.pdf");

	return 0;
}
