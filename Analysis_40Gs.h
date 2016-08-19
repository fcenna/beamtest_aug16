#ifndef WAVE_ANALYZER_H
#define WAVE_ANALYZER_H

#define NCHRO 4
#define MAX_CAMP 20000
#define CAMPFACT 1
#define MAX_CAMPREC MAX_CAMP/CAMPFACT+100

Int_t osc[NCHRO];
Int_t chan[NCHRO];
Double_t bck[NCHRO];
Double_t t_bck[NCHRO];
Double_t t_pul[NCHRO];
Int_t camp[NCHRO];
Int_t camprec[NCHRO];
Double_t FreqCut[NCHRO];

Double_t DT;
Double_t timeS[MAX_CAMP]; 

Double_t Itot[NCHRO][MAX_CAMP];

Double_t amp[NCHRO][MAX_CAMP];
Double_t m_amp[NCHRO][MAX_CAMP];
Double_t der_amp[NCHRO][MAX_CAMP];
Double_t amprec[NCHRO][MAX_CAMPREC];
Double_t m_amprec[NCHRO][MAX_CAMPREC];
Double_t d_amprec[NCHRO][MAX_CAMP];
Double_t FFT_real[NCHRO][MAX_CAMP];
Double_t FFT_comp[NCHRO][MAX_CAMP];
Double_t FFT_abs[NCHRO][MAX_CAMP];

Double_t FFT_real_cut[NCHRO][MAX_CAMP];
Double_t FFT_comp_cut[NCHRO][MAX_CAMP];
Double_t FFT_abs_cut[NCHRO][MAX_CAMP];

Double_t freq[MAX_CAMPREC];
Double_t timerec[MAX_CAMPREC];


  //Analysis variables

Double_t day;
Double_t month;
Double_t year;
Double_t hour;
Double_t minute;

Double_t MaxDeriv;
Double_t sec[NCHRO];
Double_t ampl[NCHRO];
Double_t area[NCHRO];
Double_t dVdt1030[NCHRO];
Double_t dVdt3070[NCHRO];
Double_t dVdt2080[NCHRO];
Double_t dVdt1090[NCHRO];

  Double_t totcha[NCHRO];
  Double_t maxS[NCHRO];
  Double_t maxD[NCHRO];
  Double_t max_bck[NCHRO];
  Double_t rms_bck[NCHRO];
  Double_t max_fwd_bck[NCHRO];
  Double_t rms_fwd_bck[NCHRO];
  Double_t charge[NCHRO];
  Int_t strip[NCHRO];
  Double_t FEcharge[NCHRO];
  Double_t FEchargeerr[NCHRO];
  Double_t totcharge[NCHRO];
  Double_t p_charge;
  Double_t p_chargeerr;
  Int_t sat[NCHRO];
  Double_t maxamp;
  Int_t maxstrip;
  Double_t eff_dis08[NCHRO];
  Double_t eff_dis10[NCHRO];
  Double_t eff_dis15[NCHRO];
  Double_t eff_dis20[NCHRO];
  Double_t eff_dis30[NCHRO];
  Double_t eff_dis40[NCHRO];
  Double_t eff_dis50[NCHRO];
  Double_t eff_dis60[NCHRO];
  Double_t eff_dis70[NCHRO];
Double_t eff_dis100[NCHRO];
  Double_t eff_dis150[NCHRO];
  Double_t eff_dis200[NCHRO];
  Double_t eff_dis250[NCHRO];
  Double_t eff_dis300[NCHRO];
  Double_t eff_dis400[NCHRO];
  Double_t eff_dis500[NCHRO];
  Double_t eff_dis800[NCHRO];
  Double_t eff_dis1000[NCHRO];
  Double_t eff_dis2000[NCHRO];
  Double_t eff_dis5000[NCHRO];
  Double_t eff_dis10000[NCHRO];
  Double_t t_level10[NCHRO];
  Double_t t_level15[NCHRO];
  Double_t t_level18[NCHRO];
  Double_t t_level20[NCHRO];
  Double_t t_level30[NCHRO];
  Double_t t_level40[NCHRO];
  Double_t t_level50[NCHRO];
  Double_t t_level60[NCHRO];
  Double_t t_level80[NCHRO];
  Double_t t_level100[NCHRO];
  Double_t t_level200[NCHRO];
  Double_t t_level300[NCHRO];

  Double_t t_corr50[NCHRO];

  Double_t trail_t_level4[NCHRO];
  Double_t trail_t_level5[NCHRO];
  Double_t trail_t_level10[NCHRO];
  Double_t trail_t_level15[NCHRO];
  Double_t trail_t_level18[NCHRO];
  Double_t trail_t_level20[NCHRO];
  Double_t trail_t_level30[NCHRO];
  Double_t trail_t_level40[NCHRO];
  Double_t trail_t_level80[NCHRO];
  Double_t trail_t_level100[NCHRO];

  Double_t t_frac05[NCHRO];
  Double_t t_frac10[NCHRO];
  Double_t t_frac15[NCHRO];
  Double_t t_frac20[NCHRO];
  Double_t t_frac25[NCHRO];
  Double_t t_frac30[NCHRO];
  Double_t t_frac50[NCHRO];
  Double_t t_frac70[NCHRO];
  Double_t t_frac80[NCHRO];
  Double_t t_frac90[NCHRO];
  Double_t t_rms3[NCHRO];
  Double_t t_rms5[NCHRO];
  Double_t trail_t_frac30[NCHRO];
  Double_t trail_t_frac50[NCHRO];
  Double_t trail_t_frac70[NCHRO];
  Double_t trail_t_frac90[NCHRO];
  Double_t t_max[NCHRO];
Double_t t_zero[NCHRO];
Double_t fitchi[NCHRO];
Double_t ampl_chi2[NCHRO];
Double_t rise_lin0[NCHRO];
Double_t rise_lin1[NCHRO];
Double_t rise_lin_chi2[NCHRO];
  Double_t rise_exp0[NCHRO];
  Double_t rise_exp1[NCHRO];
Double_t rise_exp_chi2[NCHRO];
  Double_t AMax[NCHRO];
  Double_t TMax[NCHRO];
Double_t NMax[NCHRO];
int polarity[NCHRO];


  // Analysis branches

Double_t TotCharge;
Int_t MaxEvt;
Int_t CS;
Int_t eff;
Double_t noise=3.;
Double_t bar0;
Double_t bar1;
Double_t A0;
Double_t A1;
Double_t delta0;
Double_t delta1;
Double_t chinorm0;
Double_t chinorm1;


  // Double_t gum_amp[MAX_CAMP];


Int_t nchro = NCHRO;  //Maximum number of Read out channels, to be filled manually in the InputFolder  
bool doFFT = false;  
bool showFFT = false;

std::string pip;
std::string InputNAME[100];
int i;
int j;
int k;
int np;
int nprec;


Double_t hv;
Double_t hvcorr;
Double_t temp;
Double_t pres;

Int_t nrun=0;
std::string String_nrun;
Int_t ntrig=0;
Int_t event=0;


std::ostringstream leaf(std::ostringstream::out);
std::ostringstream leafl(std::ostringstream::out);

  // Program variables
std::ostringstream datafile(std::ostringstream::out);

  float Tpeak = 0;
  float BallDef = 0;

int nmedia[NCHRO];
  double tempsum;
  int running=1;
  int emptycount=0;
  double time0temp;
  int segsize;
  char number[5];
  char date[12];
  char timeofday[8];
  double FFT_ratio;



//Convert little endian int to big endian int
unsigned int intToBigEndianess(const unsigned int x)
{
  return  ( x >> 24 ) |  // Move first byte to the end,
    ( ( x << 8 ) & 0x00FF0000 ) | // move 2nd byte to 3rd,
    ( ( x >> 8 ) & 0x0000FF00 ) | // move 3rd byte to 2nd,
    ( x << 24 ); // move last byte to start.
}


//returns the float value corresponding to a given bit represention of a scalar int value or vector of int values
float intBitsToFloat(const int x)
{
  union 
  {
    float f;  // assuming 32-bit IEEE 754 single-precision
    int i;    // assuming 32-bit 2's complement int
  } 
  u;
  
  u.i = x;
  return u.f;
}


//Always run background function before doing further analyses

//Compute average amplitude for background
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, max background amplitude)
void Background(Int_t camp, Double_t amp[], Double_t timeS[],Double_t t0, Double_t t_bck, Double_t *bck, Double_t *maxbck, Double_t *rmsbck)
{
  double max=-1000000;
  float sumamp=0;
  float sumampsq=0;
  int pointsavg=0;
  float meanbck;

  *bck=0;
  for(int j=0;j<camp;j++)
    {
      if((timeS[j]>t0 && timeS[j]<t_bck))
	{
	  pointsavg++;
	  sumamp=sumamp+amp[j];
	  sumampsq=sumampsq+amp[j]*amp[j];
	}
      meanbck=sumamp/float(pointsavg);
      *rmsbck=sqrt((1./float(pointsavg))*(sumampsq-float(pointsavg)*meanbck*meanbck));
    }

  *bck=meanbck;
  for(int j=0;j<camp;j++)
    {
      if((timeS[j]>t0 && timeS[j]<t_bck))
	{
	  if(std::fabs(amp[j]-meanbck)>max)
	    max=amp[j];
	}
    }

  *maxbck=max;
  //hlite.Reset();
  return;
}


//Inversione degli impulsi negativi
void Inversion(Int_t camp, Double_t neg_amp[])
{
  for(int j=0;j<camp;j++)
	{
	  neg_amp[j]=-neg_amp[j];
        }
}


//DEVELOPEMENT
// Compute maximum amplitude and peak time for waveform using a gaussian fit to the peak
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse maximum amplitude (averaged), pulse maximum amplitude)
void Amplitudes(Int_t camp, Double_t amp[], Double_t timeS[],  Double_t bck, Double_t t_bck,Double_t t_pul, Double_t *pamp, Double_t *tamp, Double_t *chi)
{
  
  Int_t np=0;
  //  int cont=0;
  int j_max= 0;
  int j_start = 0;
  int n_j = 0;
  // Limit of parabola fitting around j_max
  float AmpFrac = 0.8;
  // double tempmean;
  // double tempsig;
  double_t   max=-100000;
  double_t   tmax=-100000;

  //Maximum calculation
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  np++;
	  if(amp[j]-bck>max)
	    {
	      max=amp[j]-bck;
	      tmax=(timeS[j+1]+timeS[j])/2.;
	      j_max=j;
	    }
	}
    }

  //Fit amplitude calculation
  // Find limits at AmpFrac
  
  //from the maximum backwards 
  for(int j=j_max;j>0;j--)
    {
      if(amp[j]-bck<AmpFrac*(max))
	{
	  j_start=j;
	  break;
	}
    }

  //from the minimum limit, up the maximum, and then down till the limit on the falling edge  
  for(int j=j_start+1;j<camp-1;j++)
    {
      n_j++;
      if(amp[j]-bck<AmpFrac*(max))
	{
	  break;
	}
    }

  
  //  std::cout << j_start << " " <<j_max << " Number of points in the fit= " << n_j << std::endl; 
  
  double temp0;
  double temp1;
  double temp2;
  TF1 *f1 = new TF1("gausfit","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",timeS[j_start-1],timeS[j_start+n_j+1]);
  TGraph g(n_j);
  for(int j=0;j<n_j;j++)
    {
      g.SetPoint(j,timeS[j+j_start],amp[j+j_start]-bck);
    }

  temp0=max;
  temp1=tmax;
  temp2=(timeS[j_start]-timeS[j_start+n_j])/2.;
  f1->SetParameters(temp0,temp1,temp2);
  g.Fit("gausfit","QN","",timeS[j_start]-1,timeS[j_start+n_j-1]+1);
  //  g.Fit("gausfit","QN","",timeS[j_start]-1,timeS[j_start+n_j-1]+1);
  *pamp=f1->GetParameter(0);
  *tamp=f1->GetParameter(1);

  *chi=f1->GetChisquare()/f1->GetNDF();
  f1->Delete();

  return;
}

//Amplitude without the Fit
void Amplitudes_NF(Int_t camp, Double_t amp[], Double_t timeS[],  Double_t bck, Double_t t_bck,Double_t t_pul, Double_t *max,Double_t *tmax)
{
  

  *max=-100000;

  //Maximum calculation
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  np++;
	  if(amp[j]-bck>*max)
	    {
	      *max=amp[j]-bck;
	      *tmax=(timeS[j+1]+timeS[j])/2.;
	    }
	}
    }
  return;
}



//DEVELOPEMENT
// Compute maximum, amplitude and peak time for waveform
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse maximum amplitude (averaged), pulse maximum amplitude)
void Neg_Amplitudes(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck, Double_t *pamp, Double_t *pmax, Double_t *tmax)
{
  *pamp = 0;
  Int_t np=0;
  // int cont=0;
  double maxtemp=-100000;
  //  double tempmean;
  //  double tempsig;
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  np++;
	  if(-(amp[j]-bck)>maxtemp)
	    {
	      maxtemp=-(amp[j]-bck);
	      *tmax=(timeS[j+1]+timeS[j])/2.;
	    }
	}
    }
  *pmax=maxtemp;

  /*
  TGraph *g = new TGraph(np);
  TF1 *f= new TF1("gauss","[0]*exp(-((x-[1])*(x-[1]))/(2.*[2]*[2]))",t_bck,t_pul);
  f->SetParameters(*pmax,(t_pul+t_bck)/2.,(t_pul-t_bck)/2.);
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  cont++;
	  g->SetPoint(cont,timeS[j],amp[j]-bck);
	}
    }
  g->Fit("gauss","QN","",t_bck,t_pul);
  
  tempmean=f->GetParameter(1);
  tempsig=f->GetParameter(2);
  g->Fit("gauss","QN","",tempmean-tempsig,tempmean+tempsig);
  
  *pamp=f->GetParameter(0);
  */
  return;
}

// Calibration using curve [0](1-e^(x/[1]))
// Set the function with the proper front end each time
// (Parameter zero, parameter 1, amplitude of amplified pulse, corresponding charge/amplitude of input pulse, error on input charge/amplitude, error on output amplitude - must be known)
void Calibration(Double_t cal_low_par0,Double_t cal_low_par1,Double_t cal_high_par0,Double_t cal_high_par1, Double_t max,Double_t *charge, Double_t *chargeerr, Double_t amperr)
{
  TF1 f_low("cal_low","[0]*log(1.-x/[1])",0.,cal_low_par1-0.1);
  TF1 f_high("cal_high","[0]*log(1.-pow((x/[1]),2))",0.,cal_high_par1-0.1);
  f_low.SetParameters(cal_low_par0,cal_low_par1);
  f_high.SetParameters(cal_high_par0,cal_high_par1);
  if(max<444)
    {
      *charge=f_low(max);
      *chargeerr=(-cal_low_par0/(cal_low_par1-(max)))*amperr;
    }
  
  if(max>=444 && max<cal_high_par1)
    {
      *charge=f_high(max);
      *chargeerr=(-cal_high_par0)*(1./(1.-pow(1.-max/cal_high_par1,2)))*2.*(1.-max/cal_high_par1)*amperr;
    }
  if(max>=cal_high_par1)
    {
      *charge=1200.;
      *chargeerr=1000.;
    }
  
  return;
}

// Calibration using curve [0](1-e^(x/[1]))
// Set the function with the proper front end each time
// (Parameter zero, parameter 1, amplitudesof amplified pulse, corresponding charge/amplitude of input pulse, error on input charge/amplitude, error on output amplitude - must be known)
void Calibration_ip(Double_t cal_A0,Double_t cal_x0,Double_t cal_B0,Double_t max,Double_t *charge, Double_t *chargeerr, Double_t amperr)
{
  TF1 f("cal_ip","[1]+pow([0]/(x-[2]),0.25)",-10.,900);
  f.SetParameters(cal_A0,cal_x0,cal_B0);
  
  *charge=f(max);
  *chargeerr=-0.25*(pow(cal_A0/(max-cal_B0),1.25))*amperr/cal_A0;
    
   
  if(max>=cal_B0)
    {
      *charge=1500.;
      *chargeerr=1000.;
    }
  
  return;
}


// Integrates the pulse to give charge - readout impedance must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse charge, total charge, readout resistance)
void Charges(Int_t camp, Double_t amp[], Double_t timeS[], Double_t bck, Double_t t_bck, Double_t t_pul, Double_t *area, Double_t *charge, double R)
{
  double tempminch=0;
  //  double tempmaxch=0;
  // double tempmintotch=0;
  // double tempmaxtotch=0;
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_bck && timeS[j]<t_pul)	
	  tempminch += (amp[j+1]+amp[j])/2-bck;
    }
      *charge=tempminch*(timeS[1]-timeS[0])/R;
      *area=tempminch*(timeS[1]-timeS[0]);
  return;
}


//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, threshold to comupute efficiency, efficiency)
void Efficiencies(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck,Double_t threshold, Double_t *eff)
{
  *eff=0;
  for(int j=0;j<camp;j++)
    {
      if(timeS[j]<t_bck)
	continue;
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  if((amp[j]-bck-threshold)>0 && threshold>0)
	    {
	      *eff=1;
	      break;
	    }
	  if((amp[j]-bck-threshold)<0 && threshold<0)
	    {
	      *eff=1;
	      break;
	    }
	}
    }
  return;
}

// Time at which a fixed threshold is passed. Makes a linear fit on 6 points and inverts the equation
// time  = (threshold-q)/q
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, threshold to measure time, time at which the threshold is passed)
void Thrtime_old(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck,Double_t threshold, Double_t *thrtime)
{
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=0;j<camp;j++)
    {
      if(timeS[j]<t_bck)
	continue;
      if(timeS[j]>t_pul)
	continue;
      if(timeS[j]>t_bck && timeS[j]<t_pul)
	{
	  //	  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
	  if((amp[j]-bck-threshold)>0 && threshold>0)
	    {
	      //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	      //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	      //std::cout << "Der = " << Der << std::endl;
	      //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //   std::cout << "th1 = " << *thrtime << std::endl;
	      for ( int kk = 0;kk<npoints;kk++)
		{
		  g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
		  //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
		  //  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
		}

	      break;

	    }
	  //	  if((amp[j]-bck-threshold)<0 && threshold<0)
	  // {
	  //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
	}
    }
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
    //    Der= lin->GetParameter(1);
  if (lin->GetParameter(1)>0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
  // std::cout << "th1bis = " << *thrtime << std::endl;
  return;
}

// Time at which a threshold relative to the amplitude of the pulse is passed. Amplitude of the pulse must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, relative threshold to measure time, time at which the threshold is passed)
void ConstFractime(Double_t amp[], Double_t timeS[], Double_t bck,Double_t threshold, Double_t max, Int_t NMax,Double_t *thrtime)
{
  
  
  *thrtime=0;

  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  //  *thrtime=0;
  //  float Der = 0;
  //  for(int j=0;j<camp;j++)
    for(int j=NMax;j>1;j--)
    {
      //      if(timeS[j]>t_max || timeS[j]<t_bck )
      //	continue;
      // std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],(amp[j-npoints/2+kk]-bck)/max);
	      //  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << (amp[j-npoints/2+kk]-bck)/max << " " << threshold << std::endl;
	    }

	  break;
	  
	}
	
    }
  //  std::cout << "th1 = " << *thrtime << std::endl;
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  if (lin->GetParameter(1) !=0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
    
  //std::cout << "th2 = " << *thrtime << std::endl;
  g->Delete();
  lin->Delete();
  return;
}

// Time at which a threshold relative to the amplitude of the pulse is passed. Amplitude of the pulse must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, relative threshold to measure time, time at which the threshold is passed)
void Zerofittime(Int_t camp, Double_t amp[], Double_t timeS[], Double_t t_bck, Double_t t_pul, Double_t bck,Double_t threshold1, Double_t threshold2, Double_t t_max,Double_t *thrtime)
{
  t_pul = 0;
  *thrtime=0;
  double N=0;
  double sumV=0;
  double sumT=0;
  double sumVT=0;
  double sumTT=0;
  for(int j=camp-1;j>=0;j--)
    {
      if(timeS[j]>t_max)
	continue;
      if(timeS[j]<t_bck)
	continue;
      if(amp[j]-bck>threshold1 && amp[j]-bck<threshold2 && timeS[j]<t_max && N<3)
	{
	  N=N+1;
	  sumV=sumV+amp[j]-bck;
	  sumT=sumT+timeS[j];
	  sumVT=sumVT+(amp[j]-bck)*timeS[j];
	  sumTT=sumTT+timeS[j]*timeS[j];
	}
    }
  if(N>0)
    {
      double Vm=(N*sumVT-sumV*sumT)/(sumTT-sumT*sumT);
      double V0=(sumV-Vm*sumT)/N;
      *thrtime=-V0/Vm;
    }
  if(N<2)
    *thrtime=-99;
  return;
}

// Time at which a fixed threshold is passed (descending)
void Trailtime(Int_t camp, Double_t amp[], Double_t timeS[], Int_t NMax, Double_t bck,Double_t threshold, Double_t *thrtime)
{
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;

  for(int j=NMax;j<camp;j++)
    {
      if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	{
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
	    }
	  
	  break;
	  
	}
    }
    
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");

  if (lin->GetParameter(1)<0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
  //  std::cout << "th1bis = " << *thrtime << std::endl;
  return;
}

// Time at which a fixed threshold is passed (going to the maximum)
// The search starts from the maximum and goes backwards
void Thrtime( Double_t amp[], Double_t timeS[], Int_t NMax, Double_t bck,Double_t threshold, Double_t *thrtime)
{
  //  std::cout << "NMAX = " << NMax << std::endl;
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j>0;j--)
    {
     
      //std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	{
	  //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	  //std::cout << "Der = " << Der << std::endl;
	  //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //	       std::cout << "th1 = " << timeS[j] << std::endl;
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
	      //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
	      //	  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
	    }
	  
	  break;
	  
	}
      //	  if((amp[j]-bck-threshold)<0 && threshold<0)
      // {
      //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
    }
    
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  //    Der= lin->GetParameter(1);
  // if (lin->GetParameter(1)<0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
    //else
    // *thrtime = 0;
    //    std::cout << "th1bis = " << *thrtime << std::endl;
  return;
}


// Time at which a constant fraction threshold is passed (descending)
//(record lenght, amplitude vector, time vector, average background amplitude, relative threshold to measure time, starting time for the search (example: time of the matimum), time of the trailing edge, isteresys with respect to the leading edge)
void TrailConstFractime(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t threshold, Double_t max, Int_t NMax,Double_t *trailtime)
{
  *trailtime=0;
  /*
  for(int j=NMax;j<camp;j++)
    {
      
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  *trailtime=(timeS[j]+timeS[j-1])/2.;
	  break;
	}
      if(((amp[j]-bck)/max-threshold)>0 && threshold<0)
	{
	  *trailtime=(timeS[j]+timeS[j-1])/2.;
	  break;
	}
      
    }
*/
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  //  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j<camp;j++)
    {
      //	  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],(amp[j-npoints/2+kk]-bck)/max);
	      //		  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << (amp[j-npoints/2+kk]-bck)/max << " " << threshold << std::endl;
	    }
	  
	  break;
	  
	}
      
    }
    // std::cout << "th1 = " << *trailtime << std::endl;

    // TF1 *quad = new TF1("quad","[0]+[1]*x+[2]*x*x");
    TF1 *lin = new TF1("lin","[0]+[1]*x");
    g->Fit("lin","QN","goff");
    //g->Fit("lin");
    // g->Fit("quad");
  if (lin->GetParameter(1) !=0)
    *trailtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *trailtime = 0;
    
  //  std::cout << "th2 = " << *trailtime << std::endl;
  g->Delete();
  lin->Delete();

  
  return;
}

// Says if the event suffered from saturation of the waveform readout or not. Must be calibrated each time to evaluate the sensitivity of the routine
//(record lenght, amplitude vector, time vector,saturation flag)
void Saturation(Int_t camp, Double_t amp[], Double_t timeS[], Double_t bck, Int_t *satur)
{
  *satur=0;
  timeS = 0;
  int cont=0;
  int tempcont=0;
  double value=100000;
  for(int j=1;j<camp;j++)
    {
      if(fabs(amp[j]-amp[j-1])>0.0001)
	{
	  tempcont=0;
	}
      if(fabs(amp[j]-amp[j-1])<0.0001)
	{
	  tempcont++;
	  if(tempcont>cont)
	    {
	      cont=tempcont;
	      value=amp[j];
	    }
	}
    }
  if(fabs(value-bck)>1. && cont>int(camp/500)+5)
    *satur=1;
  return;
}

// Computed the charge centroid using the A0/cosh((x-x0)/delta0) function.
//(Charge, Error of charge, number of RO channels, distribution width, distribution width 2nd method, Distribution const, Distribution const 2nd method, Centroid , Centroid 2nd method , ChiSquare, ChiSquare 2nd method,pitch, Cluster size, threshold)
void Centroid(Double_t ch[], Double_t cherr[],Int_t nchro, Double_t *delta0, Double_t *delta1, Double_t *A0, Double_t *A1, Double_t *x0, Double_t *x1, Double_t *chinorm0, Double_t *chinorm1, Double_t pitch, Int_t CS, Double_t threshold)
{
  TF1 cint("cint",Form("[0]*atan(exp((x+%f-[1])/[2]))-[0]*atan(exp((x-%f-[1])/[2]))",pitch/2.,pitch/2.),-100.,100.);
  TF1 cdis("cdis","[0]*(1./cosh((x-[1])/[2]))",-100.,100.);
  TF1 cintf("cintf",Form("[0]*atan(exp((x+%f-[1])/[2]))-[0]*atan(exp((x-%f-[1])/[2]))",pitch/2.,pitch/2.),-100.,100.);
  TGraphErrors gbar(CS);
  double tmpmax=-100;
  double tmpbar=0;
  int k=0;
  *delta0=-1;
  *delta1=-1;
  *A0=-1;
  *A1=-1;
  *x0=-10;
  *x1=-10;
  *chinorm0=-1;
  *chinorm1=-1;
  cintf.FixParameter(2,3.9);
  if(CS>2)
    {
      for(int i=0;i<nchro;i++)
	{
	  if(ch[i]>tmpmax)
	    {
	      tmpmax=ch[i];
	      tmpbar=double(i)*pitch+pitch/2.;
	    }
	}
      
      cint.SetParameters(tmpmax,tmpbar,3.5);
      cint.SetParLimits(1,-pitch*5.,double(nchro)*pitch+pitch*5.);
      cintf.SetParameters(tmpmax,tmpbar,3.9);
      cintf.FixParameter(2,3.9);
      cintf.SetParLimits(1,-pitch*5.,double(nchro)*pitch+pitch*5.);

      for(int i=0;i<nchro;i++)
	{
	  if(ch[i]>threshold)
	    {
	      k++;
	      gbar.SetPoint(k,double(i)*pitch+pitch/2.,ch[i]);
	      gbar.SetPointError(k,0.,cherr[i]);
	    }
	}
      
      gbar.Fit("cint","QN","",0.,nchro*pitch+pitch/2.);
      *delta0=cint.GetParameter(2);
      *x0=cint.GetParameter(1);
      *A0=cint.GetParameter(0);
      *chinorm0=cint.GetChisquare()/cint.GetNDF();
      gbar.Fit("cintf","QN","",0.,nchro*pitch+pitch/2.);
      *delta1=cintf.GetParameter(2);
      *x1=cintf.GetParameter(1);
      *A1=cintf.GetParameter(0);
      *chinorm1=cintf.GetChisquare()/cintf.GetNDF();
    }
  if(CS==2)
    {
      *x0=0;
      float tempnum=0;
      float tempden=0;
      float bar0temp=0;
      float tempor=0;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      bar0temp=bar0temp+(float(j)*pitch+pitch/2.);
	    }
	}
      bar0temp=bar0temp/2.;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      tempnum=tempnum+ch[j]*(float(j)*pitch+pitch/2.-bar0temp);
	      tempden=tempden+ch[j];
	      //cout << FEcharge[j] << "\t";
	    }
	}
      tempor=tempnum/tempden;
      *x0=bar0temp+7.5*tempor;  //7.5 conversion factor from charge centroid to mm
    }
  if(CS==1)
    {
      *x0=0;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      *x0=*x0+(float(j)*pitch+pitch/2.);
	    }
	}
    }
  return;
}


//  Computes the slope of the rise of a signal as an exponential
void riseexpfit(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t t_low,  Double_t t_high, Double_t *rise_exp0, Double_t *rise_exp1, Double_t *chi)
{
  TGraph *g = new TGraph();
  int npoints=0;
  int nnpoints=0;
  int Nlow=0;

  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)  continue;
      if(timeS[i]>t_high) break;
      npoints++;
      if (npoints == 1) 
	{
	  Nlow = i;
	}
    }

  g->Set(npoints);
  nnpoints=0;
  for(int i=Nlow;i<Nlow+npoints-1;i++)
    {
            g->SetPoint(nnpoints,timeS[i]-t_low,amp[i]-bck);
	    //      g->SetPoint(nnpoints,timeS[i],amp[i]-bck);
	    //  std::cout << timeS[i]<< " " << amp[i]-bck << std::endl;
      nnpoints++;
    }
  
  TF1 *exp = new TF1("exp","exp([0])*exp([1]*x) +[2]");
  

  g->Draw();


  exp->SetParLimits(1, 0.,5.);
  //  exp->SetParLimits(2, -5.,5.);
  exp->SetParameters(-1,1.,bck);
  
  g->Fit("exp","QN","goff",0.,t_high-t_low);
  // g->Fit("exp"," "," ",0.,t_high-t_low);
  *rise_exp0=exp->GetParameter(0);
  *rise_exp1=exp->GetParameter(1);
 
  *chi=exp->GetChisquare()/exp->GetNDF();

  //  std::cout << " rise = " << *rise_exp1 << " chi2=  " << *chi << std::endl;
  exp->Delete();
  g->Delete();

  return;
}


void riselinfit(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t t_low,  Double_t t_high, Double_t *rise_lin0, Double_t *rise_lin1, Double_t *chi)
{
  TGraph *g = new TGraph();
  int npoints=0;
  int nnpoints=0;
  int Nlow=0;
  double_t amp_low =0;
  double_t amp_high =0;
  //  t_low = t_high-0.2;
  // std::cout << " t low " << t_low << " t_high " << t_high << std::endl;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)  continue;
      if(timeS[i]>t_high) break;
      npoints++;
      if (npoints == 1) 
	{
	  Nlow = i;
	}
    }
  amp_low = amp[Nlow]-bck;
  amp_high = amp[npoints]-bck;
  g->Set(npoints);
  nnpoints=0;
 
  for(int i=Nlow;i<Nlow+npoints-1;i++)
    {
      g->SetPoint(nnpoints,timeS[i],amp[i]-bck);
      nnpoints++;
    }
  
  TF1 *lin = new TF1("lin","[0]+[1]*x",t_low,t_high);
  // TF1 *exp = new TF1("exp","[0]*exp([1]*x)",t_low,t_high);

  lin->SetParameters(amp_low-t_low*(amp_high-amp_low)/(t_high-t_low),(amp_high-amp_low)/(t_high-t_low));
  g->Draw();
  g->Fit("lin","QN","goff",t_low,t_high);
  // g->Fit("lin"," ","g",t_low,t_high);
  *rise_lin0=lin->GetParameter(0);
  *rise_lin1=lin->GetParameter(1);
  *chi=lin->GetChisquare()/lin->GetNDF();
  lin->Delete();
  g->Delete();

  return;
}
//void CSA(Int_t camp, Double_t amp[], Int_t navg, Double_t m_amp[])
// {

  

// Computes the mobile average of the signal
// (record lenght,       amplitude vector, time vector, number of points for the mobile average (odd number), averaged amplitude vector)
void mobileAVG(Int_t camp, Double_t amp[], Int_t navg, Double_t m_amp[])
{
  float tempsum=0;
  for(int j=int((navg-1)/2);j<camp-int((navg-1)/2)-1;j++)
    {
      tempsum=0;
      for(int k=-int((navg-1)/2);k<int((navg-1)/2)+1;k++)
	{
	  tempsum=tempsum+amp[j+k];
	}
      m_amp[j]=tempsum/float(navg);
    }
  for(int j=0;j<int((navg-1)/2);j++)
    {
      m_amp[j]=m_amp[int((navg-1)/2)];
    }
  for(int j=camp-int((navg-1)/2)-1;j<camp;j++)
    {
      m_amp[j]=m_amp[camp-int((navg-1)/2)-2];
    }

  return;
}

void CSA( Double_t TimeUnit, Double_t TransImp, Double_t TauRC, Double_t TauFall, double_t TauRise, Int_t camp, Double_t Itot[],Double_t CSAout[])
{

  int IMaxSh = 60./TimeUnit; //TimeUnit in [ns]
  int IMax = 0;
  float Qtot = 0;
  float CSAMax = 0;
  float  PreAmp_Q[IMaxSh];
  

  //  std::cout << " Time Unit = " << TimeUnit << " Tau Fall " << TauFall << " Tau Rise " << TauRise << " TransImp = " << TransImp << std::endl;  


  for(int k=0;k<IMaxSh;k++)
    {
      PreAmp_Q[k]=0.0;
      CSAout[k] = 0.0;
    }

  
  for(int j=0;j<camp;j++)
    {
      if ( j>2 && Itot[j] == 0)
	{
	  IMax = j;
	  break;
	}
    }


 
  
  for(int i=1;i<IMax;i++)
    {
      
      //      PreAmp_Q[i] = Itot[i]*TimeUnit; // uA*ns ==> fC
	
      for (int ll = 0; ll<IMaxSh-i;ll++)  // valid only up to IMaxSh 
	{
	  PreAmp_Q[i+ll] += Itot[i]*TimeUnit*(1.-TMath::Exp(-ll*TimeUnit/TauRC)); // HS pag 101 transfer of charge to CSA
	  // if (i==2) std::cout << "itot = " << Itot[i] << " " << i+ll << " " << PreAmp_Q[i+ll] << " " << TMath::Exp(-ll*TimeUnit/TauRC) << " " << TimeUnit/TauRC << std::endl;
	}
      for (int ll = 0; ll<IMaxSh-i;ll++)  // valid only up to IMaxSh 
	{
	  CSAout[i+ll] += TauFall/(TauFall+TauRise)*(PreAmp_Q[i]-PreAmp_Q[i-1])*TransImp*
	    (TMath::Exp(-ll*TimeUnit/TauFall)-TMath::Exp(-ll*TimeUnit/TauRise)); // [Q] HS eq 4.3 This is the shaper	  
      	}
	  
      if (CSAout[i] > CSAMax)  CSAMax = CSAout[i]; 
      Qtot += Itot[i]*TimeUnit;
      //      std::cout << " i = " << 5+i*TimeUnit << "  Itot = " << Itot[i] <<  " CSA= " << CSAout[i] << " PreAmpQ= " <<  PreAmp_Q[i] << std::endl;
      
     }
  //  for (int ll = 0; ll<IMaxSh;ll++) CSAout[ll] *= TransImp*Qtot/CSAMax;

		     
  // for (int ll = 0; ll<IMaxSh;ll++)  // valid only up to IMaxSh 
  // std::cout <<  " Total Charge = " << Qtot << " Tpeak = " << Tpeak << " Ball. Deficit = " << BallDef << " CSAMax = "<< CSAMax << std::endl;       

  
  return;
}


//
//Calcola la derivata a quattro punti
// (record lenght, amplitude vector,reshaping vector)
void Derivative( Double_t TimeUnit,Int_t camp, Double_t amp[],Double_t resh[])
{
  Int_t Steps = 10;
  
  for(int j=Steps;j<camp-Steps;j++)
	{
	  resh[j]=(amp[j+10]-amp[j-10])/(2*Steps*TimeUnit);
	}
return;
}



// Generates FFT of a real vector (inreal) of lenght camp. Real part of output FFT is outre, immaginary part is outim. Module output is outmod.
void FFTrealtocomplex(Int_t camp, Double_t inreal[], Double_t outre[], Double_t outim[], Double_t outmod[])
{
  Double_t re;
  Double_t im;
  Double_t *in = new Double_t[camp];
  for(int i=0;i<camp;i++)
    {
      *(&in[i])=inreal[i];
    }
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &camp, "R2C");
  fftr2c->SetPoints(in);
  fftr2c->Transform();
  for (Int_t i=0; i<camp; i++)
    {
      fftr2c->GetPointComplex(i, re, im);
      outre[i]=re;
      outim[i]=im;
      outmod[i]=sqrt(re*re+im*im);
    }
  fftr2c->Delete();
  delete[] in;

  return;
}


// Generates the inverse FFT from a complex vector of lenght camp (two vectors: inre real part, inim immaginary part). The output vector is outreal.
void FFTcomplextoreal(Int_t camp, Double_t inre[], Double_t inim[], Double_t outreal[])
{
  Double_t *imma = new Double_t[camp];
  Double_t *real = new Double_t[camp];
  
  for(int i=0;i<camp;i++)
    {
      *(&real[i])=inre[i];
      *(&imma[i])=inim[i];
    }
  TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &camp, "C2R");
  fftc2r->SetPointsComplex(real,imma);
  fftc2r->Transform();
  for (Int_t i=0; i<camp; i++)
    {
      outreal[i]=fftc2r->GetPointReal(i)/double(camp);
    }
  fftc2r->Delete();
  delete[] imma;
  delete[] real;

  return;
}


// Redines coordinates correction distorsion of centroid. polinomial correction of coordinates (expansion at center of strips and compressio at borders)
//(centroid variable address, curve parameter)
void positionbias(Double_t *bar,Double_t sca)
{
  sca = 0;
  float barfz=*bar;
  float xloop=(1/100000.)*(int((100000.)*(barfz-4))%800000);
  
  *bar=((1/(70.+16.))*((70.)*(xloop-4)+pow((xloop-4),3))+(barfz-xloop-4)+8);

  return;
}


//Function to fit the waveform with Gumbel function
void Gumbel_fit(Int_t camp, Double_t amp[], Double_t timeS[],Double_t bck,Double_t t_low,Double_t t_up,Double_t t_amp, Double_t t_width, Double_t max, Double_t *bck_fit, Double_t *t_amp_fit, Double_t *t_width_fit, Double_t *amp_fit)
{
  TF1 *fgum = new TF1("gumbel","[0]+[3]*(exp(-(x-[1])/[2]+exp(-(x-[1])/[2])))",t_low,t_up);
  fgum->SetParameters(bck,t_amp,t_width,max);

  TGraph *g = new TGraph();
  int npoints=0;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)
	continue;
      if(timeS[i]>t_up)
	break;
      npoints++;
    }
  g->Set(npoints);
  npoints=0;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)
	continue;
      if(timeS[i]>t_up)
	break;
      g->SetPoint(npoints,timeS[i],amp[i]);
      npoints++;
    }

  g->Draw();
  g->Fit("gumbel","QN","goff",t_low,t_up);
  *bck_fit=fgum->GetParameter(0);
  *t_amp_fit=fgum->GetParameter(1);
  *t_width_fit=fgum->GetParameter(2);
  *amp_fit=fgum->GetParameter(3);

  g->Delete();
  fgum->Delete();

  return;
}

double string_to_double( const std::string& s )
 {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 } 




#endif
