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

int main()
{

  //  bool doFFT=false;  
  bool showFFT=true;
  
  //  int campfact=1;  //divider to sample rate for waveform memorization

  TFile *OutputFile = new TFile("RunXX.root","recreate");

  ifstream InputCARD("Input_Folder.txt");
  
  if(!InputCARD)
    {
      cout << "Error: could not find the InputCARD" << endl;
      return (0);
    }


  // read the active number of channels
  int ntmp;
  int np_Max = 10000;
  TRandom *xi = new TRandom();
  InputCARD >> pip >> ntmp;
  bool FWF2 = false;
  if (ntmp > 0)
    {
      cout << "LECROY data " << endl;
      nchro = ntmp;
    }
  
  else if (ntmp < 0)
    {
      cout << "WF2 data " << endl;
      nchro = -ntmp;
      FWF2 = true;
    }
  //Waveform variables



  //  const Int_t max_camp=20000;
  // const Int_t max_camprec=int(max_camp/campfact)+100;



  // Tree definition
  TTree *w = new TTree("Wave","Wave");
  w->SetDirectory(OutputFile);
  TTree *OutTree = new TTree("Analysis","Analysis");
  OutTree->SetDirectory(OutputFile);
 

  // Analysis branches
  OutTree->Branch("ntrig",&ntrig,"ntrig/I");
  OutTree->Branch("event",&event,"event/I");
  OutTree->Branch("nrun",&nrun,"nrun/I");
  OutTree->Branch("nchro",&nchro,"nchro/I");
  OutTree->Branch("pres",&pres,"pres/D");
  OutTree->Branch("temp",&temp,"temp/D");
  OutTree->Branch("noise",&noise,"noise/D");
  //  OutTree->Branch("maxstrip",&maxstrip,"maxstrip/I");
  // OutTree->Branch("maxamp",&maxamp,"maxamp/D");
  //  OutTree->Branch("p_charge",&p_charge,"p_charge/D");
  // OutTree->Branch("p_chargeerr",&p_chargeerr,"p_chargeerr/D");
  //  OutTree->Branch("CS",&CS,"CS/I");
  // OutTree->Branch("eff",&eff,"eff/I");
  // OutTree->Branch("bar0",&bar0,"bar0/D");
  // OutTree->Branch("bar1",&bar1,"bar1/D");
  // OutTree->Branch("A0",&A0,"A0/D");
  //OutTree->Branch("A1",&A1,"A1/D");
  // OutTree->Branch("delta0",&delta0,"delta0/D");
  // OutTree->Branch("delta1",&delta1,"delta1/D");
  // OutTree->Branch("hv",&hv,"hv/D");
  OutTree->Branch("hvcorr",&hvcorr,"hvcorr/D");
  //  OutTree->Branch("nchro",&nchro,"nchro/I");
  OutTree->Branch("t_bck",t_bck,"t_bck[nchro]/D");
  OutTree->Branch("t_pul",t_pul,"t_pul[nchro]/D");
  OutTree->Branch("max",maxS,"maxS[nchro]/D");
  OutTree->Branch("maxD",maxD,"maxD[nchro]/D");
  OutTree->Branch("area",area,"area[nchro]/D");
  OutTree->Branch("ampl",ampl,"ampl[nchro]/D");
  OutTree->Branch("ampl_chi2",ampl_chi2,"ampl_chi2[nchro]/D");
  OutTree->Branch("dVdt3070",dVdt3070,"dVdt3070[nchro]/D");
  OutTree->Branch("dVdt1030",dVdt1030,"dVdt1030[nchro]/D");
  OutTree->Branch("dVdt2080",dVdt2080,"dVdt2080[nchro]/D");
  OutTree->Branch("strip",strip,"strip[nchro]/I");
  OutTree->Branch("bck",bck,"bck[nchro]/D");
  OutTree->Branch("max_bck",max_bck,"max_bck[nchro]/D");
  OutTree->Branch("rms_bck",rms_bck,"rms_bck[nchro]/D");
  OutTree->Branch("max_fwd_bck",max_fwd_bck,"max_fwd_bck[nchro]/D");
  OutTree->Branch("rms_fwd_bck",rms_fwd_bck,"rms_fwd_bck[nchro]/D");


  OutTree->Branch("t_level10",t_level10,"t_level10[nchro]/D");
  OutTree->Branch("t_level15",t_level15,"t_level15[nchro]/D");
  OutTree->Branch("t_level18",t_level18,"t_level18[nchro]/D");
  OutTree->Branch("t_level20",t_level20,"t_level20[nchro]/D");
  OutTree->Branch("t_level30",t_level30,"t_level30[nchro]/D");
  OutTree->Branch("t_level40",t_level40,"t_level40[nchro]/D");
  OutTree->Branch("t_level50",t_level50,"t_level50[nchro]/D");
  OutTree->Branch("t_level60",t_level60,"t_level60[nchro]/D");
  OutTree->Branch("t_level80",t_level80,"t_level80[nchro]/D");
  OutTree->Branch("t_level100",t_level100,"t_level100[nchro]/D");
  OutTree->Branch("t_level200",t_level200,"t_level200[nchro]/D");
  OutTree->Branch("t_level300",t_level300,"t_level300[nchro]/D");

  OutTree->Branch("trail_t_level10",trail_t_level10,"trail_t_level10[nchro]/D");
  OutTree->Branch("trail_t_level20",trail_t_level20,"trail_t_level20[nchro]/D");
  OutTree->Branch("trail_t_level30",trail_t_level30,"trail_t_level30[nchro]/D");
  OutTree->Branch("trail_t_level40",trail_t_level40,"trail_t_level40[nchro]/D");
  OutTree->Branch("trail_t_level80",trail_t_level80,"trail_t_level800[nchro]/D");
  OutTree->Branch("trail_t_level100",trail_t_level100,"trail_t_level100[nchro]/D");

  OutTree->Branch("t_frac05",t_frac05,"t_frac05[nchro]/D");
  OutTree->Branch("t_frac10",t_frac10,"t_frac10[nchro]/D");
  OutTree->Branch("t_frac15",t_frac15,"t_frac15[nchro]/D");
  OutTree->Branch("t_frac20",t_frac20,"t_frac20[nchro]/D");
  OutTree->Branch("t_frac25",t_frac25,"t_frac25[nchro]/D");
  OutTree->Branch("t_frac30",t_frac30,"t_frac30[nchro]/D");
  OutTree->Branch("t_frac40",t_frac40,"t_frac40[nchro]/D");
  OutTree->Branch("t_frac50",t_frac50,"t_frac50[nchro]/D");
  OutTree->Branch("t_frac60",t_frac60,"t_frac60[nchro]/D");
  OutTree->Branch("t_frac70",t_frac70,"t_frac70[nchro]/D");
  OutTree->Branch("t_frac80",t_frac80,"t_frac80[nchro]/D");
  OutTree->Branch("t_frac90",t_frac90,"t_frac90[nchro]/D");
  OutTree->Branch("t_rms3",t_rms3,"t_rms3[nchro]/D");
  OutTree->Branch("t_rms5",t_rms5,"t_rms5[nchro]/D");
  OutTree->Branch("trail_t_frac30",trail_t_frac30,"trail_t_frac30[nchro]/D");
  OutTree->Branch("trail_t_frac50",trail_t_frac50,"trail_t_frac50[nchro]/D");
  OutTree->Branch("trail_t_frac70",trail_t_frac70,"trail_t_frac70[nchro]/D");
  OutTree->Branch("trail_t_frac90",trail_t_frac90,"trail_t_frac90[nchro]/D");
  OutTree->Branch("t_max",t_max,"t_max[nchro]/D");
  OutTree->Branch("t_zero",t_zero,"t_zero[nchro]/D");
  OutTree->Branch("rise_lin0",rise_lin0,"rise_lin0[nchro]/D");
  OutTree->Branch("rise_lin1",rise_lin1,"rise_lin1[nchro]/D");
  OutTree->Branch("rise_lin_chi2",&rise_lin_chi2,"rise_lin_chi2[nchro]/D");
  OutTree->Branch("rise_exp0",rise_exp0,"rise_exp0[nchro]/D");
  OutTree->Branch("rise_exp1",rise_exp1,"rise_exp1[nchro]/D");
  OutTree->Branch("rise_exp_chi2",rise_exp_chi2,"rise_exp_chi2[nchro]/D");
  // OutTree->Branch("year",&year,"year/D");
  // OutTree->Branch("month",&month,"month/D");
  // OutTree->Branch("day",&day,"day/D");
  // OutTree->Branch("hour",&hour,"hour/D");
  // OutTree->Branch("minute",&minute,"minute/D");
  // OutTree->Branch("sec",&sec,"sec/D");

 
  
  for(i=0;i<nchro;i++)
    {
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "camp_" << i;
      leafl << "camp_" << i <<"/I";
      OutTree->Branch(leaf.str().c_str(),&camp[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "camprec_" << i;
      leafl << "camprec_" << i <<"/I";
      OutTree->Branch(leaf.str().c_str(),&camprec[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "t_bck" << i;
      leafl << "t_bck" << i <<"/D";
      OutTree->Branch(leaf.str().c_str(),&t_bck[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "t_pul" << i;
      leafl << "t_pul" << i <<"/D";
      OutTree->Branch(leaf.str().c_str(),&t_pul[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "bck" << i;
      leafl << "bck" << i <<"/D";
      OutTree->Branch(leaf.str().c_str(),&bck[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "osc" << i;
      leafl << "osc" << i <<"/I";
      OutTree->Branch(leaf.str().c_str(),&osc[i],leafl.str().c_str());
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "chan" << i;
      leafl << "chan" << i <<"/I";
      OutTree->Branch(leaf.str().c_str(),&chan[i],leafl.str().c_str());
      
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "amp" << i;
      leafl << "amp" << i <<"[camprec_" << i << "]/D";
      OutTree->Branch(leaf.str().c_str(),&amprec[i][0],leafl.str().c_str());
      
      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "m_amp" << i;
      leafl << "m_amp" << i <<"[camprec_" << i << "]/D";
      OutTree->Branch(leaf.str().c_str(),&m_amprec[i][0],leafl.str().c_str());

      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
      leaf << "der_amp" << i;
      leafl << "der_amp" << i <<"[camprec_" << i << "]/D";
      OutTree->Branch(leaf.str().c_str(),&d_amprec[i][0],leafl.str().c_str());


  //OutTree->Branch("gum_amp",gum_amp,"gum_amp[camp_1]/D");

      if(showFFT)
	{
	  
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "FFT_abs" << i;
	  leafl << "FFT_abs" << i <<"[camprec_" << i << "]/D";
	  OutTree->Branch(leaf.str().c_str(),&FFT_abs[i][0],leafl.str().c_str());
	  
	  
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "FFT_real" << i;
	  leafl << "FFT_real" << i <<"[camprec_" << i << "]/D";
	  OutTree->Branch(leaf.str().c_str(),&FFT_real[i][0],leafl.str().c_str());
	  
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "FFT_comp" << i;
	  leafl << "FFT_comp" << i <<"[camprec_" << i << "]/D";
	  OutTree->Branch(leaf.str().c_str(),&FFT_comp[i][0],leafl.str().c_str());
	}


      
    }

  OutTree->Branch("time",timerec,"time[camprec_0]/D"); 
  OutTree->Branch("freq",freq,"freq[camp_0]/D");
  

  
  cout << "Branch settings done" << endl;

  //  TF1 f_gum("gumbel","[0]+[3]*(exp(-(x-[1])/[2]+exp(-(x-[1])/[2])))",0,200);

  float TransImp = 9.8;   
  float TauRC = .36;   
  float TauFall = 5.1/2.2;
  float TauRise = 0.5/2.2;
  //  cout << TauRise << " TFall = " << TauFall << endl; 
  event = 0;
  //Loop on different runs
  while(1)
    {
      pip="";
      if(InputCARD.eof())
	break;
      
      //Read run variables
      InputCARD >> pip >> nrun >> pip >> ntrig >> pip >> nchro >> pip >> MaxEvt;
      if (nrun == -1 ) break;
      cout << endl << endl << "Run " << nrun << "\t " << nchro << " channels" << " Events to be analized = " << MaxEvt - ntrig << endl;
      //Note: ntrig is the number of the FIRST trigger of the run
      
      //      cout << DecimalToBinaryString(nrun) << endl;
      //Inizialization of data run variables
      ntrig=ntrig-1;
     
      ostringstream convert;
      convert << nrun ;
      String_nrun = convert.str();

      //  cout << " String nrun " << String_nrun << endl;
      

      //Inizialization of program run variables
      running=1;
      
      //Inizialization of channel variables for this run
      for(i=0;i<nchro;i++)
	{
	  InputCARD >> pip >> InputNAME[i] >> pip >>nmedia[i];  //NOTE: InputNAME is the complete path + file name without trigger number and final .txt
	  FreqCut[i] = 0. ;
	  if (nmedia[i]>0)
	    {
	      cout << endl << "Channel " << i << " File name:" << endl << InputNAME[i] << " averaging " << nmedia[i] << " points " <<endl;
	    }
	  else
	    {
	      FreqCut[i] = -nmedia[i]*1.e6;
	      cout << endl << "Channel " << i << " File name:" << endl << InputNAME[i] <<  " with Freq. Cut = " << FreqCut[i] <<endl;
	    }
	    
	  camp[i]=0;
	}
      string mystr;

      //      return 0;

     
	  
      TH1F * BkgHistCSA;
      TH1F * BkgHist;
      TFile Run4f("Run4_Bkg.root");
      TFile Run11f("Run11_500MHz_Bkg.root");
      
      BkgHistCSA =  (TH1F*)Run4f.Get("Bkg1");      
      BkgHist =  (TH1F*)Run11f.Get("Bkg1");
      int HistoChoice = 0;	  

      int NbinBkgCSA = BkgHistCSA->GetNbinsX();	        
      int NbinBkg = BkgHist->GetNbinsX();	        	 

      
	  /*
	  
	  stringstream ss1;
	  ss1 << "Bkg" <<  HistoChoice ;
	  string  fileNameT1 = ss1.str();
	  
	  
	  stringstream ss2;
	  ss2 << "Bkg" <<  HistoChoice ;
	  string  fileNameT4 = ss2.str();
*/

      string RunBit;
      string RunRoot;
      RunRoot = String_nrun[0];
      if (atoi(RunRoot.c_str())==9)
	{
	  cout << " Run number starts with 9: " << nrun << endl;
	  cout << " you requested input file shifting" << endl;
	  cout << " Are you sure?" << endl;
	}
      

      
      
	  //   Run11f.ls();
      // cout << "Bin Content = " << BkgHist->GetBinContent(2) << endl;
      //Loop on different triggers of the run
      //      cout << endl << "Event number: " << endl;
      emptycount=0;
      while(running)  //running is set on 1 as long as all the channels have data for a trigger event
	{
	  if (FWF2)
	    {
	      HistoChoice = 8*xi->Uniform()+1;
	      stringstream ss1;
	      ss1 << "Bkg" <<  HistoChoice ;
	      string  fileNameT1 = ss1.str();
	      BkgHist =  (TH1F*)Run11f.Get(fileNameT1.c_str());
	      
	      
	      stringstream ss2;
	      ss2 << "Bkg" <<  HistoChoice ;
	      string  fileNameT4 = ss2.str();
	      BkgHistCSA =  (TH1F*)Run4f.Get(fileNameT4.c_str());	  
	    }
	  ntrig++;
	  event++;
	  if (MaxEvt > 0 && ntrig >MaxEvt) goto Write;
	  if(event==1)
	    cout << "Event = " << event << " at trigger = " << ntrig << endl;
	  if(ntrig%200==0)
	    cout << "Event = " << event << " at trigger = " << ntrig << endl;

	  
	  //Fill each channel of the event
	  //	  cout << "File with " << nchro << " active channels" << endl;

	  for(i=0;i<nchro;i++)
	    {
	      ampl_chi2[i]=0;
	      AMax[i]=0;
	      NMax[i]=0;
	      RunBit = String_nrun[i+1];
	      datafile.str("");
	      datafile.clear();
	      //cout << "check 1" << endl;

	      // fix file number for certain runs 
	      if (atoi(RunRoot.c_str())==9 && atoi( RunBit.c_str()))
		ntrig ++;

	      
	      if (ntmp>0) sprintf(number,"%05d", ntrig);
	      if (ntmp<0) sprintf(number,"%d", ntrig);
	      if (nrun ==777) 	      datafile << InputNAME[i] << number << ".csv";  //Read file containing the event on that channel

	      else datafile << InputNAME[i] << number << ".txt";  //Read file containing the event on that channel
	      pip=datafile.str();

	      if (atoi(RunRoot.c_str())==9 && atoi( RunBit.c_str())) ntrig --;	      
	      if (event==1)
		cout << " First file for channel " << i << " is " << datafile.str() << endl;


	      
	      // cout << "datafile: " << datafile.str() << endl;
	      ifstream data(datafile.str().c_str());
	      if(!data)
		{
		  if(emptycount<1)
		    {
		      cout << "Data file " << datafile.str() << "  not found! Run " << nrun << "  event " << ntrig << endl;
		      cout << "End of run" << endl;
		      //cout << "File name: " << datafile.str() << endl << endl;
		      goto Write;
		    }
		  emptycount++;
		  running=0;
		  break;
		}
	      //cout << "check 1fin" << endl;
	      camp[i]=0;
	      tempsum=0;
	      np=0; 
	      nprec=0;
	    	      
	      std::string token;
	      if (!FWF2) // LECROY FORMAT
		{
		  getline (data,mystr); //cout << mystr << endl;
		  getline (data,mystr);//  cout << mystr << endl;
		  istringstream iss(mystr);
		  std::getline(iss, token, ',');
		  std::getline(iss, token, ',');
		  std::getline(iss, token, ',');
		  std::getline(iss, token, ',');	      
		  segsize = string_to_double( token );
		  
		  getline (data,mystr); // cout << mystr << endl;
		  getline (data,mystr); // cout << mystr << endl;

		  istringstream isss(mystr);
		  std::getline(isss, token, ':');
		  std::getline(isss, token, ':');
		  sec[i] = string_to_double( token );
		  
		  //		  day=double(int(date[0] - '0'))*10+double(int(date[1] - '0'));
		  // year=double(int(date[8] - '0'))*1000+double(int(date[9] - '0'))*100+double(int(date[10] - '0'))*10+double(int(date[11] - '0'));
		  // hour=double(int(timeofday[0] -'0')*10)+double(int(timeofday[1] - '0') *1);
		  // minute=double(int(timeofday[3] - '0')*10)+double(int(timeofday[4] - '0')*1);
		  // sec=double(int(timeofday[6] - '0')*10)+double(int(timeofday[7] - '0')*1);
		  //cout << "Segsize: " << segsize << endl;
		  //		  cout << "Sec: " << sec << endl;

		  if (i>0)
		    {
		      //		      if (fabs(sec[i]-sec[i-1])>0 && (fabs(sec[i]-sec[i-1])!= 59 ))
		      if (fabs(sec[i]-sec[i-1])>0)
			{
			  cout << "Sincronization error on event "  << ntrig << endl;
			  cout << "Acquisition time for channel " << i << " is " << sec[i] << " seconds" << endl;
			  cout << "Acquisition time for channel " << i-1 << " is " << sec[i-1] << " seconds" << endl; 
			  // cout << "The program stops " << endl;
			  continue;
			  //			  return(-1);
			}
		    }
		  getline (data,mystr); // cout << mystr << endl;
		  
		}

	      else  // WF2 FORMAT
		{
		  getline (data,mystr); // cout << mystr << endl;
		}
	    
	      int NoiseBase = 200*xi->Uniform();
	      // cout << " NoiseBase = " << NoiseBase << endl;
	      while(np<np_Max) //upper limit for safety
		//	      while(1)
		{
		  if(data.eof())
		    break;
		  
		  
		  if (!FWF2) // LECROY
		    {
		      getline (data,mystr); 
		      istringstream iss(mystr);		   
		      std::getline(iss, token, ',');
		      timeS[np]= string_to_double( token );
		      std::getline(iss, token, ',');
		      amp[i][np]= string_to_double( token );
		    }
		  else //WF2 format
		    {
		      //	      cout << np << endl;
		      if (i == 0 ) data >>  timeS[np] >> Itot[i][np] >> pip>> pip>> pip>> pip>> amp[i][np] >>pip >> pip;
		      if (i == 1 ) data >>  timeS[np] >> Itot[i][np] >> pip>> pip>> pip>> pip>> pip >>  amp[i][np] >> pip;
		      timeS[np] *=1.e-9;
		      //  if (i == 0 && np< NbinBkg-NoiseBase) amp[i][np] += BkgHist->GetBinContent(NoiseBase+np)*1.e-3;
		      // if (i == 1 && np< NbinBkgCSA-NoiseBase) amp[i][np] += BkgHistCSA->GetBinContent(NoiseBase+np)*1.e-3;

		    }
		  
		  //  if (ntrig == 0 ) cout << "Index = " << np << " Time = " <<  timeS[np]*1e9 << ", Amplitude = " << amp[i][np]*1e3 << endl;
		  if(np==0)
		    {
		      time0temp=timeS[np];   //Sets first sample time to zero
		    }
		  timeS[np]=(timeS[np]-time0temp)*1.e9; // time in ns
		  if (!FWF2) amp[i][np]=amp[i][np]*1.e3; // ampl. in mV
		  DT = (timeS[10]-timeS[9]); // delta T in nanosecond
		  //	  cout << DT << endl;
		  if (fabs(amp[i][np])>AMax[i]) 
		    {
		      AMax[i] = fabs(amp[i][np]); 
		      TMax[i] = timeS[np];
		      polarity[i] = amp[i][np]/fabs(amp[i][np]);
		      NMax[i] = np;
		    }
		  
		  np++;
		} // end of reading the file

	      //	      if (i == 0 ) BB;
	     
	      if (event == 1 || event == 50) np_Max = np-10; // this is a safety measure, in case the end of file is not found
	      //   cout << "np max = "  << np_Max << endl;
	      if (event == 2  && i == 0) cout << "Number of points: " << np  << " with a time step of " << DT << " [ns] " << endl;
	      if (event == 1  && i == 1)
		{
		  Tpeak = TauRise/(1.-TauRise/TauFall)*std::log(TauFall/TauRise);
		  BallDef = std::pow(TauRise/TauFall, TauRise/(TauFall-TauRise));
		  cout << " Ballistic deficit " << BallDef  << " Time Peak = " << Tpeak << " [ns] " << endl;

		}
	      if (event <20  ) cout << "Event: " << event << " channel " << i << " has maximum value of " << AMax[i] << " [mV] at  " <<  TMax[i] << " [ns] " << " on sample " << NMax[i] << " polarity = " << polarity[i] << endl;
	
	      if (!FWF2)
		{
		  t_bck[i] = (TMax[i]-20>0) ? TMax[i]-20 : 7 ;
		  t_pul[i] = (TMax[i]+30 < timeS[np-2]) ? TMax[i]+30 : timeS[np-2] ;
		}
	      else
		{
		  t_bck[i] = 0. ;
		  t_pul[i] = 15.;
		}
	      camp[i]=np-1;
	      //cout << "check 2" << endl;
	      // running average of nmedia point 
	      //cout << nmedia[i] << endl;
	      
	      //if we don't do a FFT, then we do an average

	      //	      float TauRise = pow(TauRiseCSA*TauRiseCSA+TauRC*TauRC,0.5);
	      
	      if (i == -11 )
		{
		  //DT, Tau... all in [ns]

		  CSA( DT, TransImp, TauRC, TauFall, TauRise, camp[i], &Itot[i][0],&m_amp[i][0]);

		  for (int il=0; il<camp[i] ;il++)	 cout << il << ", "  << amp[i][il] << ", " << m_amp[i][il] << endl;		  
		}

	      
	      if (FreqCut[i] == 0)
		{
		  mobileAVG(camp[i],&amp[i][0],nmedia[i],&m_amp[i][0]);
		}
	      else
		{
		  if (i==0)
		    for(int l=0;l<camp[i];l++)		      
		      freq[l]=(1.e9/DT)*(1./camp[i])*l; //http://it.mathworks.com/help/matlab/math/fast-fourier-transform-fft.html			  
		  
		  FFTrealtocomplex(camp[i],&amp[i][0],&FFT_real[i][0],&FFT_comp[i][0],&FFT_abs[i][0]);
		  for(int l=0;l<camp[i];l++)
		    {
		      if (freq[l]< FreqCut[i])
			{
			  FFT_real_cut[i][l] = FFT_real[i][l];
			  FFT_comp_cut[i][l] = FFT_comp[i][l];
			}
		      FFTcomplextoreal(camp[i],&FFT_real_cut[i][0],&FFT_comp_cut[i][0],&m_amp[i][0]);
		    }
		}
	      
	      
	      
	      
	      for(j=0;j<camp[i];j++)
	      	{
	      	  if(j%CAMPFACT==0) // option to undersample the data by CAMPFACT
		    {
		      timerec[nprec]=timeS[j]; // fills the time axis
		      amprec[i][nprec]=amp[i][j]; // fills the amplitude
		      m_amprec[i][nprec]=m_amp[i][j]; //fill the average amplitude
		      if(polarity[i]<0)
			{
			  amprec[i][nprec]=-amprec[i][nprec];
		      	  m_amprec[i][nprec]=-m_amprec[i][nprec];
		      	}
		      // cout << nprec << "\t" << timerec[nprec] << "\t" << amprec[i][nprec] << endl;
		      nprec++;
		    }
		}
	      
	      camprec[i]=nprec;
	      
	      // compute the derivative on the signal shape
	      Derivative( DT, camp[i],&m_amprec[i][0],&der_amp[i][0]);
	      

	
	      //	      continue;	      	      
	      //Place your single channel analysis here!
	      // analysis
	      // for(int j=0;j<camp[i];j++) amp[i][j]=polarity[i]*m_amp[i][j];
	      //  cout << "t_bck[" <<i <<"] = "<< t_bck[i] << endl;
 


	      
	      //calculate the baseline rms and amplitude after the signal, the from 10 to 15 ns after the max

	     
	      Background(&m_amprec[i][0], NMax[i]+10./DT,NMax[i]+15./DT,&bck[i],&max_fwd_bck[i],&rms_fwd_bck[i]);

	      //calculate the baseline rms and amplitude before the signal, from -10 to -5 ns before the signal

	      //	      Background(&m_amprec[i][0], NMax[i]-10./DT, NMax[i]-5./DT,&bck[i],&max_bck[i],&rms_bck[i]);
	      Background(&m_amprec[i][0], 0./DT, 5./DT,&bck[i],&max_bck[i],&rms_bck[i]);

	      if (ntmp<0)
		{
		  bck[i] = 0;
		  max_bck[i] = 0;
		  rms_bck[i] = 0;
		}
	      //	      std::cout << t_bck[i]-5 << " " << t_bck[i] << std::endl;
	      
	      //	      std::cout << "Time bckg = " <<  t_bck[i] << " - " << t_pul[i] << std::endl;
	      
	      Amplitudes(camp[i],&m_amprec[i][0],timeS,bck[i],t_bck[i],t_pul[i],&ampl[i],&t_max[i],&ampl_chi2[i]);	      

	      // std::cout << bck[i] << " Ch = " << i << " Ampl: " << ampl[i] << std::endl;
	      
	      Charges(&m_amprec[i][0],DT,NMax[i]-5./DT, NMax[i]+5./DT,&area[i],&totcha[i],50);
	      //  cout << "area = " << area[i] << " amp = " << ampl[i] << " ratio Q/A = " << area[i]/ampl[i] <<endl;

	      // Amplitudes_NF(camp[i],&m_amprec[i][0],timeS,bck[i],t_bck[i],t_pul[i],&ampl[i],&t_max[i]);

	      // time at a fixed threshold  
	      if ( ampl[i]< 10) continue;

	      Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],3.*rms_bck[i],&t_rms3[i]);
	      Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],5.*rms_bck[i],&t_rms5[i]);		      


	      if (ampl[i]>10) Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],10.,&t_level10[i]);
	      if (ampl[i]>20) Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],20.,&t_level20[i]);
	      if (ampl[i]>30) Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],30.,&t_level30[i]);
	      if (ampl[i]>40) Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],40.,&t_level40[i]);
	      if (ampl[i]>50) Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],50.,&t_level50[i]);
	      if (ampl[i]>60) Thrtime(&m_amprec[i][0],timeS,NMax[i],bck[i],60.,&t_level60[i]);		      

	      if (ampl[i]>10) Trailtime(camp[i],&m_amprec[i][0],timeS,NMax[i],bck[i],10.,&trail_t_level10[i]);
	      if (ampl[i]>20) Trailtime(camp[i],&m_amprec[i][0],timeS,NMax[i],bck[i],20.,&trail_t_level20[i]);
	      if (ampl[i]>30) Trailtime(camp[i],&m_amprec[i][0],timeS,NMax[i],bck[i],30.5,&trail_t_level30[i]);
	      if (ampl[i]>40) Trailtime(camp[i],&m_amprec[i][0],timeS,NMax[i],bck[i],40.,&trail_t_level40[i]);
	      if (ampl[i]>80) Trailtime(camp[i],&m_amprec[i][0],timeS,NMax[i],bck[i],80.,&trail_t_level80[i]);
	      if (ampl[i]>100) Trailtime(camp[i],&m_amprec[i][0],timeS,NMax[i],bck[i],100.,&trail_t_level100[i]);

	
	      
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.05,ampl[i],NMax[i],&t_frac05[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.1,ampl[i],NMax[i],&t_frac10[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.15,ampl[i],NMax[i],&t_frac15[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.2,ampl[i],NMax[i],&t_frac20[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.25,ampl[i],NMax[i],&t_frac25[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.3,ampl[i],NMax[i],&t_frac30[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.4,ampl[i],NMax[i],&t_frac40[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.5,ampl[i],NMax[i],&t_frac50[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.6,ampl[i],NMax[i],&t_frac60[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.7,ampl[i],NMax[i],&t_frac70[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.8,ampl[i],NMax[i],&t_frac80[i]);
	      ConstFractime(&m_amprec[i][0],timeS,bck[i],0.9,ampl[i],NMax[i],&t_frac90[i]);
	      //	      std::cout << "i = " << i << " t_max = " << t_max[i] << "  t@ 30% " << t_frac30[i] << " 50% " << t_frac50[i] << std::endl;
	      
	      TrailConstFractime(camp[i],&m_amprec[i][0],timeS,bck[i],0.3,ampl[i],NMax[i],&trail_t_frac30[i]);
	      TrailConstFractime(camp[i],&m_amprec[i][0],timeS,bck[i],0.5,ampl[i],NMax[i],&trail_t_frac50[i]);
	      TrailConstFractime(camp[i],&m_amprec[i][0],timeS,bck[i],0.7,ampl[i],NMax[i],&trail_t_frac70[i]);
	      TrailConstFractime(camp[i],&m_amprec[i][0],timeS,bck[i],0.9,ampl[i],NMax[i],&trail_t_frac90[i]);
	      //	      cout << int (t_frac30[i]/DT) << " " << int (t_frac10[i]/DT) << endl;
	
	      if (( (t_frac10[i]/DT) <camp[i] &&   (t_frac20[i]/DT) < camp[i] && (t_frac30[i]/DT)<camp[i]
		  &&  (t_frac70[i]/DT) <camp[i]&&   (t_frac80[i]/DT) <camp[i]) &&
		( (t_frac10[i]/DT) >0 &&   (t_frac20[i]/DT) >0 && (t_frac30[i]/DT)>0
		  &&  (t_frac70[i]/DT) >0 &&   (t_frac80[i]/DT) >0))
		{
	      
		  if ((t_frac30[i]- t_frac10[i]) != 0)
		    dVdt1030[i] =( m_amprec[i][ int (t_frac30[i]/DT)]- m_amprec[i][int (t_frac10[i]/DT)])/(t_frac30[i]- t_frac10[i]);
	
		  if  ((t_frac70[i]- t_frac30[i]) != 0)
		    dVdt3070[i] =( m_amprec[i][ int (t_frac70[i]/DT)]- m_amprec[i][int (t_frac30[i]/DT)])/(t_frac70[i]- t_frac30[i]);
	
		  if  ((t_frac80[i]- t_frac20[i]) != 0)
		    dVdt2080[i] =( m_amprec[i][ int (t_frac80[i]/DT)]- m_amprec[i][int (t_frac20[i]/DT)])/(t_frac80[i]- t_frac20[i]);
		}
	      // cout << dVdt3070[i] << endl;
	      

	      //    riselinfit(camp[i],&m_amprec[i][0],timeS, bck[i], t_level20[i],t_level40[i], &rise_lin0[i], &rise_lin1[i],&rise_lin_chi2[i]) ;

	      riselinfit(camp[i],&m_amprec[i][0],timeS, bck[i], t_frac30[i],t_frac50[i], &rise_lin0[i], &rise_lin1[i],&rise_lin_chi2[i]) ;
	      if (rise_lin1[i] != 0) t_zero[i]= (-rise_lin0[i])/rise_lin1[i];
	      riseexpfit(camp[i],&m_amprec[i][0],timeS, bck[i], t_frac20[i]-5,t_frac20[i], &rise_exp0[i], &rise_exp1[i],&rise_exp_chi2[i]) ;
	      
	      nprec=0;
	      //
	      maxD[i] = 0;
	      for(int l=0;l<camp[i];l++)
		{
		  m_amprec[i][l] = m_amprec[i][l]-bck[i];
		  if(l%CAMPFACT==0)
		    {
		      d_amprec[i][nprec]=der_amp[i][l];
		      if (fabs(d_amprec[i][nprec])>fabs(maxD[i])) 		       
			maxD[i] = fabs(d_amprec[i][nprec]);
			  
		      nprec++;
		    }
		}
	      

	    } // end loop on channels

	  
	  if(running==0 && emptycount<100)
	    {
	      running=1;
	      continue;
	    }
	  if(running==0 && emptycount>=100)
	    {
	      break;
	    }

	  
	  //cout << "check 3" << endl;
	  if (ampl[0] > 10 || ampl[1] > 10) 
	    {
	      w->Fill();
	      OutTree->Fill();
	    }
	  // cout << "check 4" << endl;
	}

 Write:

      continue;
    

    }

   OutputFile->cd();
   w->Write();
   OutTree->Write();
   OutputFile->Close();

   return(EXIT_SUCCESS);
}


