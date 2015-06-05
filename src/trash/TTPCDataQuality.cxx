#include "TTPCDataQuality.hxx"
#include <TEventLoopFunction.hxx>
//#include "TG4TrajectoriesModule.hxx"
//#include "TG4VerticesModule.hxx"
//#include <TGeometryInfo.hxx>
#include <TDigitContainer.hxx>
#include <TPulseDigit.hxx>
#include <TTPCChannelId.hxx>

#include <TRootInput.hxx>
#include <TRootOutput.hxx>
#include <TManager.hxx>

#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

#include <cmath>
#include <math.h>
#include <iostream>
#include <exception>
#include <list>

/////////////////////////////////////////////////////////////////
// Initialize any class specific variables, but most of the work can be
// done in Initialize.  Don't create histograms here!
CP::TTPCDataQuality::TTPCDataQuality() {
  std::cout<<" Constructing TtpcDataQuality ..."<<std::endl;
}

CP::TTPCDataQuality::~TTPCDataQuality() {
  std::cout<<" Destructing TtpcDataQuality ..."<<std::endl;
}

/////////////////////////////////////////////////////////////////
// Print a usage message.  This is generally called when there is a command
// line input error.
/*void CP::TTPCDataQuality::Usage(void) { }

/////////////////////////////////////////////////////////////////
// Set an option and return true if it is valid.  This is called by the
// event loop command line argument processing code for each "-O
// [name]=[value]" option found on the command line.  If the command line
// has "-O [name]" without a value, then the value string will be equal to
// "".  This must return false if the option was not correctly processed.
// Any -O name=value options where name is not explicitly given here will
// be checked against all the modules, and the value sent to
// module::Configure() in the case of a match.
bool CP::TTPCDataQuality::SetOption(std::string option,std::string value) {
    if (option == "save") {
        fSaveOriginalFullEvent = true;
        fSaveOutputEventTree = true;
        return true;
    }
    return false;
}
*/
/////////////////////////////////////////////////////////////////
// Called for each event inside the event loop, and returns true if the
// event should be saved to the output file.  If the remainder of the
// current file should be skipped, this should through the
// ENextEventLoopFile exception.
bool CP::TTPCDataQuality::operator () (CP::TEvent& event) {
    // Need to make sure we can get geoemtry before going through analysis
    // modules.
        CaptLog("Event " << event.GetContext());
        std::cout << event.GetContext().GetTimeStamp()<<std::endl;
        std::cout << event.GetContext().GetNanoseconds()<<std::endl;
        std::cout << event.GetContext().GetSpill()<<std::endl;
        std::cout << "################" <<std::endl;

fTotalEvents=event.GetEventId();
int eventid=event.GetEventId();
std::cout<<fTotalEvents<<std::endl;
int counter=0;

	
//TH1F *h[1152];
        CP::THandle<CP::TDigitContainer> drift = event.Get<CP::TDigitContainer>("~/digits/drift");

        if (!drift) {
            CaptError("No tpc digits");
        }


        for (CP::TDigitContainer::iterator d = drift->begin();d != drift->end(); ++d) {
            const CP::TPulseDigit* pulse = dynamic_cast<const CP::TPulseDigit*>(*d);
	    //TH1F *waveform_ch = new TH1F(Form("waveform_ch%d_event%d",counter,eventid), Form("Event %d, Waveform for channel %d;T (#mus);ADC",eventid,counter),pulse->GetSampleCount(),0,pulse->GetSampleCount()/2);
	    TH1F *waveform_proj_ch = new TH1F(Form("waveform_proj_ch%d_event%d",counter,eventid), Form("Event %d, Waveform projection for channel %d;ADC;Entries",eventid,counter),5000,0,5000);
	    
            for (std::size_t i = 1; i<= pulse->GetSampleCount(); ++i) {
            //waveform_ch->SetBinContent(i,pulse->GetSample(i));
	    waveform_proj_ch->Fill(pulse->GetSample(i));
            }
fBaselineMean=waveform_proj_ch->GetMean();
fPedestalRMS=waveform_proj_ch->GetRMS();
fBaselineMean_h[counter]->Fill(fBaselineMean);
fPedestalRMSMean_h[counter]->Fill(fPedestalRMS);

          if(fApplyPulseSearchCheck) {
            for (std::size_t i = 1; i<= pulse->GetSampleCount(); ++i) {
            //waveform_ch->SetBinContent(i,pulse->GetSample(i));
	      if(pulse->GetSample(i)-fBaselineMean > fPedestalRMS*7) {
		fPulseCheckPositive_h->Fill(counter);
		break;
	      }
            }
            
            for (std::size_t i = 1; i<= pulse->GetSampleCount(); ++i) {
            //waveform_ch->SetBinContent(i,pulse->GetSample(i));
	      if(pulse->GetSample(i)-fBaselineMean < -fPedestalRMS*7) {
		fPulseCheckNegative_h->Fill(counter);
		break;
	      }
            }
	  }
	  
delete waveform_proj_ch;
//h[counter]->Write();
//delete h[counter];
        counter++;

        }
if(fApplyEventTimingCheck) {
  if(eventid>1)fEventTiming_h->Fill(event.GetContext().GetTimeStamp()+event.GetContext().GetNanoseconds()*1e-9-fLastTimeStamp-fLastNanoSeconds*1e-9);
  fLastTimeStamp=event.GetContext().GetTimeStamp();
  fLastNanoSeconds=event.GetContext().GetNanoseconds();
}
//tree->Fill();
        return false;
}

/////////////////////////////////////////////////////////////////
// Called after the arguments are processes by before reading the first
// event.  The output file is open so any histograms will be added to the
// output file.
void CP::TTPCDataQuality::Initialize(void) {
  
  std::cout<<" Initialsing TtpcDataQuality ..."<<std::endl;
  /*int wireCountX = CP::TGeometryInfo::Get().GetWireCount(0);
  int wireCountV = CP::TGeometryInfo::Get().GetWireCount(1);
  int wireCountU = CP::TGeometryInfo::Get().GetWireCount(2);
  int totalWireCount=wireCountX+wireCountV+wireCountU;*/
  //fTotalWireCount=1152;
  int frequency_up = 1000000;
  
  if(fApplyBaselineCheck) {
    fBaseline_h = new TH1F("baseline", "Baseline",fTotalWireCount,0,fTotalWireCount);
    //fBaselineMean_h[fTotalWireCount];
    for (int ch=0;ch<fTotalWireCount;ch++) {
    fBaselineMean_h[ch] = new TH1F(Form("baselineMean%d",ch), Form("Baseline Mean%d",ch) , fBaseline_uplimit, 0, fBaseline_uplimit);
    }
  }
  
  if(fApplyPedestalRMSCheck) {
    fPedestalRMS_h = new TH1F("pedestalRMS", "Pedestal RMS",fTotalWireCount,0,fTotalWireCount);
    //TH1F* fPedestalRMSMean_h[fTotalWireCount];
    for (int ch=0;ch<fTotalWireCount;ch++) {
    fPedestalRMSMean_h[ch] = new TH1F(Form("pedestalRMSMean%d",ch), Form("Pedestal RMS Mean%d",ch) , fPedestalRMS_uplimit*10, 0, fPedestalRMS_uplimit);
    }
  }
  
  if(fApplyShutOffCheck) {
    fShutoff_h = new TH1F("shutoff", "Shutoff",fTotalWireCount,0,fTotalWireCount);
  }
  
  if(fApplyPulseSearchCheck) {
    fPulseCheckPositive_h = new TH1F("pulseCheckPositive", "Pulse Check Positive",fTotalWireCount,0,fTotalWireCount);
    fPulseCheckNegative_h = new TH1F("pulseCheckNegative", "Pulse Check Negative",fTotalWireCount,0,fTotalWireCount);
  }
  
  if(fApplyFFTCheck) {
    fFFT_h = new TH2F("FFT", "FFT",frequency_up,0,frequency_up,fTotalWireCount,0,fTotalWireCount);
  }
  
  if(fApplyPulserRunCheck) {
    fPulserHeightPositive_h = new TH1F("pulserHeightPositive", "Pulser Height Positive",fTotalWireCount,0,fTotalWireCount);
    fPulserHeightNegative_h = new TH1F("pulserHeightNegative", "Pulser Height Negative",fTotalWireCount,0,fTotalWireCount);
    fPulserRiseTimePositive_h = new TH1F("pulserRiseTimePositive", "Pulser Rise Time Positive",fTotalWireCount,0,fTotalWireCount);
    fPulserRiseTimeNegative_h = new TH1F("pulserRiseTimeNegative", "Pulser Rise Time Negative",fTotalWireCount,0,fTotalWireCount);
    fPulserWidthPositive_h = new TH1F("pulserWidthPositive", "Pulser Width Positive",fTotalWireCount,0,fTotalWireCount);
    fPulserWidthNegative_h = new TH1F("pulserWidthNegative", "Pulser Width Negative",fTotalWireCount,0,fTotalWireCount);
  }

  if(fApplyEventTimingCheck) {
    fEventTiming_h = new TH1F("eventTiming", "Event Timing;Event Time Interval - 1/PulserFrequency (ms);Entries",int (4*fEventTimeBinning/fPulserFrequency),-2/fPulserFrequency,2/fPulserFrequency);
  }
}

/////////////////////////////////////////////////////////////////
// Called before the first event of a file is read, but you should prefer
// Initialize() for general initialization.  This method will be called
// once for each input file.
/*void CP::TTPCDataQuality::BeginFile(CP::TVInputFile *vinput) {
    CP::TRootInput *input = dynamic_cast<CP::TRootInput *>(vinput);
    if (!input) return;
    TFile* file = input->GetFilePointer();
    if (!file) return;
    for (ModuleList::iterator m  = fAnalysisModules.begin();
        m != fAnalysisModules.end(); ++m) {
        (*m)->SetBeginFile(file);
    }
}*/

/////////////////////////////////////////////////////////////////
// Called after the last event of a file is read, but you should
// prefer Finalize() for general finalization.  This method will
// be called once for each input file.
//void CP::TTPCDataQuality::EndFile(CP::TVInputFile *) {}

/////////////////////////////////////////////////////////////////
// Called after reading the last event.  The output file is still open, so
// you can add extra information.  Because of an idiosyncrasy in the way
// root handles histograms, objects created in Initialize() will already be
// stored in the output file.
void CP::TTPCDataQuality::Finalize(CP::TRootOutput* const output) {
std::cout<<" Finalising TtpcDataQuality ..."<<std::endl;
/// baseline check
if(fApplyBaselineCheck) {
  
  for (int ch=0;ch<fTotalWireCount;ch++)//starts at 1
  {
     fBaseline_h->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean());
     fBaseline_h->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
  }
}  
  
/// pedestal RMS check
if(fApplyPedestalRMSCheck) {
  
  for (int ch=0;ch<fTotalWireCount;ch++)//starts at 1
  {
     fPedestalRMS_h->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean());
     fPedestalRMS_h->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
  }
}

/// shut off effect check
if(fApplyShutOffCheck) {
  
  for (int ch=0;ch<fTotalWireCount;ch++)//starts at 1
  {
    TF1 *f1;
    double par[3]={fPedestalRMSMean_h[ch]->GetMaximum(),fPedestalRMSMean_h[ch]->GetMaximumBin()/10.,fPedestalRMSMean_h[ch]->GetRMS()};
  
    if(par[2]<0.1)
    {
      f1 = new TF1("f1","gaus",par[1]-3*0.1,par[1]+3*0.1);
      f1->SetParameters(par);
      f1->SetParLimits(1,par[1]-3*0.1,par[1]+3*0.1);
      f1->SetParLimits(2,0,0.1);
    }
    else
    {
      f1 = new TF1("f1","gaus",par[1]-3*par[2],par[1]+3*par[2]);
      f1->SetParameters(par);
      f1->SetParLimits(1,par[1]-3*par[2],par[1]+3*par[2]);
      f1->SetParLimits(2,0,par[2]);
    }
  
    //if(fPedestalRMSMean_h[ch]->GetEntries()>0)
      fPedestalRMSMean_h[ch]->Fit(f1,"R");
    int n_shutoff=0;
  
    double sevensigmarange=10*7*f1->GetParameter(2);
    sevensigmarange=std::max(sevensigmarange,5.);

    for (int bin=1;bin<10*f1->GetParameter(1)-sevensigmarange;bin++)
    {
      n_shutoff+=fPedestalRMSMean_h[ch]->GetBinContent(bin);
    }
    if(n_shutoff/fPedestalRMSMean_h[ch]->GetEntries()<1)fShutoff_h->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
  
  }
}

/// pulse search check
if(fApplyPedestalRMSCheck) {
  fPulseCheckPositive_h->Scale(1/fTotalEvents);
  fPulseCheckNegative_h->Scale(1/fTotalEvents);
}

/// FFT check

if(fApplyFFTCheck) {
  
}

/// pulser height, rise time and pulse width check

if(fApplyPulserRunCheck) {
  
}

/// event timing check


}


