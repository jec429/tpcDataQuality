#include "TTPCDataQuality.hxx"
#include <TEventLoopFunction.hxx>
#include <TChannelInfo.hxx>

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
#include <TChannelInfo.hxx>
#include <TChannelId.hxx>
#include <CaptGeomId.hxx>
#include <TGeometryInfo.hxx>

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
// Called for each event inside the event loop, and returns true if the
// event should be saved to the output file.  If the remainder of the
// current file should be skipped, this should through the
// ENextEventLoopFile exception.
bool CP::TTPCDataQuality::operator () (CP::TEvent& event) {

        CaptLog("Event " << event.GetContext());
        std::cout << event.GetContext().GetTimeStamp()<<std::endl;
        std::cout << event.GetContext().GetNanoseconds()<<std::endl;
        std::cout << event.GetContext().GetSpill()<<std::endl;
        std::cout << "################" <<std::endl;
	fTotalEvents++;
	
//fTotalEvents=event.GetEventId();
int eventid=event.GetEventId();
//std::cout<<fTotalEvents<<std::endl;
int counter=0;
//int wireplane, wirenumber;

        CP::THandle<CP::TDigitContainer> drift = event.Get<CP::TDigitContainer>("~/digits/drift");

        if (!drift) {
            CaptError("No tpc digits");
        }


        for (CP::TDigitContainer::iterator d = drift->begin();d != drift->end(); ++d) {
            const CP::TPulseDigit* pulse = dynamic_cast<const CP::TPulseDigit*>(*d);
	    
	    if(!fLoadWireID) {
	    CP::TChannelId chan(pulse->GetChannelId());
	    fWirePlane.push_back(CP::GeomId::Captain::GetWirePlane(CP::TChannelInfo::Get().GetGeometry(chan)));
	    fWireNumber.push_back(CP::GeomId::Captain::GetWireNumber(CP::TChannelInfo::Get().GetGeometry(chan)));
	    fLoadWireID=true;
	    }
	    //TH1F *waveform_ch = new TH1F(Form("waveform_ch%d_event%d",counter,eventid), Form("Event %d, Waveform for channel %d;T (#mus);ADC",eventid,counter),pulse->GetSampleCount(),0,pulse->GetSampleCount()/2);
	    TH1F *waveform_proj_ch = new TH1F(Form("waveform_proj_ch%d_event%d",counter,eventid), Form("Event %d, Waveform projection for channel %d;ADC;Entries",eventid,counter),5000,0,5000);
	    std::cout<<CP::TChannelInfo::Get().GetGeometry(pulse->GetChannelId(),0)<<std::endl;
	    //std::cout<<CP::TGeometryId CP::TChannelInfo::GetGeometry(pulse->GetChannelId(), 0)<<std::endl;
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
		fPulseCheckPositive_h[1]->Fill(counter);
		if(fWirePlane.at(counter)!=-1 && fWireNumber.at(counter)!=-1) {
		  fPulseCheckPositive_h[0]->Fill((2-fWirePlane.at(counter))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(counter))+fWireNumber.at(counter)+1);
		  fPulseCheckPositive_h[2]->Fill(counter);
		}
		else {
		  fPulseCheckPositive_h[3]->Fill(counter);	  
		}
		break;
	      }
            }
            
            for (std::size_t i = 1; i<= pulse->GetSampleCount(); ++i) {
            //waveform_ch->SetBinContent(i,pulse->GetSample(i));
	      if(pulse->GetSample(i)-fBaselineMean < -fPedestalRMS*7) {
		fPulseCheckNegative_h[1]->Fill(counter);
		if(fWirePlane.at(counter)!=-1 && fWireNumber.at(counter)!=-1) {
		  fPulseCheckNegative_h[0]->Fill((2-fWirePlane.at(counter))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(counter))+fWireNumber.at(counter)+1);
		  fPulseCheckNegative_h[2]->Fill(counter);
		}
		else {
		  fPulseCheckNegative_h[3]->Fill(counter);	  
		}
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
  //int wireCountX = CP::TGeometryInfo::Get().GetWireCount(0);
  //int wireCountV = CP::TGeometryInfo::Get().GetWireCount(1);
  //int wireCountU = CP::TGeometryInfo::Get().GetWireCount(2);
  
  fTotalWireCount = CP::TGeometryInfo::Get().GetWireCount(0)+CP::TGeometryInfo::Get().GetWireCount(1)+CP::TGeometryInfo::Get().GetWireCount(2);
  fTotalChannelCount = 1152;
  int frequency_up = 1000000;
  fTotalEvents = 0;
  const std::string suffix[4] = {"wire", "ch", "conn_ch", "disc_ch"};
  const int totalbins[4] = {fTotalWireCount, fTotalChannelCount, fTotalChannelCount, fTotalChannelCount};
  
  if(fApplyBaselineCheck) {
    
    for (int i=0;i<4;i++) {
    fBaseline_h[i] = new TH1F(("baseline_"+suffix[i]).c_str(), ("Baseline VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    }
    
    for (int ch=0;ch<fTotalWireCount;ch++) {
    fBaselineMean_h[ch] = new TH1F(Form("baselineMean%d",ch), Form("Baseline Mean%d",ch) , fBaseline_uplimit, 0, fBaseline_uplimit);
    }
  }
  
  if(fApplyPedestalRMSCheck) {
    
    for (int i=0;i<4;i++) {
    fPedestalRMS_h[i] = new TH1F(("pedestalRMS_"+suffix[i]).c_str(), ("Pedestal RMS VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    }
    
    for (int ch=0;ch<fTotalWireCount;ch++) {
    fPedestalRMSMean_h[ch] = new TH1F(Form("pedestalRMSMean%d",ch), Form("Pedestal RMS Mean%d",ch) , fPedestalRMS_uplimit*10, 0, fPedestalRMS_uplimit);
    }
  }
  
  if(fApplyShutOffCheck) {
    
    for (int i=0;i<4;i++) {
    fShutoff_h[i] = new TH1F(("shutoff_"+suffix[i]).c_str(), ("Shutoff VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    }

  }
  
  if(fApplyPulseSearchCheck) {
    
    for (int i=0;i<4;i++) {
    fPulseCheckPositive_h[i] = new TH1F(("pulseCheckPositive_"+suffix[i]).c_str(), ("Pulse Positive VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    fPulseCheckNegative_h[i] = new TH1F(("pulseCheckNegative_"+suffix[i]).c_str(), ("Pulse Negative VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    }

  }
  
  if(fApplyFFTCheck) {
    
    for (int i=0;i<4;i++) {
    fFFT_h[i] = new TH2F(("FFT_"+suffix[i]).c_str(), ("FFT VS "+suffix[i]).c_str(),frequency_up,0,frequency_up,totalbins[i],0,totalbins[i]);
    }
  }
  
  if(fApplyPulserRunCheck) {
    
    for (int i=0;i<4;i++) {
    fPulserHeightPositive_h[i] = new TH1F(("pulserHeightPositive_"+suffix[i]).c_str(), ("Pulser Height Positive VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    fPulserHeightNegative_h[i] = new TH1F(("pulserHeightNegative_"+suffix[i]).c_str(), ("Pulser Height Negative VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);    
    fPulserRiseTimePositive_h[i] = new TH1F(("pulserRiseTimePositive_"+suffix[i]).c_str(), ("Pulser Rise Time Positive VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);    
    fPulserRiseTimeNegative_h[i] = new TH1F(("pulserRiseTimeNegative_"+suffix[i]).c_str(), ("Pulser Rise Time Negative VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);    
    fPulserWidthPositive_h[i] = new TH1F(("pulserWidthPositive_"+suffix[i]).c_str(), ("Pulser Width Positive VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);    
    fPulserWidthNegative_h[i] = new TH1F(("pulserWidthNegative_"+suffix[i]).c_str(), ("Pulser Width Negative VS "+suffix[i]).c_str(),totalbins[i],0,totalbins[i]);
    }
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
  for (int ch=0;ch<fTotalChannelCount;ch++)//starts at 1
  {
     fBaseline_h[1]->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean());
     fBaseline_h[1]->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
     
     if(fWirePlane.at(ch)!=-1 && fWireNumber.at(ch)!=-1) {
       fBaseline_h[0]->SetBinContent((2-fWirePlane.at(ch))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(ch))+fWireNumber.at(ch)+1,fBaselineMean_h[ch]->GetMean());
       fBaseline_h[0]->SetBinError((2-fWirePlane.at(ch))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(ch))+fWireNumber.at(ch)+1,fBaselineMean_h[ch]->GetRMS());
       fBaseline_h[2]->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean());
       fBaseline_h[2]->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
     }
     else {
       fBaseline_h[3]->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean()); 
       fBaseline_h[3]->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
     }
  }
}

  
/// pedestal RMS check
if(fApplyPedestalRMSCheck) {
  
  for (int ch=0;ch<fTotalWireCount;ch++)//starts at 1
  {
     fPedestalRMS_h[1]->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean());
     fPedestalRMS_h[1]->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
     
     if(fWirePlane.at(ch)!=-1 && fWireNumber.at(ch)!=-1) {
       fPedestalRMS_h[0]->SetBinContent((2-fWirePlane.at(ch))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(ch))+fWireNumber.at(ch)+1,fPedestalRMSMean_h[ch]->GetMean());
       fPedestalRMS_h[0]->SetBinError((2-fWirePlane.at(ch))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(ch))+fWireNumber.at(ch)+1,fPedestalRMSMean_h[ch]->GetRMS());
       fPedestalRMS_h[2]->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean());
       fPedestalRMS_h[2]->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
     }
     else {
       fPedestalRMS_h[3]->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean()); 
       fPedestalRMS_h[3]->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
     }
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
    if(n_shutoff/fPedestalRMSMean_h[ch]->GetEntries()<1)
    {
       fShutoff_h[1]->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
       if(fWirePlane.at(ch)!=-1 && fWireNumber.at(ch)!=-1) {
       fShutoff_h[0]->SetBinContent((2-fWirePlane.at(ch))*CP::TGeometryInfo::Get().GetWireCount(fWirePlane.at(ch))+fWireNumber.at(ch)+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
       fShutoff_h[2]->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
       }
       else {
       fShutoff_h[3]->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
       }
    }
  
  }
}

/// pulse search check
if(fApplyPedestalRMSCheck) {
  for (int i=0;i<4;i++)
  {
    fPulseCheckPositive_h[i]->Scale(1/fTotalEvents);
    fPulseCheckNegative_h[i]->Scale(1/fTotalEvents);
  }
}

/// FFT check

if(fApplyFFTCheck) {
  
}

/// pulser height, rise time and pulse width check

if(fApplyPulserRunCheck) {
  
}

/// event timing check


}


