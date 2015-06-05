#include "TTPCDataQuality.hxx"
#include <TEventLoopFunction.hxx>
#include <TChannelInfo.hxx>

#include <TDigitContainer.hxx>
#include <TPulseDigit.hxx>
#include <TTPCChannelId.hxx>
#include <TChannelId.hxx>
#include <CaptGeomId.hxx>
#include <TGeometryInfo.hxx>

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
CP::TTPCDataQuality::TTPCDataQuality() {}

/////////////////////////////////////////////////////////////////
CP::TTPCDataQuality::~TTPCDataQuality() {}

/////////////////////////////////////////////////////////////////
// Called for each event inside the event loop, and returns true if the
// event should be saved to the output file.  If the remainder of the
// current file should be skipped, this should through the
// ENextEventLoopFile exception.
bool CP::TTPCDataQuality::operator () (CP::TEvent& event) {

    CP::TChannelInfo::Get().SetContext(event.GetContext());
    int eventid = event.GetEventId();
    int counter = 0;
    fTotalEvents++;
  
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
	    }

	    TH1F *waveform_proj_ch = new TH1F(Form("waveform_proj_ch%d_event%d",counter,eventid), Form("Event %d, Waveform projection for channel %d;ADC;Entries",eventid,counter),5000,0,5000);
            for (std::size_t i = 1; i< pulse->GetSampleCount(); ++i) {
	        waveform_proj_ch->Fill(pulse->GetSample(i));
            }

            fBaselineMean=waveform_proj_ch->GetMean();
            fPedestalRMS=waveform_proj_ch->GetRMS();

//Baseline check
          if(fApplyBaselineCheck) {
                fBaselineMean_h[counter]->Fill(fBaselineMean);
            }

//Pedestal RMS check
          if(fApplyPedestalRMSCheck) {
                fPedestalRMSMean_h[counter]->Fill(fPedestalRMS);
            }



//Pulse search check
            if(fApplyPulseSearchCheck) {
                for (std::size_t i = 1; i< pulse->GetSampleCount(); ++i) {
	            if(pulse->GetSample(i)-fBaselineMean > fPedestalRMS*7) {
		        fPulseCheckPositive_h[1]->Fill(counter);

		        if(fWirePlane.at(counter) != -1 && fWireNumber.at(counter) != -1) {
		            fPulseCheckPositive_h[0]->Fill((2 - fWirePlane.at(counter)) * 335 + fWireNumber.at(counter) + 1);
		            fPulseCheckPositive_h[2]->Fill(counter);
		        }
		        else {
		            fPulseCheckPositive_h[3]->Fill(counter);	  
		        }
	                break;
	            }
                }

                for (std::size_t i = 1; i< pulse->GetSampleCount(); ++i) {
		    if(pulse->GetSample(i)-fBaselineMean < -fPedestalRMS*7) {
		        fPulseCheckNegative_h[1]->Fill(counter);

		        if(fWirePlane.at(counter) != -1 && fWireNumber.at(counter) != -1) {
		            fPulseCheckNegative_h[0]->Fill((2 - fWirePlane.at(counter)) * 335 + fWireNumber.at(counter) + 1);
		            fPulseCheckNegative_h[2]->Fill(counter);
		        }
		        else {
		            fPulseCheckNegative_h[3]->Fill(counter);	  
		        }
	                break;
		    }
                }
	    }

//Pulse FFT check	  
	    if(fApplyFFTCheck) {
	    
	    }

//Pulse pulser height, rise time and pulse width check	  
	    if(fApplyPulserRunCheck) {
	    
	    }
	  
            delete waveform_proj_ch;
            counter++;

        }

//Event timing check
        if(fApplyEventTimingCheck) {
            if(eventid>1) fEventTiming_h->Fill(event.GetContext().GetTimeStamp() + event.GetContext().GetNanoseconds() * 1e-9 - fLastTimeStamp - fLastNanoSeconds * 1e-9);
            fLastTimeStamp = event.GetContext().GetTimeStamp();
            fLastNanoSeconds = event.GetContext().GetNanoseconds();
        }

    fLoadWireID=true;
    return false;
}

/////////////////////////////////////////////////////////////////
// Called after the arguments are processes by before reading the first
// event.  The output file is open so any histograms will be added to the
// output file.
void CP::TTPCDataQuality::Initialize(void) {
  
    fTotalChannelCount = 1152;
    fTotalWireCount = 1005;
    int frequency_up = 1000000;
    fTotalEvents = 0;
    fLoadWireID = false;
    const std::string name_suffix[4] = {"wire", "ch", "conn_ch", "disc_ch"};
    const std::string title_suffix[4] = {"Wire", "Channel", "Connected Channel", "Disconnected Channel"};
    const std::string titlex_suffix[4] = {"Wire [U (0,334), V (335,669), X(670,1004)]", "Channel", "Channel", "Channel"};
    const int totalbins[4] = {fTotalWireCount, fTotalChannelCount, fTotalChannelCount, fTotalChannelCount};

    if(fApplyBaselineCheck) {
        for (int i=0;i<4;i++) fBaseline_h[i] = new TH1F(("baseline_" + name_suffix[i]).c_str(), ("Baseline VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Baseline (ADC)").c_str(),totalbins[i],0,totalbins[i]);
        for (int ch=0;ch<fTotalChannelCount;ch++) fBaselineMean_h[ch] = new TH1F(Form("baselineMean%d",ch), Form("Baseline Mean%d;Baseline Mean (ADC);Entries",ch) , fBaseline_uplimit, 0, fBaseline_uplimit);
    }

    if(fApplyPedestalRMSCheck) {
        for (int i=0;i<4;i++) fPedestalRMS_h[i] = new TH1F(("pedestalRMS_" + name_suffix[i]).c_str(), ("Pedestal RMS VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Pedestal RMS (ADC)").c_str(),totalbins[i],0,totalbins[i]);
        for (int ch=0;ch<fTotalChannelCount;ch++) fPedestalRMSMean_h[ch] = new TH1F(Form("pedestalRMSMean%d",ch), Form("Pedestal RMS Mean%d;Pedestal RMS (ADC);Entries",ch) , fPedestalRMS_uplimit*10, 0, fPedestalRMS_uplimit);
    }
  
    if(fApplyShutOffCheck) {
        for (int i=0;i<4;i++) fShutoff_h[i] = new TH1F(("shutoff_" + name_suffix[i]).c_str(), ("Shutoff VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Shut off fraction (100%)").c_str(),totalbins[i],0,totalbins[i]);
    }
  
    if(fApplyPulseSearchCheck) {
        for (int i=0;i<4;i++) {
            fPulseCheckPositive_h[i] = new TH1F(("pulseCheckPositive_" + name_suffix[i]).c_str(), ("Pulse Positive VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Positive pulse (7#sigma) event fraction (100%)").c_str(),totalbins[i],0,totalbins[i]);
            fPulseCheckNegative_h[i] = new TH1F(("pulseCheckNegative_" + name_suffix[i]).c_str(), ("Pulse Negative VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Negative pulse (7#sigma) event fraction (100%)").c_str(),totalbins[i],0,totalbins[i]);
        }
    }
  
    if(fApplyFFTCheck) {
        for (int i=0;i<4;i++) fFFT_h[i] = new TH2F(("FFT_" + name_suffix[i]).c_str(), ("FFT VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Frequency (kHz)").c_str(),frequency_up,0,frequency_up,totalbins[i],0,totalbins[i]);
    }
  
    if(fApplyPulserRunCheck) {
        for (int i=0;i<4;i++) {
            fPulserHeightPositive_h[i] = new TH1F(("pulserHeightPositive_"+name_suffix[i]).c_str(), ("Pulser Height Positive VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Positive pulser height (ADC)").c_str(),totalbins[i],0,totalbins[i]);
            fPulserHeightNegative_h[i] = new TH1F(("pulserHeightNegative_"+name_suffix[i]).c_str(), ("Pulser Height Negative VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Negative pulser height (ADC)").c_str(),totalbins[i],0,totalbins[i]);    
            fPulserRiseTimePositive_h[i] = new TH1F(("pulserRiseTimePositive_"+name_suffix[i]).c_str(), ("Pulser Rise Time Positive VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Positive pulser height (ADC)").c_str(),totalbins[i],0,totalbins[i]);    
            fPulserRiseTimeNegative_h[i] = new TH1F(("pulserRiseTimeNegative_"+name_suffix[i]).c_str(), ("Pulser Rise Time Negative VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Negative pulser height (ADC)").c_str(),totalbins[i],0,totalbins[i]);    
            fPulserWidthPositive_h[i] = new TH1F(("pulserWidthPositive_"+name_suffix[i]).c_str(), ("Pulser Width Positive VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Positive pulser height (ADC)").c_str(),totalbins[i],0,totalbins[i]);    
            fPulserWidthNegative_h[i] = new TH1F(("pulserWidthNegative_"+name_suffix[i]).c_str(), ("Pulser Width Negative VS " + title_suffix[i] + ";" + titlex_suffix[i] + ";Negative pulser height (ADC)").c_str(),totalbins[i],0,totalbins[i]);
        }
    }

    if(fApplyEventTimingCheck) {
        fEventTiming_h = new TH1F("eventTiming", "Event Timing;Event Time Interval - 1/PulserFrequency (ms);Entries",int (4*fEventTimeBinning/fPulserFrequency),-2/fPulserFrequency,2/fPulserFrequency);
    }
}

/////////////////////////////////////////////////////////////////
// Called after reading the last event.  The output file is still open, so
// you can add extra information.  Because of an idiosyncrasy in the way
// root handles histograms, objects created in Initialize() will already be
// stored in the output file.
void CP::TTPCDataQuality::Finalize(CP::TRootOutput* const output) {
  
/// baseline check

    if(fApplyBaselineCheck) {
        for (int ch=0;ch<fTotalChannelCount;ch++)//starts at 1
        {
            fBaseline_h[1]->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean());
            fBaseline_h[1]->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
            if(ch<fWireNumber.size()) {
                if(fWirePlane.at(ch)!=-1 && fWireNumber.at(ch)!=-1) {
                    fBaseline_h[0]->SetBinContent((2-fWirePlane.at(ch))*335+fWireNumber.at(ch)+1,fBaselineMean_h[ch]->GetMean());
                    fBaseline_h[0]->SetBinError((2-fWirePlane.at(ch))*335+fWireNumber.at(ch)+1,fBaselineMean_h[ch]->GetRMS());
                    fBaseline_h[2]->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean());
                    fBaseline_h[2]->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
                }
                else {
                    fBaseline_h[3]->SetBinContent(ch+1,fBaselineMean_h[ch]->GetMean()); 
                    fBaseline_h[3]->SetBinError(ch+1,fBaselineMean_h[ch]->GetRMS());
                }
            }
        }
    }
 
/// pedestal RMS check
    if(fApplyPedestalRMSCheck) {
  
        for (int ch=0;ch<fTotalChannelCount;ch++) {
            fPedestalRMS_h[1]->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean());
            fPedestalRMS_h[1]->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
            if(ch<fWireNumber.size()) {
                if(fWirePlane.at(ch)!=-1 && fWireNumber.at(ch)!=-1) {
                    fPedestalRMS_h[0]->SetBinContent((2-fWirePlane.at(ch))*335+fWireNumber.at(ch)+1,fPedestalRMSMean_h[ch]->GetMean());
                    fPedestalRMS_h[0]->SetBinError((2-fWirePlane.at(ch))*335+fWireNumber.at(ch)+1,fPedestalRMSMean_h[ch]->GetRMS());
                    fPedestalRMS_h[2]->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean());
                    fPedestalRMS_h[2]->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
                }
                else {
                    fPedestalRMS_h[3]->SetBinContent(ch+1,fPedestalRMSMean_h[ch]->GetMean()); 
                    fPedestalRMS_h[3]->SetBinError(ch+1,fPedestalRMSMean_h[ch]->GetRMS());
                }
            }
        }
    }

/// shut off effect check
    if(fApplyShutOffCheck) {
        for (int ch=0;ch<fTotalChannelCount;ch++) {
            TF1 *f1;
            double par[3]={fPedestalRMSMean_h[ch]->GetMaximum(),fPedestalRMSMean_h[ch]->GetMaximumBin()/10.,fPedestalRMSMean_h[ch]->GetRMS()};
            if(par[2]<0.1) {
                f1 = new TF1("f1","gaus",par[1]-3*0.1,par[1]+3*0.1);
                f1->SetParameters(par);
                f1->SetParLimits(1,par[1]-3*0.1,par[1]+3*0.1);
                f1->SetParLimits(2,0,0.1);
            }
            else {
                f1 = new TF1("f1","gaus",par[1]-3*par[2],par[1]+3*par[2]);
                f1->SetParameters(par);
                f1->SetParLimits(1,par[1]-3*par[2],par[1]+3*par[2]);
                f1->SetParLimits(2,0,par[2]);
            }

            fPedestalRMSMean_h[ch]->Fit(f1,"R");
            int n_shutoff=0;
            double sevensigmarange=10*7*f1->GetParameter(2);
            sevensigmarange=std::max(sevensigmarange,5.);

            for (int bin=1;bin<10*f1->GetParameter(1)-sevensigmarange;bin++) n_shutoff+=fPedestalRMSMean_h[ch]->GetBinContent(bin);
	    
            if(n_shutoff/fPedestalRMSMean_h[ch]->GetEntries()<1) {
                fShutoff_h[1]->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
                if(ch<fWireNumber.size()) {
                    if(fWirePlane.at(ch)!=-1 && fWireNumber.at(ch)!=-1) {
                        fShutoff_h[0]->SetBinContent((2-fWirePlane.at(ch))*335+fWireNumber.at(ch)+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
                        fShutoff_h[2]->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
                    }
                    else {
                        fShutoff_h[3]->SetBinContent(ch+1,n_shutoff/fPedestalRMSMean_h[ch]->GetEntries());
                    }
                }
            }
        }
    }

/// pulse search check
    if(fApplyPedestalRMSCheck) {
        for (int i=0;i<4;i++) {
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
    if(!fApplyEventTimingCheck) {
        delete fEventTiming_h;
    }
    
/// write out all the relevant histograms on demand
    
    if(!fWriteBaselineChannels) {
        for (int ch=0;ch<fTotalChannelCount;ch++) {
            delete fBaselineMean_h[ch];
        }
    }

    if(!fWritePedestalRMSChannels) {
        for (int ch=0;ch<fTotalChannelCount;ch++) {
            delete fPedestalRMSMean_h[ch];
        }
    }

}


