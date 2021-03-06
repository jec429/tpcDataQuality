#include <TEventLoopFunction.hxx>
#include <eventLoop.hxx>
#include <cmath>
#include <iostream>
#include <exception>
#include <list>
#include <TH1F.h>
#include <TH2F.h>

namespace CP {
    class TTPCDataQuality;
}

class CP::TTPCDataQuality: public CP::TEventLoopFunction {
public:
    /// Initialize any class specific variables, but most of the work can be
    /// done in Initialize.  Don't create histograms here!
    TTPCDataQuality();
    virtual ~TTPCDataQuality();

    /// Called for each event inside the event loop, and returns true if the
    /// event should be saved to the output file.  If the remainder of the
    /// current file should be skipped, this should through the
    /// ENextEventLoopFile exception.
    bool operator () (CP::TEvent& event);

    /// Called after the arguments are processes by before reading the first
    /// event.  The output file is open so any histograms will be added to the
    /// output file.
    virtual void Initialize(void);

    /// Called after reading the last event. The output file is still open,
    /// so you can add extra information.  Because of an idiosyncrasy in the
    /// way root handles histograms, objects created in Initialize() will
    /// already be stored in the output file.
    virtual void Finalize(CP::TRootOutput* const output);

    /// Turn on and off the baseline check.
    void ApplyBaselineCheck(bool value = true)		{fApplyBaselineCheck=value;};
    
    /// Write out basline distribution for all the channels
    void WriteBaselineChannels(bool value = false)	{fWriteBaselineChannels=value;};
    
    /// Turn on and off the pedestal RMS check.
    void ApplyPedestalRMSCheck(bool value = true)	{fApplyPedestalRMSCheck=value;};
    
    /// Write out pedestal RMS distribution for all the channels
    void WritePedestalRMSChannels(bool value = false)	{fWritePedestalRMSChannels=value;};
    
    /// Turn on and off the shut off effect check.
    void ApplyShutOffCheck(bool value = true)		{fApplyShutOffCheck=value;};
    
    /// Turn on and off the pulse search check.
    void ApplyPulseSearchCheck(bool value = true)	{fApplyPulseSearchCheck=value;};

    /// Turn on and off the FFT check.
    void ApplyFFTCheck(bool value = true)		{fApplyFFTCheck=value;};

    /// Turn on and off the pulser calibration check.
    void ApplyPulserRunCheck(bool value = true)		{fApplyPulserRunCheck=value;};
    
    /// Turn on and off the pulser calibration check.
    void ApplyEventTimingCheck(bool value = true)		{fApplyEventTimingCheck=value;};
    
private:

    /// Total number of channels
    int fTotalChannelCount;
    
    /// Total number of wires
    int fTotalWireCount;
    
    /// Store the correspondence of channel to wire
    std::vector<int> fWirePlane;
    std::vector<int> fWireNumber;
    
    /// Load wire ID in the first event processed
    bool fLoadWireID;
    
    /// Uplimit of the basline distribution histograms for all the channels
    const static int fBaseline_uplimit    = 5000;///ADC
    
    /// Uplimit of the pedestal RMS distribution histograms for all the channels
    const static int fPedestalRMS_uplimit = 1000;///ADC
    
    /// Pulser frequency that is used for calculating the interval between events
    const static double fPulserFrequency  = 1;///Hz
    
    /// 1 ms bin size
    const static int fEventTimeBinning=1000;///bins/s
    
    
    /// Pedestal RMS smaller than 7 sigma of their mean value gets selected as shut off
    double sevensigmarange;
    
    /// Baseline value
    double fBaselineMean;
    
    /// Pedestal RMS
    double fPedestalRMS;
    
    /// Total number of events in this run
    double fTotalEvents;
    
    /// integer part of last event time
    double fLastTimeStamp;
    
    /// decimal part of last event time
    double fLastNanoSeconds;
    
    
    bool fApplyBaselineCheck;
    bool fApplyPedestalRMSCheck;
    bool fApplyShutOffCheck;
    bool fApplyPulseSearchCheck;
    bool fApplyFFTCheck;
    bool fApplyPulserRunCheck;
    bool fApplyEventTimingCheck;
    
    bool fWriteBaselineChannels;
    bool fWritePedestalRMSChannels;
    //TTree* tree;
    
    /// Output histograms 
    TH1F* fBaseline_h[4];
    TH1F* fBaselineMean_h[1152];
    
    //std::vector<TH1F*> fBaselineMean_h;
    TH1F* fPedestalRMS_h[4];
    TH1F* fPedestalRMSMean_h[1152];
    
    //std::vector<TH1F*> fPedestalRMSMean_h;
    TH1F* fShutoff_h[4];
    
    TH1F* fPulseCheckPositive_h[4];
    
    TH1F* fPulseCheckNegative_h[4];
    
    TH2F* fFFT_h[4];
    
    TH1F* fPulserHeightPositive_h[4];
    TH1F* fPulserHeightNegative_h[4];
    TH1F* fPulserRiseTimePositive_h[4];
    TH1F* fPulserRiseTimeNegative_h[4];
    TH1F* fPulserWidthPositive_h[4];
    TH1F* fPulserWidthNegative_h[4];

    
    TH1F* fEventTiming_h;
};
