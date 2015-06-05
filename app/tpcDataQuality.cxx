#include "TTPCDataQuality.hxx"
#include <eventLoop.hxx>

int main(int argc, char **argv) {
    CP::TTPCDataQuality userCode;
    userCode.ApplyBaselineCheck(true);
    userCode.ApplyPedestalRMSCheck(true);
    userCode.ApplyShutOffCheck(true);
    userCode.ApplyPulseSearchCheck(true);
    userCode.ApplyFFTCheck(false);
    userCode.ApplyPulserRunCheck(false);
    userCode.ApplyEventTimingCheck(true);
    userCode.WriteBaselineChannels(false); 
    userCode.WritePedestalRMSChannels(false);
    CP::eventLoop(argc,argv,userCode);
}
