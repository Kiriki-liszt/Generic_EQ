//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once
#include "GNRC_EQ_shared.h"
#include "GNRC_EQ_fft.h"

#include "public.sdk/source/vst/vstaudioeffect.h"

#include <vector>
#include <deque>     // std::deque
#include <numeric>   // std::transform_reduce
#include <algorithm> // std::fill, std::for_each

namespace yg331 {
//------------------------------------------------------------------------
//  GNRC_EQ_Processor
//------------------------------------------------------------------------
class GNRC_EQ_Processor : public Steinberg::Vst::AudioEffect
{
public:
    GNRC_EQ_Processor ();
    ~GNRC_EQ_Processor () SMTG_OVERRIDE;
    
    // Create function
    static Steinberg::FUnknown* createInstance (void* /*context*/)
    {
        return (Steinberg::Vst::IAudioProcessor*)new GNRC_EQ_Processor;
    }
    
    //------------------------------------------------------------------------
    // AudioEffect overrides:
    //------------------------------------------------------------------------
    /** Called at first after constructor */
    Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
    
    /** Called at the end before destructor */
    Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;
    
    /** Switch the Plug-in on/off */
    Steinberg::tresult PLUGIN_API setActive (Steinberg::TBool state) SMTG_OVERRIDE;
    
    /** Will be called before any process call */
    Steinberg::tresult PLUGIN_API setupProcessing (Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;
    
    /** Asks if a given sample size is supported see SymbolicSampleSizes. */
    Steinberg::tresult PLUGIN_API canProcessSampleSize (Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;
    
    /** Gets the current Latency in samples. */
    Steinberg::uint32 PLUGIN_API getLatencySamples() SMTG_OVERRIDE;
    
    /** Here we go...the process call */
    Steinberg::tresult PLUGIN_API process (Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;
    
    /** For persistence */
    Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    
    //------------------------------------------------------------------------
    // IConnectionPoint overrides:
    //------------------------------------------------------------------------
    /** Called when a message has been sent from the connection point to this. */
    // Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;
    
    //------------------------------------------------------------------------
protected:
    template <typename SampleType>
    void processSVF ( SampleType** inputs, SampleType** outputs, int32 numChannels, SampleRate getSampleRate, int32 sampleFrames );
    
    void call_after_SR_changed ();
    void call_after_parameter_changed ();
    
    uint16_t numChannels {0};
    
    bool       bBypass = false;
    ParamValue fLevel  = nrmParamLevl;
    bool       bPhase  = false;
    // ParamValue fZoom   = 2.0 / 6.0; // UNUSED
    
    // store in Norm Value
    std::array<std::array<ParamValue, bandSize>, numBands> pBand = {{
        {dftParamUsed, nrmParamType, nrmBand01Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand02Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand03Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand04Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand05Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand06Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand07Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand08Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand09Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand10Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand11Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand12Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand13Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand14Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand15Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand16Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand17Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand18Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand19Freq, nrmParamGain, nrmParamQlty},
        {dftParamUsed, nrmParamType, nrmBand20Freq, nrmParamGain, nrmParamQlty}
    }};
    
    std::array<std::array<ParamValue, bandSize>, numXover> pXovr = {{
        {dftParamUsed, nrmParamPass, nrmBand01Freq, nrmParamXtyp, nrmParamOrdr},
        {dftParamUsed, nrmParamPass, nrmBand02Freq, nrmParamXtyp, nrmParamOrdr}
    }};

    std::vector<std::array<SVF_Generic, numBands>> svf;      // vector size = numChannels
    std::vector<std::array<SVF_xover,   numXover>> svfXover; // vector size = numChannels

    // plugin enviroment
    SampleRate projectSR = 48000.0;
    SampleRate targetSR  = projectSR * OS_plain[OS_1x];
    int32 currLatency = latency_fir[OS_1x];
    
    // Oversampling and Latency
    int32      fParamOS = OS_1x; // Internal
    ParamValue fTarget = OS_1x;  // External Parameter
    std::vector<std::deque<ParamValue>> latencyDelayLine;    // vector size = numChannels
    std::array<ParamValue, Kaiser::maxTap> OS_coef;
    std::vector<std::array<ParamValue, Kaiser::maxTap>> OS_buff; // vector size = numChannels
};

//------------------------------------------------------------------------
} // namespace yg331
