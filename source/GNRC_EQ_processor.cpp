//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "GNRC_EQ_processor.h"
#include "GNRC_EQ_cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"

#include "public.sdk/source/vst/vstaudioprocessoralgo.h"
#include "public.sdk/source/vst/vsthelpers.h"

using namespace Steinberg;

namespace yg331 {
//------------------------------------------------------------------------
// GNRC_EQ_Processor
//------------------------------------------------------------------------
GNRC_EQ_Processor::GNRC_EQ_Processor ()
{
    //--- set the wanted controller for our processor
    setControllerClass (kGNRC_EQ_ControllerUID);
}

//------------------------------------------------------------------------
GNRC_EQ_Processor::~GNRC_EQ_Processor ()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::initialize (FUnknown* context)
{
    // Here the Plug-in will be instantiated

    //---always initialize the parent-------
    tresult result = AudioEffect::initialize (context);
    // if everything Ok, continue
    if (result != kResultOk)
    {
        return result;
    }

    //--- create Audio IO ------
    addAudioInput  (STR16 ("Stereo In"),  Steinberg::Vst::SpeakerArr::kStereo);
    addAudioOutput (STR16 ("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

    /* If you don't need an event bus, you can remove the next line */
    //addEventInput (STR16 ("Event In"), 1);
    
    for (auto& it : OS_buff)
        std::fill(it.begin(), it.end(), 0.0);
    std::fill(OS_coef.begin(), OS_coef.end(), 0.0);

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!
    // fprintf (stdout, "GNRC_EQ_Processor::terminate\n");
    //---do not forget to call parent ------
    return AudioEffect::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::setActive (TBool state)
{
    call_after_SR_changed ();
    // fprintf (stdout, "GNRC_EQ_Processor::setActive\n");
    //--- called when the Plug-in is enable/disable (On/Off) -----
    return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::process (Vst::ProcessData& data)
{
    // fprintf (stdout, "GNRC_EQ_Processor::process\n");
    Vst::IParameterChanges* paramChanges = data.inputParameterChanges;

    if (paramChanges)
    {
        int32 numParamsChanged = paramChanges->getParameterCount();

        for (int32 index = 0; index < numParamsChanged; index++)
        {
            Vst::IParamValueQueue* paramQueue = paramChanges->getParameterData(index);

            if (paramQueue)
            {
                Vst::ParamValue value;
                int32 sampleOffset;
                int32 numPoints = paramQueue->getPointCount();

                if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
                {
                    int key = paramQueue->getParameterId();
                    int index  = std::clamp( (key - kParamBand01_Used) / bandSize,             0, numBands - 1); // 0 - 19
                    int xIndex = std::clamp(((key - kParamBand01_Used) / bandSize) - numBands, 0, numXover - 1); // 0 - 1
                    switch (key) {
                        case kParamBypass:  bBypass = (value > 0.5f); break;
                        // case kParamZoom:   fZoom = value; break;
                        case kParamLevel:   fLevel = value; break;
                        case kParamPhase:   bPhase = (value > 0.5f); break;
                        case kParamTarget:
                            fTarget = value;
                            call_after_SR_changed ();
                            sendTextMessage("latency_changed");
                            break;
                            
                        case kParamBand01_Used: case kParamBand02_Used: case kParamBand03_Used: case kParamBand04_Used: case kParamBand05_Used:
                        case kParamBand06_Used: case kParamBand07_Used: case kParamBand08_Used: case kParamBand09_Used: case kParamBand10_Used:
                        case kParamBand11_Used: case kParamBand12_Used: case kParamBand13_Used: case kParamBand14_Used: case kParamBand15_Used:
                        case kParamBand16_Used: case kParamBand17_Used: case kParamBand18_Used: case kParamBand19_Used: case kParamBand20_Used:
                            pBand[index][bandUsed]  = value; break;
                            
                        case kParamBand01_Type: case kParamBand02_Type: case kParamBand03_Type: case kParamBand04_Type: case kParamBand05_Type:
                        case kParamBand06_Type: case kParamBand07_Type: case kParamBand08_Type: case kParamBand09_Type: case kParamBand10_Type:
                        case kParamBand11_Type: case kParamBand12_Type: case kParamBand13_Type: case kParamBand14_Type: case kParamBand15_Type:
                        case kParamBand16_Type: case kParamBand17_Type: case kParamBand18_Type: case kParamBand19_Type: case kParamBand20_Type:
                            pBand[index][bandType]  = value; break;
                            
                        case kParamBand01_Freq: case kParamBand02_Freq: case kParamBand03_Freq: case kParamBand04_Freq: case kParamBand05_Freq:
                        case kParamBand06_Freq: case kParamBand07_Freq: case kParamBand08_Freq: case kParamBand09_Freq: case kParamBand10_Freq:
                        case kParamBand11_Freq: case kParamBand12_Freq: case kParamBand13_Freq: case kParamBand14_Freq: case kParamBand15_Freq:
                        case kParamBand16_Freq: case kParamBand17_Freq: case kParamBand18_Freq: case kParamBand19_Freq: case kParamBand20_Freq:
                            pBand[index][bandFreq]  = value; break;
                            
                        case kParamBand01_Gain: case kParamBand02_Gain: case kParamBand03_Gain: case kParamBand04_Gain: case kParamBand05_Gain:
                        case kParamBand06_Gain: case kParamBand07_Gain: case kParamBand08_Gain: case kParamBand09_Gain: case kParamBand10_Gain:
                        case kParamBand11_Gain: case kParamBand12_Gain: case kParamBand13_Gain: case kParamBand14_Gain: case kParamBand15_Gain:
                        case kParamBand16_Gain: case kParamBand17_Gain: case kParamBand18_Gain: case kParamBand19_Gain: case kParamBand20_Gain:
                            pBand[index][bandGain]  = value; break;
                            
                        case kParamBand01_Qlty: case kParamBand02_Qlty: case kParamBand03_Qlty: case kParamBand04_Qlty: case kParamBand05_Qlty:
                        case kParamBand06_Qlty: case kParamBand07_Qlty: case kParamBand08_Qlty: case kParamBand09_Qlty: case kParamBand10_Qlty:
                        case kParamBand11_Qlty: case kParamBand12_Qlty: case kParamBand13_Qlty: case kParamBand14_Qlty: case kParamBand15_Qlty:
                        case kParamBand16_Qlty: case kParamBand17_Qlty: case kParamBand18_Qlty: case kParamBand19_Qlty: case kParamBand20_Qlty:
                            pBand[index][bandQlty]  = value; break;
                            
                        case kParamBandX1_Used: case kParamBandX2_Used:
                            pXovr[xIndex][bandUsed] = value; break;
                            
                        case kParamBandX1_Pass: case kParamBandX2_Pass:
                            pXovr[xIndex][bandPass] = value; break;
                            
                        case kParamBandX1_Freq: case kParamBandX2_Freq:
                            pXovr[xIndex][bandFreq] = value; break;
                            
                        case kParamBandX1_Xtyp: case kParamBandX2_Xtyp:
                            pXovr[xIndex][bandXtyp] = value; break;
                            
                        case kParamBandX1_Ordr: case kParamBandX2_Ordr:
                            pXovr[xIndex][bandOrdr] = value; break;

                        default: break;
                    }
                }
            }
        }
        call_after_parameter_changed ();
    } // end if (paramChanges)

    if (data.numSamples <= 0)
        return kResultOk; // nothing to do
    
    if (data.numInputs == 0 || data.numOutputs == 0)
        return kResultOk; // nothing to do

    // (simplification) we suppose in this example that we have the same input channel count than
    // the output
    int32 numChannels = data.inputs[0].numChannels;

    //---get audio buffers----------------
    uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
    void** in  = getChannelBuffersPointer(processSetup, data.inputs[0]);
    void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
    Vst::SampleRate getSampleRate = processSetup.sampleRate;

    //---check if silence---------------
    // check if all channel are silent then process silent
    if (data.inputs[0].silenceFlags == Vst::getChannelMask(data.inputs[0].numChannels))
    {
        // mark output silence too (it will help the host to propagate the silence)
        data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

        // the plug-in has to be sure that if it sets the flags silence that the output buffer are
        // clear
        for (int32 i = 0; i < numChannels; i++)
        {
            // do not need to be cleared if the buffers are the same (in this case input buffer are
            // already cleared by the host)
            if (in[i] != out[i])
            {
                memset(out[i], 0, sampleFramesSize);
            }
        }
    }
    else {

        data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;
    
        if (data.symbolicSampleSize == Vst::kSample32) {
            processSVF<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, numChannels, getSampleRate, data.numSamples);
        }
        else {
            processSVF<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, numChannels, getSampleRate, data.numSamples);
        }
    }
    return kResultOk;
}

//------------------------------------------------------------------------
uint32 PLUGIN_API GNRC_EQ_Processor::getLatencySamples()
{
    // fprintf (stdout, "getLatencySamples = %d\n", currLatency);

    return currLatency;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
    //--- called before any processing ----
    // fprintf (stdout, "setupProcessing\n");
    // setupProcessing is not called in restartComponent.
    
    /*
     createInstance
     initialize
     Processor::getState
     setupProcessing
     setState
     getLatencySamples
     */
    
    projectSR = newSetup.sampleRate;
    
    // setupProcessing should happen after setBusArr.
    Vst::SpeakerArrangement arr;
    getBusArrangement (Vst::BusDirections::kInput, 0, arr);
    numChannels = static_cast<uint16_t> (Vst::SpeakerArr::getChannelCount (arr));
    
    svf.resize(numChannels);
    svfXover.resize(numChannels);
    latencyDelayLine.resize(numChannels);
    OS_buff.resize(numChannels);
    
    call_after_SR_changed (); // includes call_after_parameter_changed ()
    
    return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
    // by default kSample32 is supported
    if (symbolicSampleSize == Vst::kSample32)
        return kResultTrue;

    // disable the following comment if your processing support kSample64
    if (symbolicSampleSize == Vst::kSample64)
        return kResultTrue;

    return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::setState (IBStream* state)
{
    // fprintf (stdout, "Processor::setState\n");
    
    // called when we load a preset, the model has to be reloaded
    IBStreamer streamer(state, kLittleEndian);
    
    // 1. Read Plain Values
    int32           savedBypass = 0;
    // Vst::ParamValue savedZoom   = 0.0; // UNUSED, left for compatibility
    Vst::ParamValue savedLevel  = 0.0;
    int32           savedPhase  = 0;
    int32           savedTarget = 0;
    
    bandParamSet savedBand[numBands];
    xovrParamSet savedXovr[numXover];

    if (streamer.readInt32 (savedBypass) == false) return kResultFalse;
    // if (streamer.readDouble(savedZoom  ) == false) return kResultFalse;
    if (streamer.readDouble(savedLevel ) == false) return kResultFalse;
    if (streamer.readInt32 (savedPhase ) == false) return kResultFalse;
    if (streamer.readInt32 (savedTarget) == false) return kResultFalse;
    
    for (int bands = 0; bands < numBands; bands++)
    {
        if (streamer.readInt32 (savedBand[bands].Used) == false) return kResultFalse;
        if (streamer.readInt32 (savedBand[bands].Type) == false) return kResultFalse;
        if (streamer.readDouble(savedBand[bands].Freq) == false) return kResultFalse;
        if (streamer.readDouble(savedBand[bands].Gain) == false) return kResultFalse;
        if (streamer.readDouble(savedBand[bands].Qlty) == false) return kResultFalse;
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        if (streamer.readInt32 (savedXovr[bands].Used) == false) return kResultFalse;
        if (streamer.readInt32 (savedXovr[bands].Pass) == false) return kResultFalse;
        if (streamer.readDouble(savedXovr[bands].Freq) == false) return kResultFalse;
        if (streamer.readInt32 (savedXovr[bands].Xtyp) == false) return kResultFalse;
        if (streamer.readInt32 (savedXovr[bands].Ordr) == false) return kResultFalse;
    }

    // 2. Save as Norm Values
    bBypass = savedBypass > 0;
    // fZoom   = savedZoom;
    fLevel  = paramGain.ToNormalized(savedLevel);
    bPhase  = savedPhase > 0;
    fTarget = paramTrgt.ToNormalized(savedTarget);
    
    for (int bands = 0; bands < numBands; bands++)
    {
        pBand[bands][bandUsed] = savedBand[bands].Used;
        pBand[bands][bandType] = paramType.ToNormalized(savedBand[bands].Type);
        pBand[bands][bandFreq] = paramFreq.ToNormalized(savedBand[bands].Freq);
        pBand[bands][bandGain] = paramGain.ToNormalized(savedBand[bands].Gain);
        pBand[bands][bandQlty] = paramQlty.ToNormalized(savedBand[bands].Qlty);
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        pXovr[bands][bandUsed] = savedXovr[bands].Used;
        pXovr[bands][bandPass] = paramPass.ToNormalized(savedXovr[bands].Pass);
        pXovr[bands][bandFreq] = paramFreq.ToNormalized(savedXovr[bands].Freq);
        pXovr[bands][bandXtyp] = paramXtyp.ToNormalized(savedXovr[bands].Xtyp);
        pXovr[bands][bandOrdr] = paramOrdr.ToNormalized(savedXovr[bands].Ordr);
    }
    
    call_after_parameter_changed ();

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Processor::getState (IBStream* state)
{
    // fprintf (stdout, "Processor::getState\n");
    
    // here we need to save the model
    IBStreamer streamer(state, kLittleEndian);
    
    // Save in Plain Values
    streamer.writeInt32(bBypass ? 1 : 0);
    // streamer.writeDouble(fZoom);   // UNUSED, left for compatibility
    streamer.writeDouble(paramGain.ToPlain(fLevel));
    streamer.writeInt32(bPhase ? 1 : 0);
    streamer.writeInt32(paramTrgt.ToPlainList(fTarget));
    
    for (int bands = 0; bands < numBands; bands++)
    {
        streamer.writeInt32 (                      pBand[bands][bandUsed] > 0.5 ? 1 : 0);
        streamer.writeInt32 (paramType.ToPlainList(pBand[bands][bandType]));
        streamer.writeDouble(paramFreq.ToPlain    (pBand[bands][bandFreq]));
        streamer.writeDouble(paramGain.ToPlain    (pBand[bands][bandGain]));
        streamer.writeDouble(paramQlty.ToPlain    (pBand[bands][bandQlty]));
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        streamer.writeInt32 (                      pXovr[bands][bandUsed] > 0.5 ? 1 : 0);
        streamer.writeInt32 (paramPass.ToPlainList(pXovr[bands][bandPass]));
        streamer.writeDouble(paramFreq.ToPlain    (pXovr[bands][bandFreq]));
        streamer.writeInt32 (paramXtyp.ToPlainList(pXovr[bands][bandXtyp]));
        streamer.writeInt32 (paramOrdr.ToPlainList(pXovr[bands][bandOrdr]));
    }
    
    return kResultOk;
}

template <typename SampleType>
void GNRC_EQ_Processor::processSVF
 (
    SampleType** inputs,
    SampleType** outputs,
    Steinberg::int32 numChannels,
    Steinberg::Vst::SampleRate getSampleRate,
    Steinberg::int32 sampleFrames
  )
{
    Vst::Sample64 level = DecibelConverter::ToGain(paramGain.ToPlain(fLevel));
    int32 oversampling = OS_plain[fParamOS]; 
    
    for (int32 channel = 0; channel < numChannels; channel++)
    {
        int32 samples = 0;
        while (samples < sampleFrames)
        {
            Vst::Sample64 inputSample = inputs[channel][samples];
            Vst::Sample64 drySample = inputSample;
            inputSample *= level;
            
            double up_x[OS_plain[OS_num]] = { 0.0, };
            double up_y[OS_plain[OS_num]] = { 0.0, };

            up_x[0] = inputSample;

            // Process
            for (int i = 0; i < oversampling; i++)
            {
                Vst::Sample64 overSampled = up_x[i];

                for (int bands = 0; bands < numBands; bands++)
                    overSampled = svf[channel][bands].computeSVF(overSampled);
                
                for (int bands = 0; bands < numXover; bands++)
                    overSampled = svfXover[channel][bands].computeXover(overSampled);

                up_y[i] = overSampled;
            }
            
            if (fParamOS != OS_1x) // So, if x4 -> just quad band filter -> so I need only one-step filtering
            {
                std::copy(OS_buff[channel].begin(), OS_buff[channel].begin() + fir_taps[fParamOS] - oversampling, OS_buff[channel].begin() + oversampling);
                std::reverse(up_y, up_y + oversampling);
                std::copy(up_y, up_y + oversampling, OS_buff[channel].begin());
                // transform_reduce works faster in double[], and slow in std::deque<double>
                // but if loop order channel->sample, cache miss happens, and std::deque<double> works faster
                // Well, it just depends case-by-case.
                inputSample = std::transform_reduce(OS_coef.begin(), OS_coef.begin() + fir_taps[fParamOS], OS_buff[channel].data() + oversampling - 1, 0.0);
            }
            else // if (fParamOS == overSample_1x)
            {
                inputSample = up_y[0];
            }

            // Latency compensate
            latencyDelayLine[channel].push_back(drySample);
            Vst::Sample64 delayed = *(latencyDelayLine[channel].end() - 1 - currLatency);
            latencyDelayLine[channel].pop_front();
            
            if (bPhase)
                inputSample = -inputSample;

            if (bBypass)
                inputSample = delayed;

            outputs[channel][samples] = (SampleType)inputSample;

            samples++;
        }
    }
    return;
}

void GNRC_EQ_Processor::call_after_SR_changed ()
{
    int table [4][4] = {
        {OS_1x, OS_1x, OS_1x, OS_1x}, // Target : OS_1x, x1 [no oversampling]
        {OS_2x, OS_1x, OS_1x, OS_1x}, // Target : OS_2x, x2 [ 96 /  88.2 kHz]
        {OS_4x, OS_2x, OS_1x, OS_1x}, // Target : OS_4x, x4 [192 / 176.4 kHz]
        {OS_8x, OS_4x, OS_2x, OS_1x}  // Target : OS_8x, x8 [384 / 352.8 kHz]
    };
    int sampleRateType = 3;
    if      (projectSR <=  48000.0) sampleRateType = 0;
    else if (projectSR <=  96000.0) sampleRateType = 1;
    else if (projectSR <= 192000.0) sampleRateType = 2;
    else                            sampleRateType = 3;
    
    fParamOS = table[paramTrgt.ToPlainList(fTarget)][sampleRateType];
    targetSR = OS_plain[fParamOS] * projectSR;
    currLatency = latency_fir[fParamOS];
    
    Kaiser::calcFilter(OS_plain[fParamOS] * 48000.0, 0.0, 24000.0, fir_taps[fParamOS], 80.0, OS_coef.data());
    std::for_each(OS_coef.begin(), OS_coef.end(), [this](double &n) { n *= OS_plain[fParamOS]; });
    
    for (auto& it : OS_buff)
        std::fill(it.begin(), it.end(), 0.0);
    
    for (auto& it : latencyDelayLine)
        it.resize(currLatency, 0.0);
    Steinberg::Vst::IMessage* message = allocateMessage();
    if (message)
    {
        message->setMessageID("GUI");
        message->getAttributes()->setFloat("targetSR", targetSR);
        sendMessage(message);
    }

    call_after_parameter_changed ();
};
void GNRC_EQ_Processor::call_after_parameter_changed ()
{
    for (auto& filter : svf)
    {
        for (int bands = 0; bands < numBands; bands++)
        {
            int    Used =                       pBand[bands][bandUsed] > 0.5 ? 1 : 0;
            int    Type = paramType.ToPlainList(pBand[bands][bandType]);
            double Freq = paramFreq.ToPlain    (pBand[bands][bandFreq]);
            double Gain = paramGain.ToPlain    (pBand[bands][bandGain]);
            double Qlty = paramQlty.ToPlain    (pBand[bands][bandQlty]);
            filter[bands].setSVF(Used, Type, Freq, Gain, Qlty, targetSR);
        }
    }
    for (auto& filter : svfXover)
    {
        for (int bands = 0; bands < numXover; bands++)
        {
            int    Used =                       pXovr[bands][bandUsed] > 0.5 ? 1 : 0;
            int    Pass = paramPass.ToPlainList(pXovr[bands][bandPass]);
            double Freq = paramFreq.ToPlain    (pXovr[bands][bandFreq]);
            int    Xtyp = paramXtyp.ToPlainList(pXovr[bands][bandXtyp]);
            int    Ordr = paramOrdr.ToPlainList(pXovr[bands][bandOrdr]);
            filter[bands].setXover(Used, Pass, Freq, Xtyp, Ordr, targetSR);
        }
    }
};

//------------------------------------------------------------------------
} // namespace yg331
