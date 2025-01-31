//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/vst/vsttypes.h"
// #include "pluginterfaces/base/futils.h"

#define _USE_MATH_DEFINES
#include <cmath>     // Decibel
#include <vector>

namespace yg331 {
//------------------------------------------------------------------------
using ParamValue    = Steinberg::Vst::ParamValue;
using SampleRate    = Steinberg::Vst::SampleRate;
using int32         = Steinberg::int32;
using uint32        = Steinberg::uint32;
using TBool         = Steinberg::TBool;

//------------------------------------------------------------------------
//  Class for converter
//------------------------------------------------------------------------
class DecibelConverter
{
public:
    static SMTG_CONSTEXPR double log2dB = 20.0 / M_LN10;
    static SMTG_CONSTEXPR double dB2log = M_LN10 / 20.0;
    static inline ParamValue ToGain (ParamValue dB)
    {
        // return std::pow(10.0, dB * 0.05);
        return std::exp(dB * dB2log);
    }

    static inline ParamValue ToDecibel (ParamValue gain)
    {
        return (gain > 0.0) ? (log2dB * std::log(gain)) : (-140.0);
    }
private:
    DecibelConverter() = delete;
};

class ParameterConverter
{
public:
    enum paramType {
        range = 0,
        log,
        list
    };
    
    ParameterConverter (ParamValue minValue, ParamValue maxValue, paramType type = paramType::range, int32 numSteps = -1)
    : minValue (minValue), maxValue (maxValue), type(type), numSteps (numSteps) {};

    ParamValue ToNormalized (ParamValue plain) const
    {
        switch (type) {
            case paramType::range: return ToNormalizedRange(plain); break;
            case paramType::log:   return ToNormalizedLog(plain);   break;
            case paramType::list:  return ToNormalizedList(plain);  break;
            default: break;
        }
        return ToNormalizedRange(plain);
    }

    ParamValue ToPlain (ParamValue normalized) const
    {
        switch (type) {
            case paramType::range: return ToPlainRange(normalized); break;
            case paramType::log:   return ToPlainLog(normalized);   break;
            case paramType::list:  return ToPlainList(normalized);  break;
            default: break;
        }
        return ToNormalizedRange(normalized);
    }
    
    ParamValue ToNormalizedRange (ParamValue plain) const
    {
        if (maxValue == minValue) return 0.0;
        return (plain - minValue) / (maxValue - minValue);
    }

    ParamValue ToPlainRange (ParamValue normalized) const
    {
        return normalized * (maxValue - minValue) + minValue;
    }
    
    /* Log Conversion
     
      x - x0    log(y) - log(y0)
     ------- = -----------------
     x1 - x0   log(y1) - log(y0)
     
     x0 = 0, x1 = 1 : normalized parameter
     
     norm  = ( log(plain) - log(minValue) ) / ( log(maxValue) - log(minValue) )
           = log(plain / minValue) / log(maxValue / minValue)
     
     log(plain / minValue) = norm * log(maxValue / minValue)
     plain / minValue = exp(norm * log(maxValue / minValue))
     plain = minValue * exp(norm * log(maxValue / minValue))
     
     */
    
    ParamValue ToNormalizedLog (ParamValue plain) const
    {
        if (minValue == 0.0) return 0.0;
        if (!((plain    / minValue) > 0.0)) return 0.0;
        if (!((maxValue / minValue) > 0.0)) return 0.0;
        return std::log(plain / minValue) / std::log(maxValue / minValue);
    }

    ParamValue ToPlainLog (ParamValue normalized) const
    {
        if (minValue == 0.0) return 0.0;
        if (!((maxValue / minValue) > 0.0)) return 0.0;
        return minValue * std::exp(normalized * std::log(maxValue / minValue));
    }
    
    ParamValue ToNormalizedList (int32 plain) const
    {
        // return (Steinberg::ToNormalized<type> (v, step));
        if (numSteps <= 0) return 0;
        return plain / ParamValue (numSteps);
    }

    int32 ToPlainList (ParamValue normalized) const
    {
        // return (Steinberg::FromNormalized<type> (v, step));
        if (numSteps <= 0) return 0;
        int32 cmp(normalized * (numSteps + 1));
        return numSteps < cmp ? numSteps : cmp;
    }

    void setRange (ParamValue _min, ParamValue _max)
    {
        minValue = _min;
        maxValue = _max;
    }
    
    void setSteps (int32 _steps)
    {
        numSteps = _steps;
    }
    
private:
    ParamValue minValue = 0.0;
    ParamValue maxValue = 1.0;
    paramType type = paramType::range;
    int32 numSteps = -1;
};

/**
 Based on;
 Simultaneous solving of all outputs of Linear SVF using trapezoidal integration in state space form
 © Andrew Simper, Cytomic, 2021, andy@cytomic.com
 last updated: 24th Sep 2021
 https://cytomic.com/files/dsp/SvfLinearTrapAllOutputs.pdf
**/
/*
 Bell       - typical
 LP,   HP   - typical
 LP6,  HP6  - 1st order
 LS,   HS   - typical
 LS6,  HS6  - 1st order, fc at -3dB
 HS12, HS12 - 2nd order, fc at -3dB
 Notch      - typical
 AP         - typical
 */
class SVF_Generic {
public:
    static constexpr int tBell          = 0;
    static constexpr int tLowPass       = 1; // No Gain
    static constexpr int tHighPass      = 2; // No Gain
    static constexpr int tLowPass_6     = 3; // No Gain, No Q
    static constexpr int tHighPass_6    = 4; // No Gain, No Q
    static constexpr int tLowShelf      = 5;
    static constexpr int tHighShelf     = 6;
    static constexpr int tLowShelf_6    = 7; // No Gain, Q // fc at -3dB
    static constexpr int tHighShelf_6   = 8; // No Gain, Q // fc at -3dB
    static constexpr int tLowShelf_12   = 9; // No Gain, Q // fc at -3dB
    static constexpr int tHighShelf_12  = 10; // No Gain, Q // fc at -3dB
    static constexpr int tNotch         = 11; // No Gain
    static constexpr int tAllPass       = 12; // No Gain
    static constexpr int tNum           = 12;
    static constexpr int tSize          = 13;

    static constexpr Steinberg::Vst::String128 Filter_Types[tSize] = {
        STR16("Bell"),
        STR16("Low Pass"),
        STR16("High Pass"),
        STR16("Low Pass 6"),
        STR16("High Pass 6"),
        STR16("Low Shelf"),
        STR16("High Shelf"),
        STR16("L Shelf 6"),
        STR16("H Shelf 6"),
        STR16("L Shelf 12"),
        STR16("H Shelf 12"),
        STR16("Notch"),
        STR16("All Pass")
    };

    static constexpr int o6dBoct  = 0;
    static constexpr int o12dBoct = 1;
    static constexpr int oNum     = 1;

    SVF_Generic()
    {
        initSVF();
    };
    
    // Initializer Method
    SVF_Generic(const SVF_Generic& svf)
    {
        setSVF (svf.In,
                svf.Hz,
                svf.Q,
                svf.dB,
                svf.Type,
                svf.Fs);
        initSVF ();
    };
    
    void   setIn(bool v) { In = v; }
    bool   getIn() const { return In; }
    
    void   setFreq(double v) { Hz = v; }
    double getFreq() const { return Hz; }
    
    void   setGain(double v) { dB = v; }
    double getGain() const { return dB; }
    
    void   setQ(double v) { Q = v; }
    double getQ() const { return Q; }
    
    void   setType(int v) { Type = v; }
    int    getType() const { return Type; }

    void   setFs(double v) { Fs = v; }
    double getFs() const { return Fs; }

    inline void initSVF()
    {
        ic1eq = 0.0;
        ic2eq = 0.0;
    };


    inline void setSVF (int plainIn, int plainType, double plainHz, double plaindB, double plainQ, double plainFs)
    {
        In    = plainIn ? 1 : 0;
        Hz    = plainHz;
        dB    = plaindB;
        Q     = plainQ;
        Type  = plainType;
        Fs    = plainFs;
        divFS = 1.0 / Fs;
        
        switch (Type)
        {
            case tLowPass_6     : Order = o6dBoct; break;
            case tHighPass_6    : Order = o6dBoct; break;
            case tLowShelf_6    : Order = o6dBoct; break; // fc at -3dB
            case tHighShelf_6   : Order = o6dBoct; break; // fc at -3dB
                
            case tHighPass      : Order = o12dBoct; break;
            case tLowPass       : Order = o12dBoct; break;
            case tNotch         : Order = o12dBoct; break;
            case tAllPass       : Order = o12dBoct; break;
            case tBell          : Order = o12dBoct; break;
            case tLowShelf      : Order = o12dBoct; break;
            case tHighShelf     : Order = o12dBoct; break;
            case tLowShelf_12   : Order = o12dBoct; break;// fc at -3dB
            case tHighShelf_12  : Order = o12dBoct; break;// fc at -3dB
            default             : Order = o12dBoct; break;
        }
        
        makeSVF();
    }

    static SMTG_CONSTEXPR double dB2A = M_LN10 / 40.0;
    inline void makeSVF()
    {
        if (Hz > Fs * 0.5) Hz = Fs * 0.5;
        w = M_PI * Hz * divFS;
        g = std::tan(w);
        k = 1.0 / Q;
        double A = std::exp(dB * dB2A);
        
        switch (Type)
        {
            case tLowPass       : m0 = 0;     m1 = 0;     m2 = 1;   break;
            case tHighPass      : m0 = 1;     m1 = 0;     m2 = 0;   break;
            case tNotch         : m0 = 1;     m1 = 0;     m2 = 1;   break;
            case tAllPass       : m0 = 1;     m1 = -k;    m2 = 1;   break;
            case tBell          : m0 = 1;     m1 = k * A; m2 = 1;     k = k / A;       break;
            case tLowShelf      : m0 = 1;     m1 = k * A; m2 = A * A; g = g / sqrt(A); break;
            case tHighShelf     : m0 = A * A; m1 = k * A; m2 = 1;     g = g * sqrt(A); break;
                
            case tLowPass_6     : m0 = 0;     m1 = 0;     m2 = 1;   break;
            case tHighPass_6    : m0 = 1;     m1 = 0;     m2 = 0;   break;

            case tLowShelf_6    : m0 = 1;     m1 = 0;     m2 = A * A; break; // fc at -3dB
            case tHighShelf_6   : m0 = A * A; m1 = 0;     m2 = 1;     break; // fc at -3dB

            case tLowShelf_12   : m0 = 1;     m1 = A * M_SQRT2; m2 = A * A; k = M_SQRT2; break; // fc at -3dB
            case tHighShelf_12  : m0 = A * A; m1 = A * M_SQRT2; m2 = 1;     k = M_SQRT2; break; // fc at -3dB
                
            default: break;
        }

        gt0 = 1.0 / (1.0 + g * (g + k));
        gk0 = (g + k) * gt0;
        return;
    };

    inline double tick_6dBoct(double vin) {
        // disable v1 stage
        t0 = vin - ic2eq;
        v0 = t0 / (1.0 + g);// gt0 * t0;
        t2 = g * v0;
        v2 = ic2eq + t2;
        ic2eq += 2.0 * t2;

        return m0 * v0 + m2 * v2;
    }
    inline double tick_12dBoct(double vin) {
        // tick serial(possibly quicker on cpus with low latencies)
        t0 = vin - ic2eq;
        v0 = gt0 * t0 - gk0 * ic1eq; // high
        t1 = g * v0;
        v1 = ic1eq + t1; // band
        t2 = g * v1;
        v2 = ic2eq + t2; // low
        ic1eq += 2.0 * t1;
        ic2eq += 2.0 * t2;

        return m0 * v0 + m1 * v1 + m2 * v2;
    }

    inline double computeSVF (double vin)
    {
        if (In != 1)
            return vin;
        
        if (Order == o12dBoct)
            return tick_12dBoct(vin);
        else
            return tick_6dBoct(vin);
    };
    
    inline void getTransferFunc(double freq, double* _ddr, double* _ddi)
    {
        // exp(complex(0.0, -2.0 * pi) * frequency / sampleRate)
        double _zr = (0.0        ) * freq * divFS;
        double _zi = (-2.0 * M_PI) * freq * divFS;

        // z = zr + zi;
        double zr = std::exp(_zr) * std::cos(_zi);
        double zi = std::exp(_zr) * std::sin(_zi);

        double nr = 0, ni = 0;
        double dr = 0, di = 0;

        if (Order == o6dBoct)
        {
            // Numerator complex
            nr = zr * (-m0 /* + m1 * (g - 1) */ + m2 * g) + (m0 /* + m1 * (g + 1) */ + m2 * g);
            ni = zi * (-m0 /* + m1 * (g - 1) */ + m2 * g);

            // Denominator complex
            dr = zr * (g - 1) + (g + 1);
            di = zi * (g - 1);
        }
        else {
            // z * z
            double zsq_r = zr * zr - zi * zi;
            double zsq_i = zi * zr + zr * zi;
            double gsq = g * g;

            // Numerator complex
            double c_nzsq = (m0        + m1 * g    + m2 * gsq);
            double c_nz   = (m0 * -2.0             + m2 * 2.0 * gsq);
            double c_n    = (m0        + m1 * -g   + m2 * gsq);
            nr = zsq_r * c_nzsq + zr * c_nz + c_n;
            ni = zsq_i * c_nzsq + zi * c_nz;

            // Denominator complex
            double c_dzsq = ( 1.0 + k *  g +       gsq);
            double c_dz   = (-2.0 +          2.0 * gsq);
            double c_d    = ( 1.0 + k * -g +       gsq);
            dr = zsq_r * c_dzsq + zr * c_dz + c_d;
            di = zsq_i * c_dzsq + zi * c_dz;
        }

        // Numerator / Denominator
        double norm = dr * dr + di * di;
        double ddr = (nr * dr + ni * di) / norm;
        double ddi = (ni * dr - nr * di) / norm;
        
        *_ddr = ddr;
        *_ddi = ddi;
    }

    inline double mag_response(double freq)
    {
        if (!In) return 1.0;

        double ddr, ddi;
        getTransferFunc(freq, &ddr, &ddi);
        
        return std::sqrt(ddr * ddr + ddi * ddi);
    }
    
    inline double phs_response(double freq)
    {
        if (!In) return 1.0;
        
        double ddr, ddi;
        getTransferFunc(freq, &ddr, &ddi);

        return std::atan2(ddi, ddr);
    }
    
    
    

    int    In   = 0;
    double Hz   = 1000.0;
    double Q    = M_SQRT1_2;
    double dB   = 0.0;
    int    Type = tBell;
    int    Order = o12dBoct;
    double Fs   = 48000.0;
    double divFS = 1.0 / Fs;

    double w = Hz * M_PI / Fs;
    double g = std::tan(w);
    double k = 1.0 / Q;
    double gt0 = 1 / (1 + g * (g + k));
    double gk0 = (g + k) * gt0;

    double m0 = 1.0, m1 = 0.0, m2 = 1.0;
    double v0 = 0.0, v1 = 0.0, v2 = 0.0;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0;
    double ic1eq = 0.0;
    double ic2eq = 0.0;
};

//------------------------------------------------------------------------
//  Min, Max, Default of Parameters
//------------------------------------------------------------------------
static SMTG_CONSTEXPR bool dftBypass          = false;

static SMTG_CONSTEXPR ParamValue minParamFreq = 20.0;
static SMTG_CONSTEXPR ParamValue maxParamFreq = 22000.0;

static SMTG_CONSTEXPR int32 numBands  = 20;

static SMTG_CONSTEXPR ParamValue dftBand01Freq = 20.0;
static SMTG_CONSTEXPR ParamValue dftBand02Freq = 25.0;
static SMTG_CONSTEXPR ParamValue dftBand03Freq = 31.5;
static SMTG_CONSTEXPR ParamValue dftBand04Freq = 40.0;
static SMTG_CONSTEXPR ParamValue dftBand05Freq = 50.0;
static SMTG_CONSTEXPR ParamValue dftBand06Freq = 63.0;
static SMTG_CONSTEXPR ParamValue dftBand07Freq = 80.0;
static SMTG_CONSTEXPR ParamValue dftBand08Freq = 100.0;
static SMTG_CONSTEXPR ParamValue dftBand09Freq = 125.0;
static SMTG_CONSTEXPR ParamValue dftBand10Freq = 160.0;
static SMTG_CONSTEXPR ParamValue dftBand11Freq = 200.0;
static SMTG_CONSTEXPR ParamValue dftBand12Freq = 250.0;
static SMTG_CONSTEXPR ParamValue dftBand13Freq = 315.0;
static SMTG_CONSTEXPR ParamValue dftBand14Freq = 400.0;
static SMTG_CONSTEXPR ParamValue dftBand15Freq = 500.0;
static SMTG_CONSTEXPR ParamValue dftBand16Freq = 630.0;
static SMTG_CONSTEXPR ParamValue dftBand17Freq = 800.0;
static SMTG_CONSTEXPR ParamValue dftBand18Freq = 1000.0;
static SMTG_CONSTEXPR ParamValue dftBand19Freq = 1250.0;
static SMTG_CONSTEXPR ParamValue dftBand20Freq = 1600.0;
static SMTG_CONSTEXPR ParamValue dftBandFreq[numBands] = {
    dftBand01Freq, dftBand02Freq, dftBand03Freq, dftBand04Freq, dftBand05Freq,
    dftBand06Freq, dftBand07Freq, dftBand08Freq, dftBand09Freq, dftBand10Freq,
    dftBand11Freq, dftBand12Freq, dftBand13Freq, dftBand14Freq, dftBand15Freq,
    dftBand16Freq, dftBand17Freq, dftBand18Freq, dftBand19Freq, dftBand20Freq};

static SMTG_CONSTEXPR ParamValue minParamGain = -12.0;
static SMTG_CONSTEXPR ParamValue maxParamGain = 12.0;
static SMTG_CONSTEXPR ParamValue dftParamGain = 0.0;

static SMTG_CONSTEXPR ParamValue minParamQlty = 0.1;
static SMTG_CONSTEXPR ParamValue maxParamQlty = 50.0;
static SMTG_CONSTEXPR ParamValue dftParamQlty = M_SQRT1_2;

static SMTG_CONSTEXPR int32      dftParamUsed = 0.0;
static SMTG_CONSTEXPR int32      dftParamType = SVF_Generic::tBell;

// zoom title-value struct
struct ZoomFactor
{
    const Steinberg::tchar* title;
    double factor;

    ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
};
typedef std::vector<ZoomFactor> ZoomFactorVector;
static const ZoomFactorVector zoomFactors {ZoomFactor(STR("50%"),   0.50),
                                           ZoomFactor(STR("75%"),   0.75),
                                           ZoomFactor(STR("100%"),  1.00),
                                           ZoomFactor(STR("125%"),  1.25),
                                           ZoomFactor(STR("150%"),  1.50),
                                           ZoomFactor(STR("175%"),  1.75),
                                           ZoomFactor(STR("200%"),  2.00)};
static SMTG_CONSTEXPR int32 zoomNum = 6;
static SMTG_CONSTEXPR int32 dftZoom = 2;

static const ParameterConverter paramGain      (minParamGain, maxParamGain, ParameterConverter::range);
static const ParameterConverter paramFreq      (minParamFreq, maxParamFreq, ParameterConverter::log);
static const ParameterConverter paramQlty      (minParamQlty, maxParamQlty, ParameterConverter::log);
static const ParameterConverter paramType      (0,            0,            ParameterConverter::list, SVF_Generic::tNum);

static const double nrmBand01Freq = paramFreq.ToNormalized(dftBand01Freq);
static const double nrmBand02Freq = paramFreq.ToNormalized(dftBand02Freq);
static const double nrmBand03Freq = paramFreq.ToNormalized(dftBand03Freq);
static const double nrmBand04Freq = paramFreq.ToNormalized(dftBand04Freq);
static const double nrmBand05Freq = paramFreq.ToNormalized(dftBand05Freq);
static const double nrmBand06Freq = paramFreq.ToNormalized(dftBand06Freq);
static const double nrmBand07Freq = paramFreq.ToNormalized(dftBand07Freq);
static const double nrmBand08Freq = paramFreq.ToNormalized(dftBand08Freq);
static const double nrmBand09Freq = paramFreq.ToNormalized(dftBand09Freq);
static const double nrmBand10Freq = paramFreq.ToNormalized(dftBand10Freq);
static const double nrmBand11Freq = paramFreq.ToNormalized(dftBand11Freq);
static const double nrmBand12Freq = paramFreq.ToNormalized(dftBand12Freq);
static const double nrmBand13Freq = paramFreq.ToNormalized(dftBand13Freq);
static const double nrmBand14Freq = paramFreq.ToNormalized(dftBand14Freq);
static const double nrmBand15Freq = paramFreq.ToNormalized(dftBand15Freq);
static const double nrmBand16Freq = paramFreq.ToNormalized(dftBand16Freq);
static const double nrmBand17Freq = paramFreq.ToNormalized(dftBand17Freq);
static const double nrmBand18Freq = paramFreq.ToNormalized(dftBand18Freq);
static const double nrmBand19Freq = paramFreq.ToNormalized(dftBand19Freq);
static const double nrmBand20Freq = paramFreq.ToNormalized(dftBand20Freq);

static const double nrmParamGain = paramGain.ToNormalized(dftParamGain);
static const double nrmParamQlty = paramQlty.ToNormalized(dftParamQlty);
static const double nrmParamType = paramQlty.ToNormalizedList(dftParamType);

static SMTG_CONSTEXPR int bandUsed = 0;
static SMTG_CONSTEXPR int bandType = 1;
static SMTG_CONSTEXPR int bandFreq = 2;
static SMTG_CONSTEXPR int bandGain = 3;
static SMTG_CONSTEXPR int bandQlty = 4;
static SMTG_CONSTEXPR int bandSize = 5;
typedef struct bandParamSet {
    int    Used = 0;
    int    Type = dftParamType;
    double Freq = dftBand01Freq;
    double Gain = dftParamGain;
    double Qlty = dftParamQlty;
} bandParamSet;

}
