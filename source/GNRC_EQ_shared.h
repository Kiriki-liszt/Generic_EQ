//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/vst/vsttypes.h"
// #include "pluginterfaces/base/futils.h"

#define _USE_MATH_DEFINES
#include <cmath>     // Decibel
#include <algorithm> // copy
#include <functional> // multipliy
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
        setSVF (svf.Used,
                svf.Type,
                svf.Freq,
                svf.Gain,
                svf.Qlty,
                svf.Fs);
        initSVF ();
    };
    
    void   setUsed(bool v) { Used = v; }
    bool   getUsed() const { return Used; }
    
    void   setFreq(double v) { Freq = v; }
    double getFreq() const { return Freq; }
    
    void   setGain(double v) { Gain = v; }
    double getGain() const { return Gain; }
    
    void   setQlty(double v) { Qlty = v; }
    double getQlty() const { return Qlty; }
    
    void   setType(int v) { Type = v; }
    int    getType() const { return Type; }

    void   setFs(double v) { Fs = v; }
    double getFs() const { return Fs; }

    inline void initSVF()
    {
        ic1eq = 0.0;
        ic2eq = 0.0;
    };


    inline void setSVF (int plainUsed, int plainType, double plainFreq, double plainGain, double plainQlty, double plainFs)
    {
        Used = plainUsed ? 1 : 0;
        Type = plainType;
        Freq = plainFreq;
        Gain = plainGain;
        Qlty = plainQlty;
        Fs   = plainFs;
        divFS = 1.0 / Fs;
        
        switch (Type)
        {
            case tLowPass_6     : Ordr = o6dBoct; break;
            case tHighPass_6    : Ordr = o6dBoct; break;
            case tLowShelf_6    : Ordr = o6dBoct; break; // fc at -3dB
            case tHighShelf_6   : Ordr = o6dBoct; break; // fc at -3dB
                
            case tHighPass      : Ordr = o12dBoct; break;
            case tLowPass       : Ordr = o12dBoct; break;
            case tNotch         : Ordr = o12dBoct; break;
            case tAllPass       : Ordr = o12dBoct; break;
            case tBell          : Ordr = o12dBoct; break;
            case tLowShelf      : Ordr = o12dBoct; break;
            case tHighShelf     : Ordr = o12dBoct; break;
            case tLowShelf_12   : Ordr = o12dBoct; break;// fc at -3dB
            case tHighShelf_12  : Ordr = o12dBoct; break;// fc at -3dB
            default             : Ordr = o12dBoct; break;
        }
        
        makeSVF();
    }

    static SMTG_CONSTEXPR double dB2A = M_LN10 / 40.0;
    inline void makeSVF()
    {
        if (Freq > Fs * 0.5) Freq = Fs * 0.5;
        w = M_PI * Freq * divFS;
        g = std::tan(w);
        k = 1.0 / Qlty;
        double A = std::exp(Gain * dB2A);
        
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
        if (!Used)
            return vin;
        
        if (Ordr == o12dBoct)
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

        if (Ordr == o6dBoct)
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
        if (!Used)
            return 1.0;

        double ddr, ddi;
        getTransferFunc(freq, &ddr, &ddi);
        
        return std::sqrt(ddr * ddr + ddi * ddi);
    }
    
    inline double phs_response(double freq)
    {
        if (!Used)
            return 1.0;
        
        double ddr, ddi;
        getTransferFunc(freq, &ddr, &ddi);

        return std::atan2(ddi, ddr);
    }

    int    Used = 0;
    double Freq = 1000.0;
    double Qlty = M_SQRT1_2;
    double Gain = 0.0;
    int    Type = tBell;
    int    Ordr = o12dBoct;
    double Fs   = 48000.0;
    double divFS = 1.0 / Fs;

    double w = Freq * M_PI / Fs;
    double g = std::tan(w);
    double k = 1.0 / Qlty;
    double gt0 = 1 / (1 + g * (g + k));
    double gk0 = (g + k) * gt0;

    double m0 = 1.0, m1 = 0.0, m2 = 1.0;
    double v0 = 0.0, v1 = 0.0, v2 = 0.0;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0;
    double ic1eq = 0.0;
    double ic2eq = 0.0;
};

// Crossover Filters

/*
 Table of Butterworth
 Cascading filters to get higher Butterworth filters
 https://www.earlevel.com/main/2016/09/29/cascading-filters/

 1st / 6dboct
 0.5 (one-pole)
 
 2nd / 12dboct
 0.70710678118654752440084436210484903928483593768847403658833986899536623923 = 1/(2cos(π/4))
 
 3rd / 18dboct
 0.5 (one-pole)
 1.0000000 = 1/(2cos(π/3))
 
 4th / 24dboct
 0.54119610014619698439972320536638942006107206337801544468129709565298897354 = 1/(2cos(1π/8)) = Sqrt(1 - 1/sqrt(2))
 1.30656296487637652785664317342718715358376118834926952754889836690808104146 = 1/(2cos(3π/8))
 
 5th / 30dboct
 0.5 (one-pole)
 0.61803398874989484820458683436563811772030917980576286213544862270526046281 = 1/(2cos(1π/5))
 1.61803398874989484820458683436563811772030917980576286213544862270526046281 = 1/(2cos(2π/5))
 
 6th / 36dboct
 0.51763809020504152469779767524809665669813780263986102762800641463011394949 = 1/(2cos(1π/12))
 0.70710678118654752440084436210484903928483593768847403658833986899536623923 = 1/(2cos(3π/12))
 1.93185165257813657349948639945779473526780967801680910080468615262084642795 = 1/(2cos(5π/12))
 
 7th / 42dboct
 0.5 (one-pole)
 0.55495813208737119142219487100641048106728886247091008937602596820515753594 = 1/(2cos(1π/7))
 0.80193773580483825247220463901489010233183832426371430010712484639886484085 = 1/(2cos(2π/7))
 2.24697960371746706105000976800847962126454946179280421073109887819370730491 = 1/(2cos(3π/7))
 
 8th / 48dboct
 0.50979557910415916894193980398784391368261849190893536972644293702585452219 = 1/(2cos(1π/16))
 0.60134488693504528054372182390922183506725588960067272212186411938004397295 = 1/(2cos(3π/16))
 0.89997622313641570463850954094188839553800172741284498523544711601538115976 = 1/(2cos(5π/16))
 2.56291544774150617879608629617769990901512335138734491212769629647065700005 = 1/(2cos(7π/16))
 */

/*
 if LR, 12dboct need 1 x 12dboct Q = 0.5                        (two  6dboct butterworth = one 12dboct w/ Q)
 if LR, 24dBoct need 2 x 12dboct Q = 0.7071, 0.7071             (two 12dboct butterworth = total of 4 6dboct stacked)
 if LR, 36dboct need 3 x 12dboct Q = 0.5,    1.0,   1.0         (two 18dboct butterworth = 2 x (6db + 12dB to make single 18dBoct butterworth)
 if LR, 48dboct need 4 x 12dboct Q = (0.54, 1.34), (0.54, 1.34) (two 24dboct Butterworth filters = two (2x12dboct w/ Q == 24dboct butterworth))
 */
/*
 Bessel Filters
 F for Low pass filters. For High pass, divide by it.
 https://gist.github.com/endolith/4982787
 * denotes onepole
 Q for N =
  2: 0.5773502691896258
  3: ------------------ 0.691046625825072
  4: 0.8055382818416649 0.5219345816689801
  5: ------------------ 0.9164773739482458 0.5635356208514573
  6: 1.0233139538267275 0.6111945468780016 0.5103178247487693
  7: ------------------ 1.1262575419830458 0.6608213892970732 0.5323556978995495
  8: 1.225669425408165  0.7108520744417189 0.5596091647957803 0.5059910693974681

 f multiplier to get -3 dB at fc, for N =
  2: 1.2720196495140683
  3: 1.3226757999104435* 1.4476171331469865
  4: 1.6033575162169689  1.4301715599939901
  5: 1.502316271447484*  1.755377776637094  1.556347122296921
  6: 1.904707612302767   1.6891682676204571 1.6039191287738017
  7: 1.6843681792740095* 2.04949090027047   1.8224174788591163 1.7163560448722
  8: 2.188726230527587   1.9531957590222577 1.8320926011985377 1.7784659117748955
 */

class SVF_xover{
public:
    // compiler should take care of rounding to nearest IEEE 754
    static constexpr double BUTTERWORTH_noQ = 0.5; // 6dBoct, unused
    static constexpr double BUTTERWORTH_2_1 = 0.707106781186547524400844362104849039;
    static constexpr double BUTTERWORTH_3_1 = 0.5; // 6dboct
    static constexpr double BUTTERWORTH_3_2 = 1.0;
    static constexpr double BUTTERWORTH_4_1 = 0.541196100146196984399723205366389420;
    static constexpr double BUTTERWORTH_4_2 = 1.30656296487637652785664317342718715;
    static constexpr double BUTTERWORTH_5_1 = 0.5; // 6dboct
    static constexpr double BUTTERWORTH_5_2 = 0.618033988749894848204586834365638118;
    static constexpr double BUTTERWORTH_5_3 = 1.61803398874989484820458683436563812;
    static constexpr double BUTTERWORTH_6_1 = 0.517638090205041524697797675248096657;
    static constexpr double BUTTERWORTH_6_2 = 0.707106781186547524400844362104849039;
    static constexpr double BUTTERWORTH_6_3 = 1.93185165257813657349948639945779474;
    static constexpr double BUTTERWORTH_7_1 = 0.5; // 6dboct
    static constexpr double BUTTERWORTH_7_2 = 0.554958132087371191422194871006410481;
    static constexpr double BUTTERWORTH_7_3 = 0.801937735804838252472204639014890102;
    static constexpr double BUTTERWORTH_7_4 = 2.24697960371746706105000976800847962;
    static constexpr double BUTTERWORTH_8_1 = 0.509795579104159168941939803987843914;
    static constexpr double BUTTERWORTH_8_2 = 0.601344886935045280543721823909221835;
    static constexpr double BUTTERWORTH_8_3 = 0.899976223136415704638509540941888396;
    static constexpr double BUTTERWORTH_8_4 = 2.562915447741506178796086296177699909;

    // These are precise enough to convert into IEEE 754 Double, format provides 15–17 significant decimal digits of precision.
    static constexpr double BESSEL_F_2_1 = 1.2720196495140683;
    static constexpr double BESSEL_F_3_1 = 1.3226757999104435; // 6dboct
    static constexpr double BESSEL_F_3_2 = 1.4476171331469865;
    static constexpr double BESSEL_F_4_1 = 1.6033575162169689;
    static constexpr double BESSEL_F_4_2 = 1.4301715599939901;
    static constexpr double BESSEL_F_5_1 = 1.502316271447484; // 6dboct
    static constexpr double BESSEL_F_5_2 = 1.755377776637094;
    static constexpr double BESSEL_F_5_3 = 1.556347122296921;
    static constexpr double BESSEL_F_6_1 = 1.904707612302767;
    static constexpr double BESSEL_F_6_2 = 1.6891682676204571;
    static constexpr double BESSEL_F_6_3 = 1.6039191287738017;
    static constexpr double BESSEL_F_7_1 = 1.6843681792740095; // 6dboct
    static constexpr double BESSEL_F_7_2 = 2.04949090027047;
    static constexpr double BESSEL_F_7_3 = 1.8224174788591163;
    static constexpr double BESSEL_F_7_4 = 1.7163560448722;
    static constexpr double BESSEL_F_8_1 = 2.188726230527587;
    static constexpr double BESSEL_F_8_2 = 1.9531957590222577;
    static constexpr double BESSEL_F_8_3 = 1.8320926011985377;
    static constexpr double BESSEL_F_8_4 = 1.7784659117748955;

    static constexpr double BESSEL_Q_2_1 = 0.5773502691896258;
    static constexpr double BESSEL_Q_3_1 = 0.5; // 6dboct
    static constexpr double BESSEL_Q_3_2 = 0.691046625825072;
    static constexpr double BESSEL_Q_4_1 = 0.8055382818416649; // Pro-Q 1.13921408831
    static constexpr double BESSEL_Q_4_2 = 0.5219345816689801; // 0.7381340428
    static constexpr double BESSEL_Q_5_1 = 0.5; // 6dboct
    static constexpr double BESSEL_Q_5_2 = 0.9164773739482458;
    static constexpr double BESSEL_Q_5_3 = 0.5635356208514573;
    static constexpr double BESSEL_Q_6_1 = 1.0233139538267275;
    static constexpr double BESSEL_Q_6_2 = 0.6111945468780016;
    static constexpr double BESSEL_Q_6_3 = 0.5103178247487693;
    static constexpr double BESSEL_Q_7_1 = 0.5; // 6dboct
    static constexpr double BESSEL_Q_7_2 = 1.1262575419830458;
    static constexpr double BESSEL_Q_7_3 = 0.6608213892970732;
    static constexpr double BESSEL_Q_7_4 = 0.5323556978995495;
    static constexpr double BESSEL_Q_8_1 = 1.225669425408165;
    static constexpr double BESSEL_Q_8_2 = 0.7108520744417189;
    static constexpr double BESSEL_Q_8_3 = 0.5596091647957803;
    static constexpr double BESSEL_Q_8_4 = 0.5059910693974681;
    
    
    static constexpr int numFlt = 4;
    
    static constexpr int pLow           = 0;
    static constexpr int pHigh          = 1;
    static constexpr int pNum           = 1;
    static constexpr int pSize          = 2;
    static constexpr Steinberg::Vst::String128 Pass_Types[pSize] = {
        STR16("Low Pass"),
        STR16("High Pass")
    };
    
    static constexpr int tButterworth   = 0; // 6, 12, 18, 24, 30, 36, 42, 48
    static constexpr int tBessel        = 1; //    12, 28, 24, 30, 36, 42, 48
    static constexpr int tLinkwitzRiley = 2; //    12      24      36      48
    static constexpr int tNum           = 2; // Plain <-> Norm
    static constexpr int tSize          = 3; // Loops, array, ...

    static constexpr Steinberg::Vst::String128 Filter_Types[tSize] = {
        STR16("Butterworth"),
        STR16("Bessel"),
        STR16("Linkwitz–Riley")
    };

    static constexpr int o6dBoct  = 0; // 1st order
    static constexpr int o12dBoct = 1; // 2nd order
    static constexpr int o18dBoct = 2; // 3
    static constexpr int o24dBoct = 3;
    static constexpr int o30dBoct = 4;
    static constexpr int o36dBoct = 5;
    static constexpr int o42dBoct = 6;
    static constexpr int o48dBoct = 7;
    static constexpr int oNum     = 7; // Plain <-> Norm
    static constexpr int oSize    = 8; // Loops, array, ...
    
    static constexpr Steinberg::Vst::String128 Filter_Order[oSize] = {
        STR16(" 6dB/oct"),
        STR16("12dB/oct"),
        STR16("18dB/oct"),
        STR16("24dB/oct"),
        STR16("30dB/oct"),
        STR16("36dB/oct"),
        STR16("42dB/oct"),
        STR16("48dB/oct")
    };
    
    static constexpr int usedFilterByOrder[oSize][numFlt] = {
        {1, 0, 0, 0}, // o6dBoct
        {1, 0, 0, 0}, // o12dBoct
        {1, 1, 0, 0}, // o18dBoct
        {1, 1, 0, 0}, // o24dBoct
        {1, 1, 1, 0}, // o30dBoct
        {1, 1, 1, 0}, // o36dBoct
        {1, 1, 1, 1}, // o42dBoct
        {1, 1, 1, 1}  // o48dBoct
    };
    
    static constexpr int lpByOrder[oSize][numFlt] = {
        {SVF_Generic::tLowPass_6, 0, 0, 0}, // o6dBoct
        {SVF_Generic::tLowPass  , 0, 0, 0}, // o12dBoct
        {SVF_Generic::tLowPass_6, SVF_Generic::tLowPass  , 0, 0}, // o18dBoct
        {SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , 0, 0}, // o24dBoct
        {SVF_Generic::tLowPass_6, SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , 0}, // o30dBoct
        {SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , 0}, // o36dBoct
        {SVF_Generic::tLowPass_6, SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , SVF_Generic::tLowPass  }, // o42dBoct
        {SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , SVF_Generic::tLowPass  , SVF_Generic::tLowPass  }  // o48dBoct
    };
    static constexpr int hpByOrder[oSize][numFlt] = {
        {SVF_Generic::tHighPass_6, 0, 0, 0}, // o6dBoct
        {SVF_Generic::tHighPass  , 0, 0, 0}, // o12dBoct
        {SVF_Generic::tHighPass_6, SVF_Generic::tHighPass  , 0, 0}, // o18dBoct
        {SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , 0, 0}, // o24dBoct
        {SVF_Generic::tHighPass_6, SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , 0}, // o30dBoct
        {SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , 0}, // o36dBoct
        {SVF_Generic::tHighPass_6, SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , SVF_Generic::tHighPass  }, // o42dBoct
        {SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , SVF_Generic::tHighPass  , SVF_Generic::tHighPass  }  // o48dBoct
    };
    
    static constexpr double qltyButterworthByOrder[oSize][numFlt] = {
        {BUTTERWORTH_noQ, BUTTERWORTH_noQ, BUTTERWORTH_noQ, BUTTERWORTH_noQ}, // o6dBoct
        {BUTTERWORTH_2_1, BUTTERWORTH_noQ, BUTTERWORTH_noQ, BUTTERWORTH_noQ}, // o12dBoct
        {BUTTERWORTH_3_1, BUTTERWORTH_3_2, BUTTERWORTH_noQ, BUTTERWORTH_noQ}, // o18dBoct
        {BUTTERWORTH_4_1, BUTTERWORTH_4_2, BUTTERWORTH_noQ, BUTTERWORTH_noQ}, // o24dBoct
        {BUTTERWORTH_5_1, BUTTERWORTH_5_2, BUTTERWORTH_5_3, BUTTERWORTH_noQ}, // o30dBoct
        {BUTTERWORTH_6_1, BUTTERWORTH_6_2, BUTTERWORTH_6_3, BUTTERWORTH_noQ}, // o36dBoct
        {BUTTERWORTH_7_1, BUTTERWORTH_7_2, BUTTERWORTH_7_3, BUTTERWORTH_7_4}, // o42dBoct
        {BUTTERWORTH_8_1, BUTTERWORTH_8_2, BUTTERWORTH_8_3, BUTTERWORTH_8_4}  // o48dBoct
    };
    
    static constexpr double freqLpBesselByOrder[oSize][numFlt] = {
        {1.0, 1.0, 1.0, 1.0}, // o6dBoct -> Butterworth
        {BESSEL_F_2_1, 1.0, 1.0, 1.0}, // o12dBoct
        {BESSEL_F_3_1, BESSEL_F_3_2, 1.0, 1.0}, // o18dBoct
        {BESSEL_F_4_1, BESSEL_F_4_2, 1.0, 1.0}, // o24dBoct
        {BESSEL_F_5_1, BESSEL_F_5_2, BESSEL_F_5_3, 1.0}, // o30dBoct
        {BESSEL_F_6_1, BESSEL_F_6_2, BESSEL_F_6_3, 1.0}, // o36dBoct
        {BESSEL_F_7_1, BESSEL_F_7_2, BESSEL_F_7_3, BESSEL_F_7_4}, // o42dBoct
        {BESSEL_F_8_1, BESSEL_F_8_2, BESSEL_F_8_3, BESSEL_F_8_4}  // o48dBoct
    };
    static constexpr double freqHpBesselByOrder[oSize][numFlt] = {
        {1.0, 1.0, 1.0, 1.0}, // o6dBoct -> Butterworth
        {1.0/BESSEL_F_2_1, 1.0, 1.0, 1.0}, // o12dBoct
        {1.0/BESSEL_F_3_1, 1.0/BESSEL_F_3_2, 1.0, 1.0}, // o18dBoct
        {1.0/BESSEL_F_4_1, 1.0/BESSEL_F_4_2, 1.0, 1.0}, // o24dBoct
        {1.0/BESSEL_F_5_1, 1.0/BESSEL_F_5_2, 1.0/BESSEL_F_5_3, 1.0}, // o30dBoct
        {1.0/BESSEL_F_6_1, 1.0/BESSEL_F_6_2, 1.0/BESSEL_F_6_3, 1.0}, // o36dBoct
        {1.0/BESSEL_F_7_1, 1.0/BESSEL_F_7_2, 1.0/BESSEL_F_7_3, 1.0/BESSEL_F_7_4}, // o42dBoct
        {1.0/BESSEL_F_8_1, 1.0/BESSEL_F_8_2, 1.0/BESSEL_F_8_3, 1.0/BESSEL_F_8_4}  // o48dBoct
    };
    static constexpr double qltyBesselByOrder[oSize][numFlt] = {
        {1.0, 1.0, 1.0, 1.0}, // o6dBoct -> Butterworth
        {BESSEL_Q_2_1, 1.0, 1.0, 1.0}, // o12dBoct
        {BESSEL_Q_3_1, BESSEL_Q_3_2, 1.0, 1.0}, // o18dBoct
        {BESSEL_Q_4_1, BESSEL_Q_4_2, 1.0, 1.0}, // o24dBoct
        {BESSEL_Q_5_1, BESSEL_Q_5_2, BESSEL_Q_5_3, 1.0}, // o30dBoct
        {BESSEL_Q_6_1, BESSEL_Q_6_2, BESSEL_Q_6_3, 1.0}, // o36dBoct
        {BESSEL_Q_7_1, BESSEL_Q_7_2, BESSEL_Q_7_3, BESSEL_Q_7_4}, // o42dBoct
        {BESSEL_Q_8_1, BESSEL_Q_8_2, BESSEL_Q_8_3, BESSEL_Q_8_4}  // o48dBoct
    };
    
    static constexpr double qltyLinkwitzReilyByOrder[oSize][numFlt] = {
        {BUTTERWORTH_noQ, 1.0, 1.0, 1.0}, // o6dBoct -> 12dboct
        {BUTTERWORTH_noQ, 1.0, 1.0, 1.0}, // o12dBoct                            // Phase flip needed
        {BUTTERWORTH_2_1, BUTTERWORTH_2_1, 0, 0}, // o18dBoct -> 24dboct
        {BUTTERWORTH_2_1, BUTTERWORTH_2_1, 0, 0}, // o24dBoct                    // no flipped
        {BUTTERWORTH_3_1, BUTTERWORTH_3_2, BUTTERWORTH_3_2, 0}, // o30dBoct -> 36dboct
        {BUTTERWORTH_3_1, BUTTERWORTH_3_2, BUTTERWORTH_3_2, 0}, // o36dBoct      // Phase flip needed
        {BUTTERWORTH_4_1, BUTTERWORTH_4_2, BUTTERWORTH_4_1, BUTTERWORTH_4_2}, // o42dBoct -> 48dboct
        {BUTTERWORTH_4_1, BUTTERWORTH_4_2, BUTTERWORTH_4_1, BUTTERWORTH_4_2}  // o48dBoct // no flipped
    };
    
    SVF_xover()
    {
        initXover();
    }
    
    // Initializer Method
    SVF_xover(const SVF_xover& xov)
    {
        for (int i = 0; i < numFlt; i++)
            flt[i].setSVF (xov.flt[i].getUsed(),
                           xov.flt[i].getType(),
                           xov.flt[i].getFreq(),
                           xov.flt[i].getGain(),
                           xov.flt[i].getQlty(),
                           xov.flt[i].getFs());
        initXover();
    }
    
    void   setUsed(bool v) { Used = v; }
    bool   getUsed() const { return Used; }
    
    void   setPass(int v) { Pass = v; }
    int    getPass() const { return Pass; }
    
    void   setFreq(double v) { Freq = v; }
    double getFreq() const { return Freq; }
    
    void   setXtpe(int v) { Type = v; }
    int    getXtpe() const { return Type; }
    
    void   setOrdr(int v) { Ordr = v; }
    int    getOrdr() const { return Ordr; }

    void   setFs(double v) { Fs = v; }
    double getFs() const { return Fs; }
    
    
    inline void initXover()
    {
        for (int i = 0; i < numFlt; i++)
            flt[i].initSVF();
    }

    inline void setXover (int _plainUsed, int plainPass, double _plainFreq, int plainXtyp, int plainOrdr, double plainFs)
    {
        Used = _plainUsed;
        Pass = plainPass;
        Freq = _plainFreq;
        Type = plainXtyp;
        Ordr = plainOrdr;
        Fs = plainFs;
        
        int plainUsed[numFlt] = {0, }; // from _plainUsed
        int plainType[numFlt] = {0, }; // from plainPass
        double plainFreq[numFlt] = {1000.0, }; // from _plainFreq
        double plainQlty[numFlt] = {1.0, }; // from plainOrdr
        
        plainPass = std::max(std::min(plainPass, pNum), 0);
        plainOrdr = std::max(std::min(plainOrdr, oNum), 0);
        plainXtyp = std::max(std::min(plainXtyp, tNum), 0);
        
        if (plainXtyp == tButterworth)
        {
            std::copy(&usedFilterByOrder[plainOrdr][0], &usedFilterByOrder[plainOrdr][numFlt], &plainUsed[0]);
            
            if (plainPass == pLow)
                std::copy(&lpByOrder[plainOrdr][0], &lpByOrder[plainOrdr][numFlt], &plainType[0]);
            else // if (plainPass == pHigh)
                std::copy(&hpByOrder[plainOrdr][0], &hpByOrder[plainOrdr][numFlt], &plainType[0]);
            
            std::fill(&plainFreq[0], &plainFreq[numFlt], _plainFreq);
            
            std::copy(&qltyButterworthByOrder[plainOrdr][0], &qltyButterworthByOrder[plainOrdr][numFlt], &plainQlty[0]);
        }
        else if (plainXtyp == tBessel)
        {            
            double freqArray[numFlt];
            std::fill(&freqArray[0], &freqArray[numFlt], _plainFreq);
            
            std::copy(&usedFilterByOrder[plainOrdr][0], &usedFilterByOrder[plainOrdr][numFlt], &plainUsed[0]);
            
            if (plainPass == pLow)
            {
                std::copy(&lpByOrder[plainOrdr][0], &lpByOrder[plainOrdr][numFlt], &plainType[0]);
                std::transform(&freqLpBesselByOrder[plainOrdr][0], &freqLpBesselByOrder[plainOrdr][numFlt], freqArray, &plainFreq[0], std::multiplies<double>());
            }
            else // if (plainPass == pHigh)
            {
                std::copy(&hpByOrder[plainOrdr][0], &hpByOrder[plainOrdr][numFlt], &plainType[0]);
                std::transform(&freqHpBesselByOrder[plainOrdr][0], &freqHpBesselByOrder[plainOrdr][numFlt], freqArray, &plainFreq[0], std::multiplies<double>());
            }
            
            std::copy(&qltyBesselByOrder[plainOrdr][0], &qltyBesselByOrder[plainOrdr][numFlt], &plainQlty[0]);
        }
        else // if (plainXtyp == tLinkwitz)
        {
            if (plainOrdr ==  o6dBoct) plainOrdr = o12dBoct;
            if (plainOrdr == o18dBoct) plainOrdr = o24dBoct;
            if (plainOrdr == o30dBoct) plainOrdr = o36dBoct;
            if (plainOrdr == o42dBoct) plainOrdr = o48dBoct;
            
            std::copy(&usedFilterByOrder[plainOrdr][0], &usedFilterByOrder[plainOrdr][numFlt], &plainUsed[0]);
            
            if (plainPass == pLow)
                std::copy(&lpByOrder[plainOrdr][0], &lpByOrder[plainOrdr][numFlt], &plainType[0]);
            else // if (plainPass == pHigh)
                std::copy(&hpByOrder[plainOrdr][0], &hpByOrder[plainOrdr][numFlt], &plainType[0]);
            
            std::fill(&plainFreq[0], &plainFreq[numFlt], _plainFreq);
            
            std::copy(&qltyLinkwitzReilyByOrder[plainOrdr][0], &qltyLinkwitzReilyByOrder[plainOrdr][numFlt], &plainQlty[0]);
        }
        
        if (!_plainUsed)
            std::fill(&plainUsed[0], &plainUsed[numFlt], 0);
        
        for (int i = 0; i < numFlt; i++)
            flt[i].setSVF(plainUsed[i], plainType[i], plainFreq[i], 0.0, plainQlty[i], plainFs);
        
        // setSVF contains makeSVF()
    }
    
    inline double computeXover (double vin)
    {
        double sample = vin;
        for (int i = 0; i < numFlt; i++)
            sample = flt[i].computeSVF(sample);
        if (Type == tLinkwitzRiley && Pass == pHigh)
            if (Ordr == o12dBoct || Ordr == o36dBoct)
                sample = -sample;
        return sample;
    }
    
    inline double mag_response(double freq)
    {
        double sample = 1.0;
        for (int i = 0; i < numFlt; i++)
            sample *= flt[i].mag_response(freq);
        return sample;
    }
    
    SVF_Generic flt[numFlt];
    int    Used = 0;
    int    Pass = pLow;
    double Freq = 1000.0;
    int    Type = tButterworth;
    int    Ordr = o12dBoct;
    double Fs   = 48000.0;
};

//------------------------------------------------------------------------
//  Min, Max, Default of Parameters
//------------------------------------------------------------------------
static SMTG_CONSTEXPR bool dftBypass          = false;

static SMTG_CONSTEXPR ParamValue minParamFreq = 20.0;
static SMTG_CONSTEXPR ParamValue maxParamFreq = 22000.0;

static SMTG_CONSTEXPR int32 numBands  = 20;
static SMTG_CONSTEXPR int32 numXover  = 2;

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
static SMTG_CONSTEXPR int32      dftParamPass = SVF_xover::pLow;
static SMTG_CONSTEXPR int32      dftParamXtyp = SVF_xover::tButterworth;
static SMTG_CONSTEXPR int32      dftParamOrdr = SVF_xover::o12dBoct;

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
static const ParameterConverter paramPass      (0,            0,            ParameterConverter::list, SVF_xover::pNum);
static const ParameterConverter paramXtyp      (0,            0,            ParameterConverter::list, SVF_xover::tNum);
static const ParameterConverter paramOrdr      (0,            0,            ParameterConverter::list, SVF_xover::oNum);

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
static const double nrmParamPass = paramPass.ToNormalizedList(dftParamPass);
static const double nrmParamXtyp = paramXtyp.ToNormalizedList(dftParamXtyp);
static const double nrmParamOrdr = paramOrdr.ToNormalizedList(dftParamOrdr);

static SMTG_CONSTEXPR int bandUsed = 0;
static SMTG_CONSTEXPR int bandType = 1;
static SMTG_CONSTEXPR int bandPass = 1; // Xover
static SMTG_CONSTEXPR int bandFreq = 2;
static SMTG_CONSTEXPR int bandGain = 3;
static SMTG_CONSTEXPR int bandXtyp = 3; // Xover
static SMTG_CONSTEXPR int bandQlty = 4;
static SMTG_CONSTEXPR int bandOrdr = 4; // Xover
static SMTG_CONSTEXPR int bandSize = 5;
typedef struct bandParamSet {
    int    Used = 0;
    int    Type = dftParamType;
    double Freq = dftBand01Freq;
    double Gain = dftParamGain;
    double Qlty = dftParamQlty;
} bandParamSet;
typedef struct xovrParamSet {
    int    Used = 0;
    int    Pass = dftParamPass;
    double Freq = dftBand01Freq;
    int    Xtyp = dftParamXtyp;
    int    Ordr = dftParamOrdr;
} xovrParamSet;

}
