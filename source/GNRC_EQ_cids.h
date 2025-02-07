//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace yg331 {
//------------------------------------------------------------------------
static const Steinberg::FUID kGNRC_EQ_ProcessorUID (0xBE9ACDC3, 0xEA6D5D86, 0xA6A76A31, 0xFAEED6B4);
static const Steinberg::FUID kGNRC_EQ_ControllerUID (0xE7DEC59C, 0x391E5034, 0x8FBCFD6A, 0x175F9AD6);

#define GNRC_EQ_VST3Category "Fx|EQ"

enum {
	kParamBypass = 0,
	kParamZoom,

	kParamLevel,
    kParamPhase,
    kParamTarget,
    kParamStereo,

	kParamBand01_Used, kParamBand01_Type, kParamBand01_Freq, kParamBand01_Gain, kParamBand01_Qlty,
    kParamBand02_Used, kParamBand02_Type, kParamBand02_Freq, kParamBand02_Gain, kParamBand02_Qlty,
    kParamBand03_Used, kParamBand03_Type, kParamBand03_Freq, kParamBand03_Gain, kParamBand03_Qlty,
    kParamBand04_Used, kParamBand04_Type, kParamBand04_Freq, kParamBand04_Gain, kParamBand04_Qlty,
    kParamBand05_Used, kParamBand05_Type, kParamBand05_Freq, kParamBand05_Gain, kParamBand05_Qlty,
    kParamBand06_Used, kParamBand06_Type, kParamBand06_Freq, kParamBand06_Gain, kParamBand06_Qlty,
    kParamBand07_Used, kParamBand07_Type, kParamBand07_Freq, kParamBand07_Gain, kParamBand07_Qlty,
    kParamBand08_Used, kParamBand08_Type, kParamBand08_Freq, kParamBand08_Gain, kParamBand08_Qlty,
    kParamBand09_Used, kParamBand09_Type, kParamBand09_Freq, kParamBand09_Gain, kParamBand09_Qlty,
    kParamBand10_Used, kParamBand10_Type, kParamBand10_Freq, kParamBand10_Gain, kParamBand10_Qlty,
    kParamBand11_Used, kParamBand11_Type, kParamBand11_Freq, kParamBand11_Gain, kParamBand11_Qlty,
    kParamBand12_Used, kParamBand12_Type, kParamBand12_Freq, kParamBand12_Gain, kParamBand12_Qlty,
    kParamBand13_Used, kParamBand13_Type, kParamBand13_Freq, kParamBand13_Gain, kParamBand13_Qlty,
    kParamBand14_Used, kParamBand14_Type, kParamBand14_Freq, kParamBand14_Gain, kParamBand14_Qlty,
    kParamBand15_Used, kParamBand15_Type, kParamBand15_Freq, kParamBand15_Gain, kParamBand15_Qlty,
    kParamBand16_Used, kParamBand16_Type, kParamBand16_Freq, kParamBand16_Gain, kParamBand16_Qlty,
    kParamBand17_Used, kParamBand17_Type, kParamBand17_Freq, kParamBand17_Gain, kParamBand17_Qlty,
    kParamBand18_Used, kParamBand18_Type, kParamBand18_Freq, kParamBand18_Gain, kParamBand18_Qlty,
    kParamBand19_Used, kParamBand19_Type, kParamBand19_Freq, kParamBand19_Gain, kParamBand19_Qlty,
    kParamBand20_Used, kParamBand20_Type, kParamBand20_Freq, kParamBand20_Gain, kParamBand20_Qlty,

    kParamBandX1_Used, kParamBandX1_Pass, kParamBandX1_Freq, kParamBandX1_Xtyp, kParamBandX1_Ordr,
    kParamBandX2_Used, kParamBandX2_Pass, kParamBandX2_Freq, kParamBandX2_Xtyp, kParamBandX2_Ordr
};
//------------------------------------------------------------------------
} // namespace yg331
