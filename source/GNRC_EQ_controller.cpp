//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "GNRC_EQ_controller.h"
#include "GNRC_EQ_cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "pluginterfaces/base/ustring.h"
#include "base/source/fstreamer.h"

//#include "vstgui/vstgui.h"
#include "vstgui/vstgui_uidescription.h"
#include "vstgui/uidescription/detail/uiviewcreatorattributes.h"

using namespace Steinberg;

namespace VSTGUI {

EQCurveView::EQCurveView(
    const CRect& size,
    IControlListener* listener,
    int32_t tag,
    CBitmap* background
)
    : CControl(size, listener, tag, background)
{
    BackColor = kWhiteCColor;
    LineColor = kBlackCColor;
    BorderColor = kBlackCColor;
    byPass = false;
    level = 0.0;
    idleRate = 60;
    setWantsIdle(true);
}
EQCurveView::EQCurveView(const EQCurveView& v)
    : CControl(v)
    , BackColor(v.BackColor)
    , LineColor(v.LineColor)
    , BorderColor(v.BorderColor)
    , level(v.level)
    , byPass(v.byPass)
{
    setWantsIdle(true);
}

void EQCurveView::setParamNorm(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue normValue)
{
    int index  = std::clamp( ((int)tag - yg331::kParamBand01_Used) / yg331::bandSize,                    0, yg331::numBands - 1); // 0 - 19
    int xIndex = std::clamp((((int)tag - yg331::kParamBand01_Used) / yg331::bandSize) - yg331::numBands, 0, yg331::numXover - 1); // 0 - 1
    switch (tag)
    {
        case yg331::kParamBand01_Used: case yg331::kParamBand02_Used: case yg331::kParamBand03_Used: case yg331::kParamBand04_Used: case yg331::kParamBand05_Used:
        case yg331::kParamBand06_Used: case yg331::kParamBand07_Used: case yg331::kParamBand08_Used: case yg331::kParamBand09_Used: case yg331::kParamBand10_Used:
        case yg331::kParamBand11_Used: case yg331::kParamBand12_Used: case yg331::kParamBand13_Used: case yg331::kParamBand14_Used: case yg331::kParamBand15_Used:
        case yg331::kParamBand16_Used: case yg331::kParamBand17_Used: case yg331::kParamBand18_Used: case yg331::kParamBand19_Used: case yg331::kParamBand20_Used:
            pband[index].Used = normValue > 0.5 ? 1 : 0; break;
            
        case yg331::kParamBand01_Type: case yg331::kParamBand02_Type: case yg331::kParamBand03_Type: case yg331::kParamBand04_Type: case yg331::kParamBand05_Type:
        case yg331::kParamBand06_Type: case yg331::kParamBand07_Type: case yg331::kParamBand08_Type: case yg331::kParamBand09_Type: case yg331::kParamBand10_Type:
        case yg331::kParamBand11_Type: case yg331::kParamBand12_Type: case yg331::kParamBand13_Type: case yg331::kParamBand14_Type: case yg331::kParamBand15_Type:
        case yg331::kParamBand16_Type: case yg331::kParamBand17_Type: case yg331::kParamBand18_Type: case yg331::kParamBand19_Type: case yg331::kParamBand20_Type:
            pband[index].Type = yg331::paramType.ToPlain(normValue); break;
            
        case yg331::kParamBand01_Freq: case yg331::kParamBand02_Freq: case yg331::kParamBand03_Freq: case yg331::kParamBand04_Freq: case yg331::kParamBand05_Freq:
        case yg331::kParamBand06_Freq: case yg331::kParamBand07_Freq: case yg331::kParamBand08_Freq: case yg331::kParamBand09_Freq: case yg331::kParamBand10_Freq:
        case yg331::kParamBand11_Freq: case yg331::kParamBand12_Freq: case yg331::kParamBand13_Freq: case yg331::kParamBand14_Freq: case yg331::kParamBand15_Freq:
        case yg331::kParamBand16_Freq: case yg331::kParamBand17_Freq: case yg331::kParamBand18_Freq: case yg331::kParamBand19_Freq: case yg331::kParamBand20_Freq:
            pband[index].Freq = yg331::paramFreq.ToPlain(normValue); break;
            
        case yg331::kParamBand01_Gain: case yg331::kParamBand02_Gain: case yg331::kParamBand03_Gain: case yg331::kParamBand04_Gain: case yg331::kParamBand05_Gain:
        case yg331::kParamBand06_Gain: case yg331::kParamBand07_Gain: case yg331::kParamBand08_Gain: case yg331::kParamBand09_Gain: case yg331::kParamBand10_Gain:
        case yg331::kParamBand11_Gain: case yg331::kParamBand12_Gain: case yg331::kParamBand13_Gain: case yg331::kParamBand14_Gain: case yg331::kParamBand15_Gain:
        case yg331::kParamBand16_Gain: case yg331::kParamBand17_Gain: case yg331::kParamBand18_Gain: case yg331::kParamBand19_Gain: case yg331::kParamBand20_Gain:
            pband[index].Gain = yg331::paramGain.ToPlain(normValue); break;
            
        case yg331::kParamBand01_Qlty: case yg331::kParamBand02_Qlty: case yg331::kParamBand03_Qlty: case yg331::kParamBand04_Qlty: case yg331::kParamBand05_Qlty:
        case yg331::kParamBand06_Qlty: case yg331::kParamBand07_Qlty: case yg331::kParamBand08_Qlty: case yg331::kParamBand09_Qlty: case yg331::kParamBand10_Qlty:
        case yg331::kParamBand11_Qlty: case yg331::kParamBand12_Qlty: case yg331::kParamBand13_Qlty: case yg331::kParamBand14_Qlty: case yg331::kParamBand15_Qlty:
        case yg331::kParamBand16_Qlty: case yg331::kParamBand17_Qlty: case yg331::kParamBand18_Qlty: case yg331::kParamBand19_Qlty: case yg331::kParamBand20_Qlty:
            pband[index].Qlty = yg331::paramQlty.ToPlain(normValue); break;
            
        case yg331::kParamBandX1_Used: case yg331::kParamBandX2_Used:
            xover[xIndex].Used = normValue > 0.5 ? 1 : 0; break;
            
        case yg331::kParamBandX1_Pass: case yg331::kParamBandX2_Pass:
            xover[xIndex].Pass = yg331::paramPass.ToPlain(normValue); break;
            
        case yg331::kParamBandX1_Freq: case yg331::kParamBandX2_Freq:
            xover[xIndex].Freq = yg331::paramFreq.ToPlain(normValue); break;
            
        case yg331::kParamBandX1_Xtyp: case yg331::kParamBandX2_Xtyp:
            xover[xIndex].Xtyp = yg331::paramXtyp.ToPlain(normValue); break;
            
        case yg331::kParamBandX1_Ordr: case yg331::kParamBandX2_Ordr:
            xover[xIndex].Ordr = yg331::paramOrdr.ToPlain(normValue); break;
            
        default: break;
    }
    
    for (int bands = 0; bands < yg331::numBands; bands++)
        svf[bands].setSVF(pband[bands].Used, pband[bands].Type, pband[bands].Freq, pband[bands].Gain, pband[bands].Qlty, EQ_SR);
    
    for (int bands = 0; bands < yg331::numXover; bands++)
        svfXover[bands].setXover(xover[bands].Used, xover[bands].Pass, xover[bands].Freq, xover[bands].Xtyp, xover[bands].Ordr, EQ_SR);
}

static constexpr CColor color_01(121, 222, 82 ); // #79de52
static constexpr CColor color_02(78,  144, 221); // #4e90dd
static constexpr CColor color_03(181, 61,  221); // #b53ddd
static constexpr CColor color_04(225, 58,  47 ); // #e13b2f
static constexpr CColor color_05(79,  33,  235); // #4f21eb
static constexpr CColor color_06(205, 253, 84 ); // #cdfd54
static constexpr CColor color_07(86,  183, 249); // #56b7f9
static constexpr CColor color_08(119, 251, 181); // #77fbb5
static constexpr CColor color_09(124, 179, 104); // #7cb368
static constexpr CColor color_10(103, 141, 179); // #678db3
static constexpr CColor color_11(156, 95,  178); // #9c5fb2
static constexpr CColor color_12(179, 91,  88 ); // #b35b58
static constexpr CColor color_13(107, 84,  85 ); // #6b5455
static constexpr CColor color_14(173, 195, 97 ); // #adc361
static constexpr CColor color_15(101, 161, 192); // #65a1c0
static constexpr CColor color_16(113, 194, 158); // #71c29e
static constexpr CColor color_17(131, 159, 120); // #839f78
static constexpr CColor color_18(121, 140, 158); // #798c9e
static constexpr CColor color_19(145, 116, 157); // #91749d
static constexpr CColor color_20(175, 113, 112); // #af7170
static constexpr CColor color_X1(122, 111, 161); // #7a6fa1
static constexpr CColor color_X2(155, 166, 116); // #9ba674
static constexpr CColor color_Line(239, 194, 82);
CColor pallet[yg331::numBands] = {
    color_01, color_02, color_03, color_04, color_05, color_06, color_07, color_08, color_09, color_10,
    color_11, color_12, color_13, color_14, color_15, color_16, color_17, color_18, color_19, color_20};
CColor pallet_xovr[yg331::numXover] = {color_X1, color_X2};

// overrides
void EQCurveView::draw(CDrawContext* pContext) {

    pContext->setLineWidth(1);
    pContext->setFillColor(getBackColor());
    pContext->setFrameColor(getBorderColor());
    pContext->drawRect(getViewSize(), VSTGUI::kDrawFilledAndStroked);

    // Given frequency, return screen x position
    auto freq_to_x = [this](double width, double freq) -> double {
        return width * log(freq / MIN_FREQ) / FREQ_LOG_MAX;
    };

    // Given screen x position, return frequency
    auto x_to_freq = [this](double width, double x) -> double {
        return std::max(std::min(MIN_FREQ * exp(FREQ_LOG_MAX * x / width), MAX_FREQ), MIN_FREQ);
    };

    // Given a magnitude, return y screen position as 0..1 with applied tilt
    auto mag_to_01 = [](double m, double freq) -> double {
        if (m == 0) m = 0.00001;
        if (freq == 0) freq = 1.0;
        return 1.0 - (( ((20.0 * log10(m)) + (4.5 * ((log(freq) / log(2.0)) - (log(1024.0) / log(2.0))))) - ceiling) / (noise_floor - ceiling));
    };

    // Given a magnitude (1.0 .... very small number), return y screen position
    auto mag_to_y = [](double height, double m) -> double {
        if (m == 0) m = 0.00001;
        return (((20.0 * log10(m)) - ceiling) / (noise_floor - ceiling)) * height;
    };

    // Given decibels, return screen y position
    auto db_to_y = [](double height, double dB) -> double {
        return (((dB - ceiling) / (noise_floor - ceiling)) * height);
    };

    // Given screen y position, return decibels
    auto y_to_db = [](double height, double y) -> double {
        return ceiling + ((y / height) * (noise_floor - ceiling));
    };

    auto dB_to_y_EQ = [](double height, double dB) -> double {
        return height * (1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5));
    };
    

    auto border = getBorderColor();
    border.setNormAlpha(0.5);
    
    const VSTGUI::CRect r(getViewSize());
    const double r_width = r.getWidth();
    const double r_height = r.getHeight();

    {
        pContext->setFrameColor(border);
        for (int x = 2; x < 10; x++) {
            VSTGUI::CCoord Hz_10 = freq_to_x(r_width, 10.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_10, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_10, r.top);
            pContext->drawLine(_p1, _p2);
        }
        for (int x = 1; x < 10; x++) {
            VSTGUI::CCoord Hz_100 = freq_to_x(r_width, 100.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_100, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_100, r.top);
            pContext->drawLine(_p1, _p2);
        }
        for (int x = 1; x < 10; x++) {
            VSTGUI::CCoord Hz_1000 = freq_to_x(r_width, 1000.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_1000, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_1000, r.top);
            pContext->drawLine(_p1, _p2);
        }

        for (int x = 1; x < 3; x++) {
            VSTGUI::CCoord Hz_10000 = freq_to_x(r_width, 10000.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_10000, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_10000, r.top);
            pContext->drawLine(_p1, _p2);
        }
    }

    {
        pContext->setFrameColor(border);
        for (int cnt = -(int)DB_EQ_RANGE; cnt < (int)DB_EQ_RANGE; cnt += 5)
        {
            VSTGUI::CCoord dB = dB_to_y_EQ(r_height, cnt);
            const VSTGUI::CPoint _p1(r.left,  r.bottom - dB);
            const VSTGUI::CPoint _p2(r.right, r.bottom - dB);
            pContext->drawLine(_p1, _p2);
        }
    }

    std::vector<double> each_band[yg331::numBands];
    for (int bands = 0; bands < yg331::numBands; bands++)
    {
        VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
        if (EQ_curve)
        {
            VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
            EQ_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
            for (double x = -0.5; x <= r.getWidth() + 1; x+=0.5)
            {
                double tmp = MIN_FREQ * std::exp(FREQ_LOG_MAX * x / r.getWidth());
                double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);
                double dB = 20 * log10(svf[bands].mag_response(freq));
                if (svf[bands].Type == yg331::SVF_Generic::tAllPass)
                    dB = svf[bands].phs_response(freq);

                each_band[bands].push_back(dB);
                
                double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
                
                if (svf[bands].Type == yg331::SVF_Generic::tAllPass)
                    m = (((dB / M_PI) / 2) + 0.5);

                double scy = m * r.getHeight();

                if (byPass) scy = 0.5 * r.getHeight();
                EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
            }
            EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
            EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
            EQ_curve->closeSubpath();

            // pContext->setFrameColor(getLineColor().setNormAlpha(0.5));
            CColor line = pallet[bands];
            // line.setNormAlpha(0.5);
            if (svf[bands].getUsed() == 0) line.setNormAlpha(0.0);
            pContext->setFrameColor(line);
            pContext->setDrawMode(VSTGUI::kAntiAliasing);
            pContext->setLineWidth(1.0);
            pContext->setLineStyle(VSTGUI::kLineSolid);
            if (svf[bands].Type == yg331::SVF_Generic::tAllPass)
            {
                pContext->setLineStyle(VSTGUI::kLineOnOffDash);
            }
            pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathStroked);
            EQ_curve->forget();
        }
    }
    for (int bands = 0; bands < yg331::numBands; bands++)
    {
        VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
        if (EQ_curve)
        {
            VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
            EQ_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
            for (double x = -0.5, i = 0; x <= r.getWidth() + 1; x+=0.5, i++)
            {
                double tmp = MIN_FREQ * std::exp(FREQ_LOG_MAX * x / r.getWidth());
                double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);
                double dB = each_band[bands][i];
                
                double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
                
                if (svf[bands].Type == yg331::SVF_Generic::tAllPass)
                    m = (((dB / M_PI) / 2) + 0.5);

                double scy = m * r.getHeight();

                if (byPass) scy = 0.5 * r.getHeight();
                EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
            }
            EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, y_mid));
            EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, y_mid));
            EQ_curve->closeSubpath();

            // pContext->setFrameColor(getLineColor().setNormAlpha(0.5));
            CColor line = pallet[bands];
            // line.setNormAlpha(0.5);
            if (svf[bands].getUsed() == 0) line.setNormAlpha(0.0);
            pContext->setFrameColor(line);
            line.setNormAlpha(0.2);
            if (svf[bands].getUsed() == 0) line.setNormAlpha(0.0);
            pContext->setFillColor(line);
            pContext->setDrawMode(VSTGUI::kAntiAliasing);
            pContext->setLineWidth(1.0);
            pContext->setLineStyle(VSTGUI::kLineSolid);
            pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathFilled);
            EQ_curve->forget();
        }
    }
    
    std::vector<double> xovr_band[yg331::numXover];
    for (int bands = 0; bands < yg331::numXover; bands++)
    {
        VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
        if (EQ_curve)
        {
            VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
            EQ_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
            for (double x = -0.5; x <= r.getWidth() + 1; x+=0.5)
            {
                double tmp = MIN_FREQ * std::exp(FREQ_LOG_MAX * x / r.getWidth());
                double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);
                double dB = 20 * log10(svfXover[bands].mag_response(freq));

                xovr_band[bands].push_back(dB);
                
                double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
                
                if (svf[bands].Type == yg331::SVF_Generic::tAllPass)
                    m = (((dB / M_PI) / 2) + 0.5);

                double scy = m * r.getHeight();

                if (byPass) scy = 0.5 * r.getHeight();
                EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
            }
            EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
            EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
            EQ_curve->closeSubpath();

            // pContext->setFrameColor(getLineColor().setNormAlpha(0.5));
            CColor line = pallet_xovr[bands];
            // line.setNormAlpha(0.5);
            if (svfXover[bands].getUsed() == 0) line.setNormAlpha(0.0);
            pContext->setFrameColor(line);
            pContext->setDrawMode(VSTGUI::kAntiAliasing);
            pContext->setLineWidth(1.0);
            pContext->setLineStyle(VSTGUI::kLineSolid);
            if (svf[bands].Type == yg331::SVF_Generic::tAllPass)
            {
                pContext->setLineStyle(VSTGUI::kLineOnOffDash);
            }
            pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathStroked);
            EQ_curve->forget();
        }
    }
    for (int bands = 0; bands < yg331::numXover; bands++)
    {
        VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
        if (EQ_curve)
        {
            VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
            EQ_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
            for (double x = -0.5, i = 0; x <= r.getWidth() + 1; x+=0.5, i++)
            {
                double tmp = MIN_FREQ * std::exp(FREQ_LOG_MAX * x / r.getWidth());
                double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);
                double dB = xovr_band[bands][i];
                
                double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);

                double scy = m * r.getHeight();

                if (byPass) scy = 0.5 * r.getHeight();
                EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
            }
            EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, y_mid));
            EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, y_mid));
            EQ_curve->closeSubpath();

            // pContext->setFrameColor(getLineColor().setNormAlpha(0.5));
            CColor line = pallet_xovr[bands];
            // line.setNormAlpha(0.5);
            if (svfXover[bands].getUsed() == 0) line.setNormAlpha(0.0);
            pContext->setFrameColor(line);
            line.setNormAlpha(0.2);
            if (svfXover[bands].getUsed() == 0) line.setNormAlpha(0.0);
            pContext->setFillColor(line);
            pContext->setDrawMode(VSTGUI::kAntiAliasing);
            pContext->setLineWidth(1.0);
            pContext->setLineStyle(VSTGUI::kLineSolid);
            pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathFilled);
            EQ_curve->forget();
        }
    }
    
    
    
    
    VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
    if (EQ_curve)
    {
        VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
        // Start from back
        EQ_curve->beginSubpath(VSTGUI::CPoint(r.right + 1, y_mid));
        for (double x = r.getWidth() + 1; x >= -0.5; x-=0.5)
        {
            double dB = 0.0;
            for (int bands = 0; bands < yg331::numBands; bands++)
            {
                if (svf[bands].Type != yg331::SVF_Generic::tAllPass)
                    dB += each_band[bands].back();
                each_band[bands].pop_back();
            }
            for (int bands = 0; bands < yg331::numXover; bands++)
            {
                dB += xovr_band[bands].back();
                xovr_band[bands].pop_back();
            }

            double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
            double scy = m * r.getHeight();

            if (byPass) scy = 0.5 * r.getHeight();
            EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
        }
        EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
        EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
        EQ_curve->closeSubpath();

        // pContext->setFrameColor(getLineColor());
        pContext->setFrameColor(color_Line);
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(2.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathStroked);
        EQ_curve->forget();
    }

    // box outline
    pContext->setLineWidth(1);
    pContext->setFrameColor(getBorderColor());
    pContext->drawRect(getViewSize(), VSTGUI::kDrawStroked);

    setDirty(false);
};

bool EQCurveView::sizeToFit() {
    if (getDrawBackground())
    {
        CRect vs(getViewSize());
        vs.setWidth(getDrawBackground()->getWidth());
        vs.setHeight(getDrawBackground()->getHeight());
        setViewSize(vs);
        setMouseableArea(vs);
        return true;
    }
    return false;
};


static const std::string kAttrBackColor = "back-color";
static const std::string kAttrBorderColor = "border-color";
static const std::string kAttrLineColor = "line-color";

class MyEQCurveViewFactory : public ViewCreatorAdapter
{
public:
    //register this class with the view factory
    MyEQCurveViewFactory() { UIViewFactory::registerViewCreator(*this); }

    //return an unique name here
    IdStringPtr getViewName() const override { return "EQ Curve View"; }

    //return the name here from where your custom view inherites.
    //    Your view automatically supports the attributes from it.
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCControl; }

    //create your view here.
    //    Note you don't need to apply attributes here as
    //    the apply method will be called with this new view
    CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
    {
        return new EQCurveView(CRect(0, 0, 100, 20), nullptr, -1, nullptr);
    }
    bool apply(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) const override
    {
        auto* v = dynamic_cast<EQCurveView*> (view);

        if (!v)
            return false;

        CColor color;
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBackColor), color, description))
            v->setBackColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBorderColor), color, description))
            v->setBorderColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrLineColor), color, description))
            v->setLineColor(color);

        return true;
    }

    bool getAttributeNames(StringList& attributeNames) const override
    {
        attributeNames.emplace_back(kAttrBackColor);
        attributeNames.emplace_back(kAttrBorderColor);
        attributeNames.emplace_back(kAttrLineColor);
        return true;
    }

    AttrType getAttributeType(const std::string& attributeName) const override
    {
        if (attributeName == kAttrBackColor)
            return kColorType;
        if (attributeName == kAttrBorderColor)
            return kColorType;
        if (attributeName == kAttrLineColor)
            return kColorType;
        return kUnknownType;
    }

    //------------------------------------------------------------------------
    bool getAttributeValue(
        CView* view,
        const string& attributeName,
        string& stringValue,
        const IUIDescription* desc) const override
    {
        auto* v = dynamic_cast<EQCurveView*> (view);

        if (!v)
            return false;

        if (attributeName == kAttrBackColor)
        {
            UIViewCreator::colorToString(v->getBackColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrBorderColor)
        {
            UIViewCreator::colorToString(v->getBorderColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrLineColor)
        {
            UIViewCreator::colorToString(v->getLineColor(), stringValue, desc);
            return true;
        }

        return false;
    }
};

MyEQCurveViewFactory __gMyEQCurveViewFactory;

}

namespace yg331 {
//------------------------------------------------------------------------
// LogRangeParameter Declaration
//------------------------------------------------------------------------
class LogRangeParameter_noUnit : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    LogRangeParameter_noUnit (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    
    /** Converts a normalized value to plain value (e.g. 0.5 to 10000.0Hz). */
    Vst::ParamValue toPlain(Vst::ParamValue _valueNormalized) const SMTG_OVERRIDE;
    
    /** Converts a plain value to a normalized value (e.g. 10000 to 0.5). */
    Vst::ParamValue toNormalized(Vst::ParamValue plainValue) const SMTG_OVERRIDE;
    
    /** Converts a normalized value to a string. */
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
    
    OBJ_METHODS (LogRangeParameter_noUnit, RangeParameter)
};
//------------------------------------------------------------------------
// LogRangeParameter Implementation
//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter_noUnit::toPlain(Vst::ParamValue _valueNormalized) const
{
    double FREQ_LOG_MAX = std::log(getMax() / getMin());
    double tmp = getMin() * std::exp(FREQ_LOG_MAX * _valueNormalized);
    double freq = (std::max)((std::min)(tmp, getMax()), getMin());
    return freq;
    //return _valueNormalized * (getMax() - getMin()) + getMin();
}

//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter_noUnit::toNormalized(Vst::ParamValue plainValue) const
{
    SMTG_ASSERT(getMax() - getMin() != 0);
    double FREQ_LOG_MAX = std::log(getMax() / getMin());
    return std::log(plainValue / getMin()) / FREQ_LOG_MAX;
    //return (plainValue - getMin()) / (getMax() - getMin());
}

void LogRangeParameter_noUnit::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
            // wrapper.append(STR16(" "));
            // wrapper.append(getInfo().units);
        }
    }
}

//------------------------------------------------------------------------
// LinRangeParameter Declaration
//------------------------------------------------------------------------
class LinRangeParameter : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    LinRangeParameter (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
};
//------------------------------------------------------------------------
// LinRangeParameter Implementation
//------------------------------------------------------------------------
void LinRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
            wrapper.append(STR16(" "));
            wrapper.append(getInfo().units);
        }
    }
}

//------------------------------------------------------------------------
// LinRangeParameter Declaration
//------------------------------------------------------------------------
class LinRangeParameter_noUnit : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    LinRangeParameter_noUnit (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
};
//------------------------------------------------------------------------
// LinRangeParameter Implementation
//------------------------------------------------------------------------
void LinRangeParameter_noUnit::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
        }
    }
}

//------------------------------------------------------------------------
// GNRC_EQ_Controller Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::initialize (FUnknown* context)
{
    // Here the Plug-in will be instantiated

    //---do not forget to call parent ------
    tresult result = EditControllerEx1::initialize (context);
    if (result != kResultOk)
    {
        return result;
    }

    // Here you could register some parameters

    int32 stepCount;
    int32 flags;
    int32 tag;
    Vst::ParamValue defaultVal;
    Vst::ParamValue defaultPlain;
    Vst::ParamValue minPlain;
    Vst::ParamValue maxPlain;

    tag = kParamBypass;
    stepCount = 1;
    defaultVal = dftBypass ? 1 : 0;
    flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsBypass;
    parameters.addParameter(STR16("Bypass"), nullptr, stepCount, defaultVal, flags, tag);

    flags = Vst::ParameterInfo::kCanAutomate;
    
    tag = kParamLevel;
    stepCount = 0;
    auto* ParamLevel = new LinRangeParameter(STR16("Level"), tag, STR16("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
    ParamLevel->setPrecision(1);
    parameters.addParameter(ParamLevel);
    
    tag = kParamPhase;
    stepCount = 1;
    defaultVal = 0;
    parameters.addParameter(STR16("Phase"), nullptr, stepCount, defaultVal, flags, tag);
    
    flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;
    
    auto* Band_Target = new Vst::StringListParameter(STR16("Target"), kParamTarget, STR16(""), flags);
    for (int i = 0; i < OS_size; i++)
        Band_Target->appendString(target_SR_names[i]);
    Band_Target->getInfo().defaultNormalizedValue = 0.0;
    parameters.addParameter(Band_Target);
    
    for (int bands = 0; bands < numBands; bands++)
    {
        UString128 base;
        base.assign("Band ");
        UString128 bandNumber;
        bandNumber.printInt(bands+1);
        base.append(bandNumber);
        
        Vst::UnitInfo unitInfo;
        Vst::Unit* unit; // create a unit for each bands
        unitInfo.id = bands+1;
        unitInfo.parentUnitId = Steinberg::Vst::kRootUnitId;    // attached to the root unit
        Steinberg::UString (unitInfo.name, USTRINGSIZE (unitInfo.name)).assign (base);
        unitInfo.programListId = Steinberg::Vst::kNoProgramListId;
        unit = new Vst::Unit (unitInfo);
        addUnit (unit);
        
        
        stepCount = 1;
        defaultVal = 0;
        flags = Vst::ParameterInfo::kCanAutomate;
        
        UString128 title_Used;
        title_Used.assign(base);
        title_Used.append(USTRING(" Used"));
        parameters.addParameter(title_Used, nullptr, stepCount, defaultVal, flags, bands * bandSize + kParamBand01_Used, unitInfo.id);
        
        
        flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;
        
        UString128 title_Type;
        title_Type.assign(base);
        title_Type.append(USTRING(" Type"));
        auto* Band_Type = new Vst::StringListParameter(title_Type, bands * bandSize + kParamBand01_Type, STR16(""), flags);
        for (int i = 0; i < SVF_Generic::tSize; i++)
            Band_Type->appendString(SVF_Generic::Filter_Types[i]);
        Band_Type->getInfo().defaultNormalizedValue = nrmParamType;
        Band_Type->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Type);
        
        
        stepCount = 0;
        flags = Vst::ParameterInfo::kCanAutomate;
        
        UString128 title_Freq;
        title_Freq.assign(base);
        title_Freq.append(USTRING(" Freq"));
        auto* Band_Hz = new LogRangeParameter_noUnit(title_Freq, bands * bandSize + kParamBand01_Freq, STR("Hz"), minParamFreq, maxParamFreq, dftBandFreq[bands], stepCount, flags);
        Band_Hz->setPrecision(1);
        Band_Hz->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Hz);
        
        UString128 title_Gain;
        title_Gain.assign(base);
        title_Gain.append(USTRING(" Gain"));
        auto* Band_dB = new LinRangeParameter_noUnit(title_Gain, bands * bandSize + kParamBand01_Gain, STR("dB"), minParamGain, maxParamGain, dftParamGain, stepCount, flags);
        Band_dB->setPrecision(2);
        Band_dB->setUnitID(unitInfo.id);
        parameters.addParameter(Band_dB);
        
        UString128 title_Qlty;
        title_Qlty.assign(base);
        title_Qlty.append(USTRING(" Qlty"));
        auto* Band_Q = new LogRangeParameter_noUnit(title_Qlty, bands * bandSize + kParamBand01_Qlty, STR16("Q"), minParamQlty, maxParamQlty, dftParamQlty, stepCount, flags);
        Band_Q->setPrecision(1);
        Band_Q->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Q);
    }
    
    for (int count = 0; count < numXover; count++)
    {
        UString128 base;
        base.assign("Xover ");
        UString128 bandNumber;
        bandNumber.printInt(count+1);
        base.append(bandNumber);
        
        Vst::UnitInfo unitInfo;
        Vst::Unit* unit; // create a unit for each bands
        unitInfo.id = numBands + count + 1;
        unitInfo.parentUnitId = Steinberg::Vst::kRootUnitId;    // attached to the root unit
        Steinberg::UString (unitInfo.name, USTRINGSIZE (unitInfo.name)).assign (base);
        unitInfo.programListId = Steinberg::Vst::kNoProgramListId;
        unit = new Vst::Unit (unitInfo);
        addUnit (unit);
        
        stepCount = 1;
        defaultVal = 0;
        flags = Vst::ParameterInfo::kCanAutomate;
        
        UString128 title_Used;
        title_Used.assign(base);
        title_Used.append(USTRING(" Used"));
        parameters.addParameter(title_Used, nullptr, stepCount, defaultVal, flags, count * bandSize + kParamBandX1_Used, unitInfo.id);
        
        flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;
        
        UString128 title_Pass;
        title_Pass.assign(base);
        title_Pass.append(USTRING(" Pass"));
        auto* Band_Pass = new Vst::StringListParameter(title_Pass, count * bandSize + kParamBandX1_Pass, STR16(""), flags);
        for (int i = 0; i < SVF_xover::pSize; i++)
            Band_Pass->appendString(SVF_xover::Pass_Types[i]);
        Band_Pass->getInfo().defaultNormalizedValue = nrmParamType;
        Band_Pass->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Pass);
        
        stepCount = 0;
        flags = Vst::ParameterInfo::kCanAutomate;
        
        UString128 title_Freq;
        title_Freq.assign(base);
        title_Freq.append(USTRING(" Freq"));
        auto* Band_Freq = new LogRangeParameter_noUnit(title_Freq, count * bandSize + kParamBandX1_Freq, STR("Hz"), minParamFreq, maxParamFreq, dftBandFreq[count], stepCount, flags);
        Band_Freq->setPrecision(1);
        Band_Freq->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Freq);
        
        flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;
        
        UString128 title_Type;
        title_Type.assign(base);
        title_Type.append(USTRING(" Type"));
        auto* Band_Type = new Vst::StringListParameter(title_Type, count * bandSize + kParamBandX1_Xtyp, STR16(""), flags);
        for (int i = 0; i < SVF_xover::tSize; i++)
            Band_Type->appendString(SVF_xover::Filter_Types[i]);
        Band_Type->getInfo().defaultNormalizedValue = nrmParamType;
        Band_Type->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Type);
        
        UString128 title_Ordr;
        title_Ordr.assign(base);
        title_Ordr.append(USTRING(" Order"));
        auto* Band_Ordr = new Vst::StringListParameter(title_Ordr, count * bandSize + kParamBandX1_Ordr, STR16(""), flags);
        for (int i = 0; i < SVF_xover::oSize; i++)
            Band_Ordr->appendString(SVF_xover::Filter_Order[i]);
        Band_Ordr->getInfo().defaultNormalizedValue = nrmParamType;
        Band_Ordr->setUnitID(unitInfo.id);
        parameters.addParameter(Band_Ordr);
    }
    
    // GUI only parameter
    Vst::StringListParameter* zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
    for (ZoomFactorVector::const_iterator it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
    {
        zoomParameter->appendString(it->title);
    }
    zoomParameter->setNormalized(zoomParameter->toNormalized(dftZoom));
    zoomParameter->getInfo().defaultNormalizedValue = zoomParameter->toNormalized(dftZoom);
    zoomParameter->addDependent(this);
    uiParameters.addParameter(zoomParameter);
    
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!

    getParameterObject(kParamZoom)->removeDependent(this);
    
    //---do not forget to call parent ------
    return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::setComponentState (IBStream* state)
{
    // Here you get the state of the component (Processor part)
    // fprintf (stdout, "GNRC_EQ_Controller::setComponentState\n");
    
    if (!state)
        return kResultFalse;
    
    IBStreamer streamer(state, kLittleEndian);

    // 1. Read Plain Values
    int32           savedBypass = 0;
    // Vst::ParamValue savedZoom   = 0.0;
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
        pBand[bands][bandUsed] =                        savedBand[bands].Used;
        pBand[bands][bandType] = paramType.ToNormalized(savedBand[bands].Type);
        pBand[bands][bandFreq] = paramFreq.ToNormalized(savedBand[bands].Freq);
        pBand[bands][bandGain] = paramGain.ToNormalized(savedBand[bands].Gain);
        pBand[bands][bandQlty] = paramQlty.ToNormalized(savedBand[bands].Qlty);
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        pXovr[bands][bandUsed] =                        savedXovr[bands].Used;
        pXovr[bands][bandPass] = paramPass.ToNormalized(savedXovr[bands].Pass);
        pXovr[bands][bandFreq] = paramFreq.ToNormalized(savedXovr[bands].Freq);
        pXovr[bands][bandXtyp] = paramXtyp.ToNormalized(savedXovr[bands].Xtyp);
        pXovr[bands][bandOrdr] = paramOrdr.ToNormalized(savedXovr[bands].Ordr);
    }
    
    // 3. Set Parameters
    setParamNormalized(kParamBypass, bBypass ? 1 : 0);
    // setParamNormalized(kParamZoom, fZoom);
    setParamNormalized(kParamLevel,  fLevel);
    setParamNormalized(kParamPhase,  bPhase ? 1 : 0);
    setParamNormalized(kParamTarget, fTarget);
    
    for (int bands = 0; bands < numBands; bands++)
    {
        setParamNormalized(kParamBand01_Used + bands*bandSize, pBand[bands][bandUsed]);
        setParamNormalized(kParamBand01_Type + bands*bandSize, pBand[bands][bandType]);
        setParamNormalized(kParamBand01_Freq + bands*bandSize, pBand[bands][bandFreq]);
        setParamNormalized(kParamBand01_Gain + bands*bandSize, pBand[bands][bandGain]);
        setParamNormalized(kParamBand01_Qlty + bands*bandSize, pBand[bands][bandQlty]);
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        setParamNormalized(kParamBandX1_Used + bands*bandSize, pXovr[bands][bandUsed]);
        setParamNormalized(kParamBandX1_Pass + bands*bandSize, pXovr[bands][bandPass]);
        setParamNormalized(kParamBandX1_Freq + bands*bandSize, pXovr[bands][bandFreq]);
        setParamNormalized(kParamBandX1_Xtyp + bands*bandSize, pXovr[bands][bandXtyp]);
        setParamNormalized(kParamBandX1_Ordr + bands*bandSize, pXovr[bands][bandOrdr]);
    }

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::setState (IBStream* state)
{
    // Here you get the state of the controller
    // fprintf (stdout, "GNRC_EQ_Controller::setState\n");
    if (!state)
        return kResultFalse;

    IBStreamer streamer(state, kLittleEndian);

    // 1. Read Plain Values
    Vst::ParamValue savedZoom   = dftZoom/zoomNum;
    Vst::ParamValue savedLevel  = 0.0;
    int32           savedPhase  = 0;
    int32           savedTarget = 0;
    
    bandParamSet savedBand[numBands];
    xovrParamSet savedXovr[numBands];

    if (streamer.readDouble(savedZoom  ) == false) savedZoom = dftZoom/zoomNum;
    if (streamer.readDouble(savedLevel ) == false) savedLevel = fLevel;
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
    fZoom   = savedZoom;
    fLevel  = paramGain.ToNormalized(savedLevel);;
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
    
    // 3. Set Parameters
    setParamNormalized(kParamZoom,  fZoom);
    setParamNormalized(kParamLevel, fLevel);
    setParamNormalized(kParamPhase, bPhase ? 1 : 0);
    setParamNormalized(kParamTarget, fTarget);
    
    for (int bands = 0; bands < numBands; bands++)
    {
        setParamNormalized(kParamBand01_Used  + bands*bandSize, pBand[bands][bandUsed]);
        setParamNormalized(kParamBand01_Type  + bands*bandSize, pBand[bands][bandType]);
        setParamNormalized(kParamBand01_Freq  + bands*bandSize, pBand[bands][bandFreq]);
        setParamNormalized(kParamBand01_Gain  + bands*bandSize, pBand[bands][bandGain]);
        setParamNormalized(kParamBand01_Qlty  + bands*bandSize, pBand[bands][bandQlty]);
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        setParamNormalized(kParamBandX1_Used + bands*bandSize, pXovr[bands][bandUsed]);
        setParamNormalized(kParamBandX1_Pass + bands*bandSize, pXovr[bands][bandPass]);
        setParamNormalized(kParamBandX1_Freq + bands*bandSize, pXovr[bands][bandFreq]);
        setParamNormalized(kParamBandX1_Xtyp + bands*bandSize, pXovr[bands][bandXtyp]);
        setParamNormalized(kParamBandX1_Ordr + bands*bandSize, pXovr[bands][bandOrdr]);
    }

    return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::getState (IBStream* state)
{
    // Here you are asked to deliver the state of the controller (if needed)
    // Note: the real state of your plug-in is saved in the processor
    // fprintf (stdout, "GNRC_EQ_Controller::getState\n");
    if (!state)
        return kResultFalse;

    IBStreamer streamer(state, kLittleEndian);

    fZoom   = getParamNormalized(kParamZoom);
    fLevel  = getParamNormalized(kParamLevel);
    bPhase  = getParamNormalized(kParamPhase);
    fTarget = getParamNormalized(kParamTarget);

    for (int bands = 0; bands < numBands; bands++)
    {
        pBand[bands][bandUsed] = getParamNormalized(kParamBand01_Used  + bands*bandSize);
        pBand[bands][bandType] = getParamNormalized(kParamBand01_Type  + bands*bandSize);
        pBand[bands][bandFreq] = getParamNormalized(kParamBand01_Freq  + bands*bandSize);
        pBand[bands][bandGain] = getParamNormalized(kParamBand01_Gain  + bands*bandSize);
        pBand[bands][bandQlty] = getParamNormalized(kParamBand01_Qlty  + bands*bandSize);
    }
    for (int bands = 0; bands < numXover; bands++)
    {
        pXovr[bands][bandUsed] = getParamNormalized(kParamBandX1_Used  + bands*bandSize);
        pXovr[bands][bandPass] = getParamNormalized(kParamBandX1_Pass  + bands*bandSize);
        pXovr[bands][bandFreq] = getParamNormalized(kParamBandX1_Freq  + bands*bandSize);
        pXovr[bands][bandXtyp] = getParamNormalized(kParamBandX1_Xtyp  + bands*bandSize);
        pXovr[bands][bandOrdr] = getParamNormalized(kParamBandX1_Ordr  + bands*bandSize);
    }

    streamer.writeDouble(fZoom);
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

    return kResultTrue;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API GNRC_EQ_Controller::createView (FIDString name)
{
    // Here the Host wants to open your editor (if you have one)
    if (FIDStringsEqual (name, Vst::ViewType::kEditor))
    {
        // create your editor here and return a IPlugView ptr of it
        auto* view = new VSTGUI::VST3Editor (this, "Main", "GNRC_EQ_editor.uidesc");
        std::vector<double> _zoomFactors;
        _zoomFactors.push_back(0.50);
        _zoomFactors.push_back(0.75);
        _zoomFactors.push_back(1.00);
        _zoomFactors.push_back(1.25);
        _zoomFactors.push_back(1.50);
        _zoomFactors.push_back(1.75);
        _zoomFactors.push_back(2.00);
        view->setAllowedZoomFactors(_zoomFactors);
        view->setZoomFactor(1.0);
        view->setIdleRate(1000.0/60.0);
        setKnobMode(Steinberg::Vst::KnobModes::kLinearMode);
        return view;
    }
    return nullptr;
}

VSTGUI::IController* GNRC_EQ_Controller::createSubController (VSTGUI::UTF8StringPtr name,
                                          const VSTGUI::IUIDescription* description,
                                          VSTGUI::VST3Editor* editor)
{
    if (VSTGUI::UTF8StringView(name) == "EQCurveViewController")
    {
        EQCurveViewController* controller = new EQCurveViewController(editor, this);
        for (int bands = 0; bands < numBands; bands++)
        {
            controller->addBandParam(getParameterObject(kParamBand01_Used + bands*bandSize),
                                     getParameterObject(kParamBand01_Type + bands*bandSize),
                                     getParameterObject(kParamBand01_Freq + bands*bandSize),
                                     getParameterObject(kParamBand01_Gain + bands*bandSize),
                                     getParameterObject(kParamBand01_Qlty + bands*bandSize) );
        }
        for (int bands = 0; bands < numXover; bands++)
        {
            controller->addXovrParam(getParameterObject(kParamBandX1_Used + bands*bandSize),
                                     getParameterObject(kParamBandX1_Pass + bands*bandSize),
                                     getParameterObject(kParamBandX1_Freq + bands*bandSize),
                                     getParameterObject(kParamBandX1_Xtyp + bands*bandSize),
                                     getParameterObject(kParamBandX1_Ordr + bands*bandSize) );
        }
        controller->addLevelParam(getParameterObject(kParamLevel));
        controller->addBypassParam(getParameterObject(kParamBypass));
        addEQCurveViewController(controller);
        return controller;
    }
    return nullptr;
};

//------------------------------------------------------------------------
void GNRC_EQ_Controller::editorAttached(Vst::EditorView* editor)
{
    editors.push_back(editor);
}

//------------------------------------------------------------------------
void GNRC_EQ_Controller::editorRemoved(Vst::EditorView* editor)
{
    editors.erase(std::find(editors.begin(), editors.end(), editor));
}

void PLUGIN_API GNRC_EQ_Controller::update(FUnknown* changedUnknown, int32 message)
{
    EditControllerEx1::update(changedUnknown, message);

    // GUI Resizing
    // check 'zoomtest' code at
    // https://github.com/steinbergmedia/vstgui/tree/vstgui4_10/vstgui/tests/uidescription%20vst3/source

    Vst::Parameter* param = FCast<Vst::Parameter>(changedUnknown);
    if (!param)
        return;

    if (param->getInfo().id == kParamZoom)
    {
        size_t index = static_cast<size_t> (param->toPlain(param->getNormalized()));

        if (index >= zoomFactors.size())
            return;

        for (EditorVector::const_iterator it = editors.begin(), end = editors.end(); it != end; ++it)
        {
            VSTGUI::VST3Editor* editor = dynamic_cast<VSTGUI::VST3Editor*>(*it);
            if (editor)
                editor->setZoomFactor(zoomFactors[index].factor);
        }
    }
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::setParamNormalized (Vst::ParamID tag, Vst::ParamValue value)
{
    // called by host to update your parameters
    tresult result = EditControllerEx1::setParamNormalized (tag, value);
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::getParamStringByValue (Vst::ParamID tag, Vst::ParamValue valueNormalized, Vst::String128 string)
{
    // called by host to get a string for given normalized value of a specific parameter
    // (without having to set the value!)
    return EditControllerEx1::getParamStringByValue (tag, valueNormalized, string);
}

//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::getParamValueByString (Vst::ParamID tag, Vst::TChar* string, Vst::ParamValue& valueNormalized)
{
    // called by host to get a normalized value from a string representation of a specific parameter
    // (without having to set the value!)
    return EditControllerEx1::getParamValueByString (tag, string, valueNormalized);
}

//------------------------------------------------------------------------
tresult GNRC_EQ_Controller::receiveText(const char* text)
{
    // received from Component
    if (text)
    {
        if (strcmp("latency_changed\n", text))
        {
            if (auto componentHandler = getComponentHandler ())
            {
                componentHandler->restartComponent (Vst::kLatencyChanged); // if Vst::kReloadComponent, it crashes
                // fprintf (stdout, "restartComponent\n");
            }
        }
    }
    return kResultOk;
}

//------------------------------------------------------------------------
// DataExchangeController Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API GNRC_EQ_Controller::notify(Vst::IMessage* message)
{
    if (!message)
        return kInvalidArgument;
    
    if (strcmp (message->getMessageID (), "GUI") == 0)
    {
        
        if (!curveControllers.empty())
                {
                    for (auto iter = curveControllers.begin(); iter != curveControllers.end(); iter++)
                    {
                        if (auto attributes = message->getAttributes ())
                        {
                            const void* data;
                            uint32 sizeInBytes;
                            ParamValue getValue = 0.0;
                            
                            if (attributes->getFloat  ("targetSR",  getValue) == kResultTrue)
                            {
                                targetSR  = getValue;
                                (*iter)->setEQsampleRate(targetSR);
                            }
                        }
                    }
                }
            return kResultOk;
    }

    return EditControllerEx1::notify(message);
}

//------------------------------------------------------------------------
} // namespace yg331
