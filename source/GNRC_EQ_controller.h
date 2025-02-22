//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "GNRC_EQ_shared.h"

#include <array>

#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "vstgui/uidescription/delegationcontroller.h"
#include "vstgui/lib/controls/cknob.h"

namespace VSTGUI {
//------------------------------------------------------------------------
// EQ Curve Display
//------------------------------------------------------------------------
class EQCurveView : public CControl {
public:
    EQCurveView(const CRect& size, IControlListener* listener, int32_t tag, CBitmap* background );
    EQCurveView(const EQCurveView& v);

    // get/set Attributes
    virtual void setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
    CColor getBackColor() const { return BackColor; }

    virtual void setBorderColor(CColor color) { if (BorderColor != color) { BorderColor = color; setDirty(true); } }
    CColor getBorderColor() const { return BorderColor; }

    virtual void setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
    CColor getLineColor() const { return LineColor; }
    
    void setParamNorm(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue normValue);
    void setLevel(Steinberg::Vst::ParamValue normValue) {level = yg331::paramLevl.ToPlain(normValue);}
    void setBypass(Steinberg::Vst::ParamValue normValue) {byPass = normValue == 0.0 ? false : true;}
    void setEQsampleRate(double SR) {
        EQ_SR = SR;
        for (int bands = 0; bands < yg331::numBands; bands++)
            svf[bands].setSVF(pband[bands].Used, pband[bands].Type, pband[bands].Freq, pband[bands].Gain, pband[bands].Qlty, EQ_SR);
        
        for (int bands = 0; bands < yg331::numXover; bands++)
            svfXover[bands].setXover(xover[bands].Used, xover[bands].Pass, xover[bands].Freq, xover[bands].Xtyp, xover[bands].Ordr, EQ_SR);
    }

    // overrides
    void setDirty(bool state) override { CView::setDirty(state); };
    void draw(CDrawContext* pContext) override;
    void setViewSize(const CRect& newSize, bool invalid = true) override { CControl::setViewSize(newSize, invalid); };
    bool sizeToFit() override;

    /** called on idle when view wants idle */
    void onIdle() override {
        invalid();
    };

    CLASS_METHODS(EQCurveView, CControl)

    //------------------------------------------------------------------------

protected:
    ~EQCurveView() noexcept override
    {
    };
    
    static SMTG_CONSTEXPR double MAX_FREQ = 22000.0;
    static SMTG_CONSTEXPR double MIN_FREQ = 10.0;
    double FREQ_LOG_MAX = log(MAX_FREQ / MIN_FREQ);
    static SMTG_CONSTEXPR double ceiling     = 0.0;   // dB
    static SMTG_CONSTEXPR double noise_floor = -72.0; // dB
    static SMTG_CONSTEXPR double DB_EQ_RANGE = 15.0;

    CColor  BackColor;
    CColor  LineColor;
    CColor  BorderColor;
    
    yg331::bandParamSet pband[yg331::numBands];
    yg331::xovrParamSet xover[yg331::numXover];
    
    yg331::SVF_Generic  svf[yg331::numBands];
    yg331::SVF_xover    svfXover[yg331::numXover];

    bool    byPass = false;
    double  level = 0.0;
    double  EQ_SR = 48000.0;
};
}

namespace yg331 {

class EQCurveViewController;

//------------------------------------------------------------------------
//  GNRC_EQ_Controller
//------------------------------------------------------------------------
class GNRC_EQ_Controller
    : public Steinberg::Vst::EditControllerEx1
    , public VSTGUI::VST3EditorDelegate
{
public:
//------------------------------------------------------------------------
    GNRC_EQ_Controller () = default;
    ~GNRC_EQ_Controller () SMTG_OVERRIDE = default;

    // Create function
    static Steinberg::FUnknown* createInstance (void* /*context*/)
    {
        return (Steinberg::Vst::IEditController*)new GNRC_EQ_Controller;
    }

    // IPluginBase
    Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;

    // EditController
    Steinberg::tresult PLUGIN_API setComponentState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::IPlugView* PLUGIN_API createView (Steinberg::FIDString name) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API setParamNormalized (Steinberg::Vst::ParamID tag,
                                                      Steinberg::Vst::ParamValue value) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getParamStringByValue (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::ParamValue valueNormalized,
                                                         Steinberg::Vst::String128 string) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getParamValueByString (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::TChar* string,
                                                         Steinberg::Vst::ParamValue& valueNormalized) SMTG_OVERRIDE;
    // VST3EditorDelegate
    /** verify a view after it was created */
    /*
    VSTGUI::CView* PLUGIN_API verifyView(VSTGUI::CView* view,
                                         const VSTGUI::UIAttributes& attributes,
                                         const VSTGUI::IUIDescription* description,
                                         VSTGUI::VST3Editor* editor) SMTG_OVERRIDE;
    */

    // EditController
    Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API receiveText(const char* text) SMTG_OVERRIDE;
    void PLUGIN_API update(Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE;
    void editorAttached(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE; ///< called from EditorView if it was attached to a parent
    void editorRemoved (Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE; ///< called from EditorView if it was removed from a parent

    //------------------------------------------------------------------------
    Steinberg::Vst::Parameter* getParameterObject(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        Steinberg::Vst::Parameter* param = EditControllerEx1::getParameterObject(tag);
        if (param == 0)
        {
            param = uiParameters.getParameter(tag);
        }
        return param;
    }
    bool isPrivateParameter(const Steinberg::Vst::ParamID paramID) SMTG_OVERRIDE
    {
        return uiParameters.getParameter(paramID) != 0 ? true : false;
    }

    // make sure that our UI only parameters doesn't call the following three EditController methods: beginEdit, endEdit, performEdit
    //------------------------------------------------------------------------
    Steinberg::tresult beginEdit(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::beginEdit(tag);
        return Steinberg::kResultFalse;
    }

    //------------------------------------------------------------------------
    Steinberg::tresult performEdit(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue valueNormalized) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::performEdit(tag, valueNormalized);
        return Steinberg::kResultFalse;
    }

    //------------------------------------------------------------------------
    Steinberg::tresult endEdit(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::endEdit(tag);
        return Steinberg::kResultFalse;
    }

    //---from VST3EditorDelegate-----------
    VSTGUI::IController* createSubController (VSTGUI::UTF8StringPtr name,
                                              const VSTGUI::IUIDescription* description,
                                              VSTGUI::VST3Editor* editor) SMTG_OVERRIDE;
    
    //---Internal functions-------
    void addEQCurveViewController(EQCurveViewController* controller)
    {
        curveControllers.push_back(controller);
    };
    void removeEQCurveViewController(EQCurveViewController* controller)
    {
        auto it = std::find(curveControllers.begin(), curveControllers.end(), controller);
        if (it != curveControllers.end())
            curveControllers.erase(it);
    };
    
    //---Interface---------
    DEFINE_INTERFACES
        // Here you can add more supported VST3 interfaces
        // DEF_INTERFACE (Vst::IXXX)
        DEF_INTERFACE (Steinberg::Vst::IUnitInfo)
    END_DEFINE_INTERFACES (EditController)
    // DELEGATE_REFCOUNT (EditController)
    DELEGATE_REFCOUNT (EditControllerEx1)

//------------------------------------------------------------------------
protected:
    // UI only parameter list
    Steinberg::Vst::ParameterContainer uiParameters;

    // editor list
    typedef std::vector<Steinberg::Vst::EditorView*> EditorVector;
    EditorVector editors;

    using UICurveControllerList = std::vector<EQCurveViewController*>;
    UICurveControllerList curveControllers;

    TBool      bBypass = false;
    ParamValue fLevel  = nrmParamLevl;
    TBool      bPhase  = false;
    ParamValue fTarget = OS_1x;
    ParamValue fZoom   = 2.0 / 6.0;
    ParamValue fStereo = ST_ST;
    
    std::array<std::array<double, bandSize>, numBands> pBand = {{
        {dftParamUsed, nrmParamType, dftBand01Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand02Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand03Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand04Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand05Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand06Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand07Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand08Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand09Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand10Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand11Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand12Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand13Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand14Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand15Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand16Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand17Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand18Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand19Freq, dftParamGain, dftParamQlty},
        {dftParamUsed, nrmParamType, dftBand20Freq, dftParamGain, dftParamQlty}
    }};
    
    std::array<std::array<ParamValue, bandSize>, numXover> pXovr = {{
        {dftParamUsed, nrmParamPass, nrmBand01Freq, nrmParamXtyp, nrmParamOrdr},
        {dftParamUsed, nrmParamPass, nrmBand02Freq, nrmParamXtyp, nrmParamOrdr}
    }};
    
    SampleRate targetSR = 48000.0;
};


//------------------------------------------------------------------------
// EQCurveViewController
//------------------------------------------------------------------------
class EQCurveViewController :
    public Steinberg::FObject,
    public VSTGUI::DelegationController
{
public:
    using Parameter          = Steinberg::Vst::Parameter;
    using ParameterContainer = Steinberg::Vst::ParameterContainer;
    
    EQCurveViewController ( IController* baseController, GNRC_EQ_Controller* mainController ) :
        DelegationController(baseController),
        mainController(mainController),
        eqCurveView(nullptr)
    {
    }

    ~EQCurveViewController() override
    {
        while (!pBand.empty()) {
            pBand.front()->removeDependent(this);
            pBand.erase(pBand.cbegin());
        }
        if (pLevel)  { pLevel ->removeDependent(this); pLevel  = nullptr; }
        if (pBypass) { pBypass->removeDependent(this); pBypass = nullptr; }
        
        mainController->removeEQCurveViewController(this);
    }
    
    void setEQsampleRate(double SR) {eqCurveView->setEQsampleRate(SR);}
    
    void addBandParam(Parameter* pUsed, Parameter* pType, Parameter* pFreq, Parameter* pGain, Parameter* pQlty)
    {
        if (pUsed) { pBand.push_back(pUsed); pUsed->addDependent(this); }
        if (pType) { pBand.push_back(pType); pType->addDependent(this); }
        if (pFreq) { pBand.push_back(pFreq); pFreq->addDependent(this); }
        if (pGain) { pBand.push_back(pGain); pGain->addDependent(this); }
        if (pQlty) { pBand.push_back(pQlty); pQlty->addDependent(this); }
    }
    void addXovrParam(Parameter* pUsed, Parameter* pPass, Parameter* pFreq, Parameter* pXtyp, Parameter* pOrdr)
    {
        if (pUsed) { pBand.push_back(pUsed); pUsed->addDependent(this); }
        if (pPass) { pBand.push_back(pPass); pPass->addDependent(this); }
        if (pFreq) { pBand.push_back(pFreq); pFreq->addDependent(this); }
        if (pXtyp) { pBand.push_back(pXtyp); pXtyp->addDependent(this); }
        if (pOrdr) { pBand.push_back(pOrdr); pOrdr->addDependent(this); }
    }
    void addLevelParam(Parameter* p)
    {
        pLevel = p; p->addDependent(this);
    }
    void addBypassParam(Parameter* p)
    {
        pBypass = p; p->addDependent(this);
    }

private:
    using CControl       = VSTGUI::CControl;
    using CView          = VSTGUI::CView;
    using EQCurveView    = VSTGUI::EQCurveView;
    using UTF8String     = VSTGUI::UTF8String;
    using UIAttributes   = VSTGUI::UIAttributes;
    using IUIDescription = VSTGUI::IUIDescription;

    // FObject
    void PLUGIN_API update (Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE
    {
        if (eqCurveView)
        {
            if (auto* p = Steinberg::FCast<Parameter>(changedUnknown); p)
            {
                Steinberg::Vst::ParamID tag = p->getInfo().id;
                if (message == kChanged)
                {
                    for (auto& elem: pBand)
                    {
                        if (elem->getInfo().id == tag) eqCurveView->setParamNorm(tag, p->getNormalized());
                    }
                    if (p == pLevel  && pLevel)  eqCurveView->setLevel (p->getNormalized());
                    if (p == pBypass && pBypass) eqCurveView->setBypass(p->getNormalized());
                }
                else if (message == kWillDestroy)
                {
                    for (auto elem = pBand.begin(); elem != pBand.end();)
                    {
                        if ((*elem)->getInfo().id == tag) { (*elem)->removeDependent(this); elem = pBand.erase(elem);}
                        else {elem++;}
                    }
                    if (p == pLevel  && pLevel)  { pLevel ->removeDependent(this); pLevel  = nullptr; }
                    if (p == pBypass && pBypass) { pBypass->removeDependent(this); pBypass = nullptr; }
                }
            }
        }
    }

    //--- from IControlListener ----------------------
    // void valueChanged(CControl* pControl) override { }
    // void controlBeginEdit(CControl* /*pControl*/) override {}
    // void controlEndEdit(CControl* pControl) override  { }

    //--- from IControlListener ----------------------
    //--- is called when a view is created -----
    CView* verifyView(
        CView* view,
        const UIAttributes& /*attributes*/,
        const IUIDescription* /*description*/) override
    {
        if (EQCurveView* control = dynamic_cast<EQCurveView*>(view); control)
        {
            eqCurveView = control;
        }
        return view;
    }

    GNRC_EQ_Controller* mainController;
    EQCurveView* eqCurveView;
    
    std::vector<Parameter*> pBand;
    Parameter* pLevel;
    Parameter* pBypass;
};

//------------------------------------------------------------------------
} // namespace yg331
