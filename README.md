# Generic-EQ

Generic EQ is a minimum phase EQ for REW Auto EQ - Generic type.  

## Features  

20+2 bands to fully support REW EQ.  

| REW                | Generic EQ   |
|--------------------|--------------|
| PK, Modal          | Bell         |
| LP, HP, LP Q, HP Q | LP, HP       |
| LP1, HP1           | LP 6, HP 6   |
| LS, LS Q, HS, HS Q | LS, HS       |
| LS 6dB, HS 6dB     | LS 6, HS 6   |
| LS 12dB, HS 12dB   | LS 12, HS 12 |
| All pass           | AP           |

> Modal controls Q in time domain, so it's technically same as PK.  
> Pass filters are butterworth, meaning Q = 0.7071   
> Shelf 6dB 12dB have fc at -3dB = 0.7071 point, while Shelf and Shelf Q have fc at middle(zero-crossing of phase)  

Exact gain-q dependency to match REW EQ's curve.  

Statefull filter design for better numerical stability.  

Runs in double precision 64-bit internal processing. Also double precision input / output if supported.  

Windows and Mac, VST3 and AU.  

[![GitHub Release](https://img.shields.io/github/v/release/kiriki-liszt/Generic_EQ?style=flat-square&label=Get%20latest%20Release)](https://github.com/Kiriki-liszt/Generic_EQ/releases/latest)
[![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/kiriki-liszt/Generic_EQ/total?style=flat-square&label=total%20downloads&color=blue)](https://tooomm.github.io/github-release-stats/?username=Kiriki-liszt&repository=Generic_EQ)  

[![Static Badge](https://img.shields.io/badge/coffee%20maybe%3F%20%3D%5D%20-gray?style=for-the-badge&logo=buy-me-a-coffee)](https://buymeacoffee.com/kirikiaris)  

<img src="https://github.com/Kiriki-liszt/Generic_EQ/blob/main/screenshot.png?raw=true"  width="600"/>  

## Compatibility  

VST3, AUv2  

## System Requirements

Audio Units  

* Mac OS X 10.13 or later (Intel or Apple Silicon Native)

VST3  

* Mac OS X 10.13 or later (Intel or Apple Silicon Native)
* Windows 10 or later

## How to use  

1. Windows

Unzip Win.zip from latest release and copy to "C:\Program Files\Common Files\VST3".  

2. MacOS(Intel tested, Apple Silicon not tested).  

Unzip MacOS.zip from latest release and copy vst3 to "/Library/Audio/Plug-Ins/VST3" and component to "/Library/Audio/Plug-Ins/Components".  

> If it doesn't go well, configure security options in console as  
>
> ``` console  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/VST3/Generic_EQ.vst3  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/Components/Generic_EQ.component
>
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/VST3/Generic_EQ.vst3  
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/Components/Generic_EQ.component
> ```  
>
> tested by @jonasborneland [here](https://github.com/Kiriki-liszt/JS_Inflator_to_VST2_VST3/issues/12#issuecomment-1616671177)

## Licensing  

Generic EQ is using GPL v3 license.  

### VST  

> Q: I would like to share the source code of my VST 3 plug-in/host on GitHub or other such platform.  
>
> - You can choose the GPLv3 license and feel free to share your plug-ins/host's source code including or referencing the VST 3 SDK's sources on GitHub.  
> - **You are allowed to provide a binary form of your plug-ins/host too, provided that you provide its source code as GPLv3 too.**  
> - Note that you have to follow the Steinberg VST usage guidelines.  
>  
> <https://steinbergmedia.github.io/vst3_dev_portal/pages/FAQ/Licensing.html>  

![VST Logo](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/142e3c12-cd5f-415d-9b72-8b4f04419633)

VSTSDK 3.7.12 used  
VSTGUI 4.14 used  

## Project Build  

Use CMake to build itself or make IDE project file.  
Check .github/workflows/Windows Build.yml.  
Remember to git clone VSTSDK, too.  
Supports Windows, Mac, Linux(same as VSTSDK).  

## Version logs

v1.0.0.b : intial try.  

## Further lookings  

[https://gearspace.com/board/showpost.php?p=15864586&postcount=730](https://gearspace.com/board/showpost.php?p=15864586&postcount=730)  
[Vicanek, Martin. Matched Second Order Digital Filters. (2016).](https://www.vicanek.de/articles/BiquadFits.pdf)  
[John Flynn & Joshua D. Reiss (2018). Improving the frequency response magnitude and phase of analogue-matched digital filters](https://www.eecs.qmul.ac.uk/~josh/documents/2018/19412.pdf)  
[D. W. Gunness, O. S. Chauhan, “Optimizing the Magnitude Response of Matched z-Transform Filters (“MZTi”) for Loudspeaker Equalization”](https://www.khabdha.org/wp-content/uploads/2008/03/optimizing-the-magnitude-response-of-mzt-filters-mzti-2007.pdf)  

[Martin Vicanek. Matched Two-Pole Digital Shelving Filters. (23. March 2024).](https://vicanek.de/articles/2poleShelvingFits.pdf)  
[Martin Vicanek. Matched Second Order Digital Filters. (14. February 2016).](https://www.vicanek.de/articles/BiquadFits.pdf)  

## Ref  

<https://dafx14.fau.de/papers/dafx14_aaron_wishnick_time_varying_filters_for_.pdf>  
<https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_52.pdf>  
<https://cytomic.com/files/dsp/SVF-vs-DF1.pdf>  
<https://www.researchgate.net/publication/282326563>  
<https://www.dsprelated.com/freebooks/filters/Implementation_Structures_Recursive_Digital.html>  
<https://forum.juce.com/t/dsp-module-discussion-iir-filter-and-statevariablefilter/23891>  

## TODO  

[] Decramp using Martin Vicanek method.  
