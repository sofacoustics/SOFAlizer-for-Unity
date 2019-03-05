# SOFAlizer for Unity 1.2.0

SOFAlizer is a SOFA-based audio spatializer for Unity. It is a Unity [native plugin](https://docs.unity3d.com/Manual/NativePlugins.html). It loads sets of spatial filters so-called [HRTFs](https://en.wikipedia.org/wiki/Head-related_transfer_function) stored in [SOFA](https://www.sofaconventions.org/) files and renders spatial audio based on a selected HRTF set. The output is a binaural signal intended to be listened to via headphones. 

Currently only Windows is supported.

Installation:
-------------

Compile SOFAlizer DLL:
* Start NativeCode\VisualStudio\AudioPluginSOFA.sln
* Select x64, Release
* Rebuild AudioPluginSOFA 

Start the scene:
* In Unity start the scene: Assets\SampleScene.unity
* Select SOFAlizer as spatializer plug in: Edit --> Project Settings --> Audio --> Spatializer Plugin --> SOFA Spatializer
* Copy some SOFA files (e.g., from https://www.sofaconventions.org/mediawiki/index.php/Files) to the root directory. Rename them to `hrtfX.sofa` with `X` from 0 to 9. 
* Start the engine by clicking on "Run"

Usage:
------
On "Run" of the scene, SOFAlizer:
* Loads SOFA files `hrtfX.sofa` with `X` being from 0 to 9. These files must be located in the root directory of this project. The HRTFs must be stored using the SOFA convention "SimpleFreeFieldHRIR". 
* Resamples each HRTF set to the sampling rate provided in Audio --> Project Settings --> Audio --> System Sample Rate. This may take some time (can be avoided by using HRTFs sampled at the same rate as that from the Unity project). 
* Normalizes each HRTF set relative to the HRTF for the frontal position. 
* Crops the HRTFs to the length of 256 samples. 
* Transforms all impulse responses to the spectral domain for fast convolution in real-time.

In the scene, the parameter "SOFA Selector" selects the index of the processed HRTF set. The selection of the loaded HRTF sets can be done at any time in real-time. If a non-loaded HRTF set is selected, the audio will be muted.

In order to change the loaded HRTFs, stop the scene rendering, replace one or more of the SOFA files, and re-run the scene. 

Log file SOFAlizer.log will be created where the loaded files and their transformations are logged.

**Enjoy!**


Acknowledgements:
-----------------

* Christian Hoene for the great libmysofa library for easy reading SOFA files: https://github.com/hoene/libmysofa

History:
--------

Version 1.2.1
=============
nothing here yet...

Version 1.2.0
=============
* Interpolation between the new and old position for a smooth movement of the listener. 
* Major clean up of the code. 
* Multiple Debug levels: 1: only load/unload, 2: real time information

Version 1.1.0
=============

Contributions from others (multiple instances, crashes, bug fixes)

Version 1.0.1
=============

Bug fix for crash when loading/unloading the plugin. 

Version 1.0.0
=============

First release