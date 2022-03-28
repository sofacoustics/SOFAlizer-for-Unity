# SOFAlizer for Unity 1.x

SOFAlizer is a SOFA-based audio spatializer for Unity. It is a Unity [native plugin](https://docs.unity3d.com/Manual/NativePlugins.html). It loads sets of spatial filters so-called [HRTFs](https://en.wikipedia.org/wiki/Head-related_transfer_function) stored in [SOFA](https://www.sofaconventions.org/) files and renders spatial audio based on a selected HRTF set. The output is a binaural signal intended to be listened to via headphones. 

Currently only Windows is supported. Tested with Unity 2021.2.8f1, compiled with MS VS 2019.

Installation:
-------------

Compile SOFAlizer DLL:
* Start NativeCode\VisualStudio\AudioPluginSOFA.sln
* Select x64, Release or Debug
* Rebuild AudioPluginSOFA 
* Adapt or install Windows SDK version if an error occurs during compilation. (Properties)

Start the scene:
* In Unity start the scene: Assets\SampleScene.unity
* Select SOFAlizer as spatializer plug in: Edit --> Project Settings --> Audio --> Spatializer Plugin --> SOFA Spatializer
* Copy some SOFA files (e.g., from https://www.sofaconventions.org/mediawiki/index.php/Files) to the root directory. Rename them to `hrtfX.sofa` with `X` from 0 to 49. 
* Start the engine by clicking on "Run"

Usage:
------
On "Run" of the scene, SOFAlizer:
* Loads SOFA files `hrtfX.sofa` with `X` being from 0 to 49. These files must be located in the root directory of this project. The HRTFs must be stored using the SOFA convention "SimpleFreeFieldHRIR". 
* Resamples each HRTF set to the sampling rate provided in Audio --> Project Settings --> Audio --> System Sample Rate. This may take some time (can be avoided by using HRTFs sampled at the same rate as that from the Unity project). 
* Normalizes each HRTF set relative to the HRTF for the frontal position. 
* Crops the HRTFs to the length of 128 samples. 
* Transforms all impulse responses to the spectral domain for fast convolution in real-time.
* Creates SOFAlizer.log with the information of the loaded HRTF sets. In order to change the loaded HRTF sets, stop the scene rendering, replace one or more of the SOFA files, and re-run the scene.

Scene parameters: 
* SOFA Selector: selects the index of the processed HRTF set. The selection of the loaded HRTF sets can be done at any time in real-time. If a non-loaded HRTF set is selected, the audio will be muted.
* Debug: determines the level of debug information shown in the console (only available if compiled in the debug mode!):
  * 0 ... loading information only
  * 1 ... show the active HRTF set
  * 2 ... show the requested and the actually used HRTF directions
  * 4 ... show the actually used trajectory 
  * 8 ... show the audio block size in samples
  The flags can be logically combined, i.e., 15 shows everything. 
* IgnoreListenerOrientation: if non-zero, the listener orientation is ignored when determining the HRTF direction, i.e., only the source position is considered. 




**Enjoy!**


Acknowledgements:
-----------------

* Christian Hoene for the great libmysofa library for easy reading SOFA files: https://github.com/hoene/libmysofa


History:
--------

Version 1.6.0
=============
* Shows the version of libmysofa
* Closes the log file after loading, re-opens on unloading
* Shows the amplification applied during loudness normalization
* To disable the gain normalization on HRTF-set load, create a file ".normalization_disabled" in the root directory of this project

Version 1.5.0
=============
* Loads up to 50 HRTF sets
* Debug version shows the information of the loaded HRTF sets only

Version 1.4.2
=============
* Automatic gain normalization on HRTF-set load re-added (was removed in v1.4.0)

Version 1.4.1
=============
* Debug 8 added to show the audio block size
* Bug fix in calculating the trajectory

Version 1.4.0
=============
* IgnoreListenerOrientation: added
* Debug: interpretation changed to be able to use as a combination of debug flags
* Compiler warnings removed

Version 1.3.0
=============
* Updated by miho to current version of libmysofa (09.2021), compiled in Windows 10, Visual Studio 2019, 64bit
	updated files: mysofa.lib, mysofa.pdb, mysofa.h
* Minor updates in Plugin_Spatializer.cpp
* SOFAlizer log file info improved

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