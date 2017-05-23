# SOFAlizer for Unity

SOFAlizer is a SOFA-based audio spatializer for Unity. It is a Unity [native plugin](https://docs.unity3d.com/Manual/NativePlugins.html). It loads sets of spatial filters so-called [HRTFs](https://en.wikipedia.org/wiki/Head-related_transfer_function) stored in [SOFA](https://www.sofaconventions.org/) files and renders spatial audio based on a selected HRTF set. The output is a binaural signal intended to be listened to via headphones. 


Installation:
-------------

Compile libmysofa:
* Download libmysofa from https://github.com/hoene/libmysofa (Clone in git)
* The directory "libmysofa" must be at the same level as the directory "SOFA-Spatializer-for-Unity".
* Start libmysofa.sln in libmysofa in Visual Studio
* Rebuild mysofa (Release, x64)

Compile SOFAlizer DLL:
* Start NativeCode\VisualStudio\AudioPluginSOFA.sln
* Select x64, Release
* Rebuild AudioPluginSOFA 

Start the scene:
* In Unity start the scene: Assets\SampleScene.unity
* Select SOFAlizer as spatializer plug in: Edit --> Project Settings --> Audio --> Spatializer Plugin --> SOFA Spatializer
* Copy some SOFA files (e.g., from https://www.sofaconventions.org/mediawiki/index.php/Files) to the root directory. Rename them hrtfX.sofa with X from 0 to 9. 
* Start the engine: Run

Usage:
------
On "run" of the Unity scene, SOFAlizer:
* loads SOFA files "hrtfX.sofa" with X being from 0 to 9. These files must be located in the root directory of this project. 
* resamples each HRTF set to the sampling rate provided in Audio --> Project Settings --> Audio --> System Sample Rate. This may take some time (can be avoided by using HRTFs sampled at the same rate as that from the Unity project). 
* normalizes each HRTF set relative to the HRTF for the frontal position. 
* crops the HRTFs to the length of 256 samples. 

In the scene, the parameter "SOFA Selector" selects the index of the processed HRTF set. The selection of the loaded HRTF sets can be done at any time in real-time. If a non-loaded HRTF set is selected, the audio will be muted.

In order to change the loaded HRTFs, stop the scene rendering, replace one or more of the SOFA files, and re-run the scene. 

'''Enjoy!'''


Acknowledgements:
-----------------

* Christian Hoene for the great mysofa library. 