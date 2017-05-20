# SOFA-Spatializer-for-Unity
SOFA Native Spatializer Plugin for Unity

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
* Start the engine: Run

HRTFs:
------
On start-up of the Unity scene, the plugin loads SOFA files "hrtfX.sofa" with X being from 0 to 9. These files must be located in the root directory of this project. Currently, all SOFA files must provide HRTFs sampled at the same rate. 

In the scene, the parameter "SOFA Selector" selects the index of the processed HRTF set. The selection of the loaded HRTF sets can be done at any time in real-time. If a non-loaded HRTF set is selected, the audio will be muted.

In order to change the loaded HRTFs, close Unity, replace one or more of the SOFA files, and reload Unity. 

'''Enjoy!'''


Acknowledgements:
-----------------

* Christian Hoene for the great mysofa library. 