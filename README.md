# SOFA-Spatializer-for-Unity
SOFA Native Spatializer Plugin for Unity

Install:
--------

Download libmysofa from https://github.com/hoene/libmysofa (Clone, Rebuild mysofa). The directory "libmysofa" must be at the same level as the directory "SOFA-Spatializer-for-Unity".

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
The plugin uses the file "hrtfs.sofa" in the root directory of this project. This file will be loaded on Unity start. To change to another HRTF, close Unity, overwrite the existing file hrtfs.sofa and reload Unity. Enjoy!