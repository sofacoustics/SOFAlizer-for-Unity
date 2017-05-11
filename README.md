# SOFA-Spatializer-for-Unity
SOFA Native Spatialzer Plugin for Unity

Install:
--------

Download libmysofa from https://github.com/hoene/libmysofa (Clone and Rebuild). The directory "libmysofa" must be at the same level as the directory "SOFA-Spatializer-for-Unity".

Compile SOFAlizer DLL (NativeCode\VisualStudio\AudioPluginSOFA.sln --> Rebuild)

Start the scene in Unity (Assets\SampleScene.unity --> Run)

HRTFs:
------
The plugin uses the file "hrtfs.sofa" in the root directory of this project. This file will be loaded on Unity start. To change to another HRTF, close Unity, overwrite the existing file hrtfs.sofa and reload Unity. Enjoy!