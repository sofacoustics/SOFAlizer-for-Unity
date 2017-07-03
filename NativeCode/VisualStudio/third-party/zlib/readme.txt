This directory stores the zlib libraries from https://zlib.net/ written by Jean-loup Gailly (compression) and Mark Adler (decompression). 
We use version 1.2.11 of 15.1.2017. SOFAlizer requires two files:
* zlibstat.lib
* zlibstat.pdb (only for debugging pursposes)

To recreate these files, 
1) download the zlib-1.2.11 package from https://zlib.net/
2) Go to contrib\vstudio\vc14
3) Open the VS project in VS
4) In the options, change the RuntimeLibrary from "MultiThreadedDLL" to "MultiThreaded"
5) Rebuild
