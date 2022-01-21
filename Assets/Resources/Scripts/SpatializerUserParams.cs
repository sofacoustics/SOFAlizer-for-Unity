// The spatialization API is only supported by the final Unity 5.2 version and newer.
// If you get script compile errors in this file, comment out the line below.
//#define ENABLE_SPATIALIZER_API

using UnityEngine;
using System.Collections;

public class SpatializerUserParams : MonoBehaviour
{
//    #if ENABLE_SPATIALIZER_API
    public bool EnableSpatialization = true;
    public float DistanceAttn = 1.0f;
    public float FixedVolume = 0.0f;
    public float CustomRolloff = 0.0f;
	public float SOFASelector = 0.0f;
    public float DebugConsole = 0.0f;
    public float IgnoreListenerOrientation = 0.0f;
//#endif

    void Start()
    {
    }

    void Update()
    {
        var source = GetComponent<AudioSource>();
//        #if ENABLE_SPATIALIZER_API
        source.SetSpatializerFloat(0, DistanceAttn);
        source.SetSpatializerFloat(1, FixedVolume);
        source.SetSpatializerFloat(2, CustomRolloff);
        source.SetSpatializerFloat(3, SOFASelector);
        source.SetSpatializerFloat(4, DebugConsole);
        source.SetSpatializerFloat(5, IgnoreListenerOrientation);
        source.GetSpatializerFloat(0, out DistanceAttn); // Get back clipped parameters from plugin
		source.GetSpatializerFloat(1, out FixedVolume);
		source.GetSpatializerFloat(2, out CustomRolloff);
		source.GetSpatializerFloat(3, out SOFASelector);        
        source.GetSpatializerFloat(4, out DebugConsole);
        source.GetSpatializerFloat(5, out IgnoreListenerOrientation);
        source.spatialize = EnableSpatialization;
 //       #endif

		if (Input.GetKeyDown (KeyCode.Alpha0)) {
			SOFASelector = 0;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha1)) {
			SOFASelector = 1;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha2)) {
			SOFASelector = 2;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha3)) {
			SOFASelector = 3;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha4)) {
			SOFASelector = 4;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha5)) {
			SOFASelector = 5;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha6)) {
			SOFASelector = 6;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha7)) {
			SOFASelector = 7;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha8)) {
			SOFASelector = 8;
		}
		else if (Input.GetKeyDown (KeyCode.Alpha9)) {
			SOFASelector = 9;
		}
    }
}
