// SOFAlizer - (c) Acoustics Research Institute

// Authors: Piotr Majdak, Claudia Jenny, Michael Mihocic
// Based on Plugin_Spatializer from the package "Native Audio Plugin SDK"
// Licensed under EUPL, see the license file
// Uses libmysofa (Symonics GmbH, Christian Hoene)
// Variable _DEBUG is true when compiler is set to Debug mode

// Please note that this will only work on Unity 5.2 or higher.

#define SOFALIZER_VERSION "1.6.0" // SOFAlizer version

#include "AudioPluginUtil.h"
#include "mysofa.h"  // include libmysofa by, Copyright (c) 2016-2021, Symonics GmbH, Christian Hoene
#include <ctime>

#pragma warning(disable : 4996) // disable warnings about using unsafe functions
#define fequals(a, b) (fabs(a - b) < 0.00001)

extern float reverbmixbuffer[];

namespace Spatializer
{
    enum
    {
        P_AUDIOSRCATTN,
        P_FIXEDVOLUME,
        P_CUSTOMFALLOFF,
		P_SOFASELECTOR,
		P_DEBUGLEVEL,
		P_INGORELISTENERORIENTATION,
		P_NUM
    };

    const int HRTFLEN = 128; // currently max size of the HRTFs
	const int MAX_SOFAS = 50; // currently max number of HRTF sets to be loaded on load.

	//const float GAINCORRECTION = 2.0f;

	static int instanceCounter = 0;

	class HRTFData
	{

	public:
		MYSOFA_HRTF* mysofa[MAX_SOFAS]{};		// stores the SOFA structure
		MYSOFA_LOOKUP* mylookup[MAX_SOFAS]{};    // for the lookup
		MYSOFA_NEIGHBORHOOD* myneighborhood[MAX_SOFAS]{};  // for the lookup
		unsigned int mysofaflag[MAX_SOFAS]{}; // 0: HRTF set not loaded, 1: HRTF loaded and ready for rendering
		UnityComplexNumber* hrtf[MAX_SOFAS]{};  // All HRTFs in frequency domain, number of elements = M * R * HRTFLEN * 2
		FILE* Log{};	// logging to the file on load/unload of HRTFs
		FILE* pConsole{};		// for debugging, use fprintf(pConsole, "my string");

		HRTFData()
		{
			for (unsigned int i = 0; i < MAX_SOFAS; i++)
			{
				// stamp all files as not loaded and not usable
				mysofaflag[i] = 0;
			}
		}

	}; // end of class HRTFData

    static HRTFData sharedData;

    struct InstanceChannel
    {
		UnityComplexNumber x[HRTFLEN * 2] = { 0 };
        UnityComplexNumber y[HRTFLEN * 2];
		float buffer[HRTFLEN * 2] = { 0 };
	};

    struct EffectData
    {
		float p[P_NUM]{};
		InstanceChannel ch[2]{};
		float oldsourcepos[3]{};
    };

	void LoadSOFAs(UnityAudioEffectState* state)
	{
		// time
		char cur_time[128];
		time_t      t;
		struct tm* ptm;
		t = time(NULL);
		ptm = localtime(&t);
		strftime(cur_time, 128, "%Y-%m-%d %H:%M:%S", ptm);
		int MySOFAmajor, MySOFAminor, MySOFApatch;
		int normalization = 1; 

		// Check if normalization disabled
		FILE* fn_normalize=NULL;
		fopen_s(&fn_normalize, ".normalization_disabled", "r");
		if (fn_normalize != NULL) fclose(fn_normalize);

		// Allocate a console for debugging. Use fprintf(sharedData.pConsole, string); for printing
#if _DEBUG
		AllocConsole();
		freopen_s(&(sharedData.pConsole), "CONOUT$", "wb", stdout);
		fprintf(sharedData.pConsole, "SOFAlizer %s: SOFA-based spatializer (c) Piotr Majdak & Michael Mihocic, ARI, OeAW\n", SOFALIZER_VERSION);
		mysofa_getversion(&MySOFAmajor, &MySOFAminor, &MySOFApatch);
		fprintf(sharedData.pConsole, "libmysofa version:%u.%u.%u\n", MySOFAmajor, MySOFAminor, MySOFApatch);
		fprintf(sharedData.pConsole, "%s\n", cur_time);
		fprintf(sharedData.pConsole, "System sampling rate: %u, Buffer size:%u samples\n", state->samplerate, state->dspbuffersize);
#endif
		fopen_s(&(sharedData.Log), "SOFAlizer.log", "a");
		fprintf(sharedData.Log, "SOFAlizer %s: SOFA-based spatializer (c) Piotr Majdak & Michael Mihocic, ARI, OeAW\n", SOFALIZER_VERSION);
		mysofa_getversion(&MySOFAmajor, &MySOFAminor, &MySOFApatch);
		fprintf(sharedData.pConsole, "libmysofa version:%u.%u.%u\n", MySOFAmajor, MySOFAminor, MySOFApatch);
		fprintf(sharedData.Log, "%s\n", cur_time);
		fprintf(sharedData.Log, "System sampling rate: %u Hz, Buffer size:%u samples\n", state->samplerate, state->dspbuffersize);

		if (state->samplerate < 8000) return; // terminate if sampling rate below 8 kHz

		// Iterate through the SOFA files
		for (unsigned int i = 0; i < MAX_SOFAS; i++)
		{
			// first, assume this file as not loaded and not usable
			sharedData.mysofaflag[i] = 0;
			// Open the SOFA file
			int err;
			char filename[50];
			sprintf_s(filename, sizeof(filename), "hrtf%u.sofa", i);
			sharedData.mysofa[i] = mysofa_load(filename, &err);
			if (sharedData.mysofa[i])
			{
				err = mysofa_check(sharedData.mysofa[i]); 	// Check the loaded structure
				if (err == MYSOFA_OK)
				{
#if _DEBUG
					fprintf(sharedData.pConsole, "%s: %u HRTFs, %u samples @ %u Hz. %s, %s, %s.\n",
							filename, sharedData.mysofa[i]->M, sharedData.mysofa[i]->N,
							(unsigned int)(sharedData.mysofa[i]->DataSamplingRate.values[0]),
							mysofa_getAttribute(sharedData.mysofa[i]->attributes, "DatabaseName"),
							mysofa_getAttribute(sharedData.mysofa[i]->attributes, "Title"),
							mysofa_getAttribute(sharedData.mysofa[i]->attributes, "ListenerShortName"));
#endif
					fprintf(sharedData.Log, "%s: %u HRTFs, %u samples @ %u Hz. %s, %s, %s.\n",
							filename, sharedData.mysofa[i]->M, sharedData.mysofa[i]->N,
							(unsigned int)(sharedData.mysofa[i]->DataSamplingRate.values[0]),
							mysofa_getAttribute(sharedData.mysofa[i]->attributes, "DatabaseName"),
							mysofa_getAttribute(sharedData.mysofa[i]->attributes, "Title"),
							mysofa_getAttribute(sharedData.mysofa[i]->attributes, "ListenerShortName"));
					// Convert to cartesian, initialize the look up
					mysofa_tocartesian(sharedData.mysofa[i]);
					sharedData.mylookup[i] = mysofa_lookup_init(sharedData.mysofa[i]);
					if (sharedData.mylookup[i])
					{
						// initialize lookup table
						sharedData.myneighborhood[i] = mysofa_neighborhood_init(sharedData.mysofa[i], sharedData.mylookup[i]);
						// resample if required
						if (state->samplerate != sharedData.mysofa[i]->DataSamplingRate.values[0])
						{
							err = mysofa_resample(sharedData.mysofa[i], (float)state->samplerate);
#if _DEBUG
							fprintf(sharedData.pConsole, "--> resampled to %5.0f Hz, new IR length: %u samples\n", sharedData.mysofa[i]->DataSamplingRate.values[0], sharedData.mysofa[i]->N);
#endif
							fprintf(sharedData.Log, "--> resampled to %5.0f Hz, new IR length: %u samples\n", sharedData.mysofa[i]->DataSamplingRate.values[0], sharedData.mysofa[i]->N);
						}

						// Scale HRTFs to have a normalized amplitude relative to the frontal position
						float gain = 1.0f;
						if (fn_normalize == NULL) gain=mysofa_loudness(sharedData.mysofa[i]);
#if _DEBUG
						if ((gain > 1.01) || (gain < 0.99)) fprintf(sharedData.pConsole, "--> amplified by %6.2f\n", gain);
#endif

						// Determine the final length of each IR
						unsigned int length; // length of the IRs to be copied 
						if (sharedData.mysofa[i]->N > HRTFLEN)
						{
							length = HRTFLEN; // too long: only the HRTFLEN part of each IR
#if _DEBUG
							fprintf(sharedData.pConsole, "--> cropped to %u samples\n", HRTFLEN);
#endif
							fprintf(sharedData.Log, "--> cropped to %u samples\n", HRTFLEN);
						}
						else
						{
							length = sharedData.mysofa[i]->N; // OK, consider the full IR
						}

						// Transform the HRIRs to complex-valued spectral filters and copy to the HRTF array called myhrtf
						sharedData.hrtf[i] = (UnityComplexNumber *)malloc(sizeof(UnityComplexNumber) * sharedData.mysofa[i]->M * sharedData.mysofa[i]->R * 2 * HRTFLEN);
						UnityComplexNumber h[HRTFLEN * 2]; // for temporary impulse response and spectrum
						float *hrir = sharedData.mysofa[i]->DataIR.values; // temporary pointer to the source array
						UnityComplexNumber *hrtf = sharedData.hrtf[i]; // temporary pointer to the destination array						
						for (unsigned int a = 0; a < sharedData.mysofa[i]->M * sharedData.mysofa[i]->R; a++)
						{
							// copy from source array
							memset(h, 0, sizeof(h)); // clear buffer 
							for (unsigned int n = 0; n < length; n++)
								h[n + HRTFLEN].re = hrir[a*sharedData.mysofa[i]->N + n];
							// FFT
							FFT::Forward(h, HRTFLEN * 2, false);
							// Copy to destination array
							for (int n = 0; n < HRTFLEN * 2; n++)
								*(hrtf++) = h[n];
						}
						sharedData.mysofaflag[i] = 1; // set this file as usable
					}
					else
					{
#if _DEBUG
						fprintf(sharedData.pConsole, "%s: Look-up init failed!\n", filename);
#endif
						fprintf(sharedData.Log, "%s: Look-up init failed!\n", filename);
					}
				}
				else
				{
#if _DEBUG
					fprintf(sharedData.pConsole, "%s: Check failed!\n", filename);
#endif
					fprintf(sharedData.Log, "%s: Check failed!\n", filename);
				}
			}
			else
			{
#if _DEBUG
				//fprintf(sharedData.pConsole, "%s: Can't load, error code: %d (%s)\n", filename, err, strerror(err));
#endif
				//fprintf(sharedData.Log, "%s: Can't load, error code: %d (%s)\n", filename, err, strerror(err));
			}
		}//iterate through all sofa files

		fclose(sharedData.Log);
	}

	void UnloadSOFAs(void)
	{
#if _DEBUG
		fprintf(sharedData.pConsole, "Unloaded!\n\n");
#endif
		fopen_s(&(sharedData.Log), "SOFAlizer.log", "a");
		if (sharedData.Log != 0)
		{
			fprintf(sharedData.Log, "Unloaded! *******************************\n\n");
			fclose(sharedData.Log);
		}

		for (unsigned int i = 0; i < MAX_SOFAS; i++)
		{
			if(sharedData.mysofaflag[i])
			{
				mysofa_free(sharedData.mysofa[i]);
				free(sharedData.hrtf[i]);
				mysofa_lookup_free(sharedData.mylookup[i]);
				mysofa_neighborhood_free(sharedData.myneighborhood[i]);
			}
		}
	}

    inline bool IsHostCompatible(UnityAudioEffectState* state)
    {
        // Somewhat convoluted error checking here because hostapiversion is only supported from SDK version 1.03 (i.e. Unity 5.2) and onwards.
        return
            state->structsize >= sizeof(UnityAudioEffectState) &&
            state->hostapiversion >= UNITY_AUDIO_PLUGIN_API_VERSION;
    }

    int InternalRegisterEffectDefinition(UnityAudioEffectDefinition& definition)
    {
        int numparams = P_NUM;
        definition.paramdefs = new UnityAudioParameterDefinition[numparams];
        RegisterParameter(definition, "AudioSrc Attn", "", 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, P_AUDIOSRCATTN, "AudioSource distance attenuation");
        RegisterParameter(definition, "Fixed Volume", "", 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, P_FIXEDVOLUME, "Fixed volume amount");
        RegisterParameter(definition, "Custom Falloff", "", 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, P_CUSTOMFALLOFF, "Custom volume falloff amount (logarithmic)");
		RegisterParameter(definition, "SOFA Selector", "", 0.0f, MAX_SOFAS-1, 0.0f, 1.0f, 1.0f, P_SOFASELECTOR, "HRTF Selector");
		RegisterParameter(definition, "Debug Console", "", 0.0f, 128.0f, 0.0f, 1.0f, 1.0f, P_DEBUGLEVEL, "Level of debugging information to the console");
		RegisterParameter(definition, "Ign List Orien", "", 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, P_INGORELISTENERORIENTATION, "0..Use actual listener orientation; 1..Ignore actual listener orientation");
		definition.flags |= UnityAudioEffectDefinitionFlags_IsSpatializer;
        return numparams;
    }

    static UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK DistanceAttenuationCallback(UnityAudioEffectState* state, float distanceIn, float attenuationIn, float* attenuationOut)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        *attenuationOut =
            data->p[P_AUDIOSRCATTN] * attenuationIn +
            data->p[P_FIXEDVOLUME] +
            data->p[P_CUSTOMFALLOFF] * (1.0f / FastMax(1.0f, distanceIn));
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK CreateCallback(UnityAudioEffectState* state)
    {
		EffectData* effectdata = new EffectData;
        memset(effectdata, 0, sizeof(EffectData));
        state->effectdata = effectdata;
        if (IsHostCompatible(state))
            state->spatializerdata->distanceattenuationcallback = DistanceAttenuationCallback;
        InitParametersFromDefinitions(InternalRegisterEffectDefinition, effectdata->p);
		// Load all SOFA files
		if (0 == instanceCounter) {
			LoadSOFAs(state);
		}
		
		instanceCounter++;

        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ReleaseCallback(UnityAudioEffectState* state)
    {		
		EffectData* data = state->GetEffectData<EffectData>();
		delete data;

		instanceCounter--;
		
		// Free the memory from all HRTFs 
		if (0 == instanceCounter) {
			UnloadSOFAs();
		}
		
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK SetFloatParameterCallback(UnityAudioEffectState* state, int index, float value)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        if (index >= P_NUM)
            return UNITY_AUDIODSP_ERR_UNSUPPORTED;
        data->p[index] = value;
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK GetFloatParameterCallback(UnityAudioEffectState* state, int index, float* value, char *valuestr)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        if (index >= P_NUM)
            return UNITY_AUDIODSP_ERR_UNSUPPORTED;
        if (value != NULL)
            *value = data->p[index];
        if (valuestr != NULL)
            valuestr[0] = 0;
        return UNITY_AUDIODSP_OK;
    }

    int UNITY_AUDIODSP_CALLBACK GetFloatBufferCallback(UnityAudioEffectState* state, const char* name, float* buffer, int numsamples)
    {
        return UNITY_AUDIODSP_OK;
    }


    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ProcessCallback(UnityAudioEffectState* state, float* inbuffer, float* outbuffer, unsigned int length, int inchannels, int outchannels)
    {
        // Check that I/O formats are right and that the host API supports this feature
        if (inchannels != 2 || outchannels != 2 ||
            !IsHostCompatible(state) || state->spatializerdata == NULL)
        {
            memcpy(outbuffer, inbuffer, (size_t)length * outchannels * sizeof(float));
            return UNITY_AUDIODSP_OK;
        }

        static const float kRad2Deg = 180.0f / kPI;

		// Get plugin parameters
		EffectData* data = state->GetEffectData<EffectData>();
		unsigned int Selper = (unsigned int)(data->p[P_SOFASELECTOR]);
		unsigned int debug = (unsigned int)(data->p[P_DEBUGLEVEL]);
		unsigned int IgnoreListenerOrientation = (unsigned int)(data->p[P_INGORELISTENERORIENTATION]);

		// Get Listener and Source spatial data
		float* m = state->spatializerdata->listenermatrix;
		float* s = state->spatializerdata->sourcematrix;

		// Calculate the direction vector of the required HRTF
		float dir_x, dir_y, dir_z; 
		if (IgnoreListenerOrientation)
		{	// Ignore the listener orientation (and currently the position as well)
			dir_x = s[12];
			dir_y = s[13];
			dir_z = s[14];
		}
		else
		{	// Consider the listener orientation, ignoring the source orientation (=directivity)
			dir_x = m[0] * s[12] + m[4] * s[13] + m[8] * s[14] + m[12];
			dir_y = m[1] * s[12] + m[5] * s[13] + m[9] * s[14] + m[13];
			dir_z = m[2] * s[12] + m[6] * s[13] + m[10] * s[14] + m[14];
		}

		if(sharedData.mysofaflag[Selper] == 0)
		{
#if _DEBUG
			if (debug & 1) fprintf(sharedData.pConsole, "Active HRTF Set: #%d (not loaded)\n", (int)Selper);
#endif
			// Set filters to zero (mute)
			memset(outbuffer, 0, (size_t) length * outchannels * sizeof(float));
			return UNITY_AUDIODSP_OK;
		}
#if _DEBUG		
		if (debug & 1) fprintf(sharedData.pConsole, "Active HRTF Set: #%d loaded. ", (int)Selper);

		// Print the source direction in spherical coordinates
		if (debug & 2) {
			float azimuth = (fabsf(dir_z) < 0.001f) ? 0.0f : atan2f(dir_x, dir_z);
			if (azimuth < 0.0f)
				azimuth += 2.0f * kPI;
			azimuth = FastClip(azimuth * kRad2Deg, 0.0f, 360.0f);
			float elevation = atan2f(dir_y, sqrtf(dir_x * dir_x + dir_z * dir_z) + 0.001f) * kRad2Deg;
			if (IgnoreListenerOrientation)
				fprintf(sharedData.pConsole, "Source Direction (%d, %d) ", (int)azimuth, (int)elevation);
			else
				fprintf(sharedData.pConsole, "Source-Listener Direction (%d, %d) ", (int)azimuth, (int)elevation);
		}
#endif
		// Calculate the source direction in cartesian coordinates
		float t[3], o[3];

		o[0] = data->oldsourcepos[0];  // get previously saved position
		o[1] = data->oldsourcepos[1];
		o[2] = data->oldsourcepos[2];

		t[0] = dir_z; // Z in Unity is X in SOFA
		t[1] = -dir_x; // X in Unity is -Y in SOFA
		t[2] = dir_y; // Y in Unity is Z in SOFA

		data->oldsourcepos[0] = t[0]; // save current position as old
		data->oldsourcepos[1] = t[1];
		data->oldsourcepos[2] = t[2];

		// Get the indicies to the nearest HRTF directions, interpolated between the old position and the new one
		int nearest[16]{}, nidx = 0;
		float pos[3], nmax = float(length) / (HRTFLEN * 2);
		for (unsigned int sampleOffset = 0; sampleOffset < length; sampleOffset += HRTFLEN)
		{
			pos[0] = (t[0] - o[0]) / nmax*(nidx + 1) + o[0];
			pos[1] = (t[1] - o[1]) / nmax*(nidx + 1) + o[1];
			pos[2] = (t[2] - o[2]) / nmax*(nidx + 1) + o[2];
			nearest[nidx] = mysofa_lookup(sharedData.mylookup[Selper], pos);
#if _DEBUG
			if (debug & 4) {
				pos[0] = sharedData.mysofa[Selper]->SourcePosition.values[3 * nearest[nidx]];
				pos[1] = sharedData.mysofa[Selper]->SourcePosition.values[3 * nearest[nidx] + 1];
				pos[2] = sharedData.mysofa[Selper]->SourcePosition.values[3 * nearest[nidx] + 2];
				mysofa_c2s(pos);
				fprintf(sharedData.pConsole, "-> (%5.1f, %5.1f)", pos[0], pos[1]);
			}
#endif
			nidx++;
		}

		// Create a pointer to the left-ear HRTF (the right-ear HRTF is right behind the left-ear)
		UnityComplexNumber *IRL;

		float reverbmix = state->spatializerdata->reverbzonemix;
		float* reverb = reverbmixbuffer;
		nidx = 0;
		for (unsigned int sampleOffset = 0; sampleOffset < length; sampleOffset += HRTFLEN)
		{
			for (int c = 0; c < 2; c++)
			{
					// Pointer to the correct HRTF: [nearest * N * R]
				IRL = sharedData.hrtf[Selper] + (size_t) nearest[nidx] * (2 * HRTFLEN) * 2 + 2*(size_t)HRTFLEN*c;
					// Pointer to the channel data
				InstanceChannel& ch = data->ch[c];

				// Implements Overlap and save
					// Save the input
				for (int n = 0; n < HRTFLEN; n++)
				{
					ch.buffer[n] = ch.buffer[n + HRTFLEN];    // save the previous buffer
					ch.buffer[n + HRTFLEN] = inbuffer[n * 2]; // always use left channel
				}
					// FT of the input
				for (int n = 0; n < HRTFLEN * 2; n++)
				{
					ch.x[n].re = ch.buffer[n];
					ch.x[n].im = 0.0f;
				}
				FFT::Forward(ch.x, HRTFLEN * 2, false);
					// spectral multiplication
				for (int n = 0; n < HRTFLEN * 2; n++)
					UnityComplexNumber::Mul<float, float, float>(ch.x[n], IRL[n], ch.y[n]);
					// IF of the result
				FFT::Backward(ch.y, HRTFLEN * 2, false);
					// Copy the output to the buffer, discard samples after HRTFLEN
				for (int n = 0; n < HRTFLEN; n++)
				{
					outbuffer[n * 2 + c] = ch.y[n].re;
					reverb[n * 2 + c] += ch.y[n].re * reverbmix;
				}
			}

			inbuffer += HRTFLEN * 2; // HRTFLEN * 2 channels
			outbuffer += HRTFLEN * 2;
			reverb += HRTFLEN * 2;
			nidx++;
		}

#if _DEBUG
		if (debug & 2) {
			if (nidx <1 ) nidx = 1;
			pos[0] = sharedData.mysofa[Selper]->SourcePosition.values[3 * nearest[nidx-1]];
			pos[1] = sharedData.mysofa[Selper]->SourcePosition.values[3 * nearest[nidx-1] + 1];
			pos[2] = sharedData.mysofa[Selper]->SourcePosition.values[3 * nearest[nidx-1] + 2];
			mysofa_c2s(pos);
			fprintf(sharedData.pConsole, " -> HRTF found for (%5.1f, %5.1f)", pos[0], pos[1]);
		}
		if (debug & 8) {
			fprintf(sharedData.pConsole, "Audio block size: %u",length);
		}
		if (debug > 0) fprintf(sharedData.pConsole, "\n");
#endif
        return UNITY_AUDIODSP_OK;
    }
}
