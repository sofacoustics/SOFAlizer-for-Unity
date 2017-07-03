// SOFAlizer - (c) Acoustics Research Institute

// Authors: Piotr Majdak, Claudia Jenny
// Based on Plugin_Spatializer from the package "Native Audio Plugin SDK"
// Licensed under EUPL, see the license file
// Uses libmysofa (Symonics GmbH, Christian Hoene)

// Please note that this will only work on Unity 5.2 or higher.

#define SOFALIZER_VERSION "1.0.1" // SOFAlizer version

#include "AudioPluginUtil.h"
#include "mysofa.h"  // include libmysofa by, Copyright (c) 2016-2017, Symonics GmbH, Christian Hoene

extern float reverbmixbuffer[];

namespace Spatializer
{
    enum
    {
        P_AUDIOSRCATTN,
        P_FIXEDVOLUME,
        P_CUSTOMFALLOFF,
		P_SOFASELECTOR,
        P_NUM
    };

    const int HRTFLEN = 256; // currently max size of the HRTFs
	const int MAX_SOFAS = 10; // currently max number of HRTF sets to be loaded on load.

	const float GAINCORRECTION = 2.0f;

	class HRTFData
	{

	public:
		MYSOFA_HRTF *mysofa[MAX_SOFAS];		// stores the SOFA structure
		MYSOFA_LOOKUP *mylookup[MAX_SOFAS];    // for the lookup
		MYSOFA_NEIGHBORHOOD *myneighborhood[MAX_SOFAS];  // for the lookup
		unsigned int mysofaflag[MAX_SOFAS]; // 0: HRTF set not loaded, 1: HRTF loaded and ready for rendering
		UnityComplexNumber *hrtf[MAX_SOFAS];  // All HRTFs in frequency domain, number of elements = M * R * HRTFLEN * 2
		FILE *Log;	// logging to the file on load/unload of HRTFs
#if _DEBUG
		FILE *pConsole;		// for debugging, use fprintf(pConsole, "my string");
#endif

		HRTFData()
		{
#if _DEBUG
			// Allocate a console for debugging. Use fprintf(sharedData.pConsole, string); for printing
			AllocConsole();
			freopen_s(&pConsole, "CONOUT$", "wb", stdout);
			fprintf(pConsole, "SOFAlizer %s: SOFA-based spatializer (c) Piotr Majdak, ARI, OeAW\n", SOFALIZER_VERSION);
#endif
			fopen_s(&Log, "SOFAlizer.log", "a");
			fprintf(Log, "SOFAlizer %s: SOFA-based spatializer (c) Piotr Majdak, ARI, OeAW\n", SOFALIZER_VERSION);

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
        UnityComplexNumber h[HRTFLEN * 2];
        UnityComplexNumber x[HRTFLEN * 2];
        UnityComplexNumber y[HRTFLEN * 2];
        float buffer[HRTFLEN * 2];
    };

    struct EffectData
    {
        float p[P_NUM];
        InstanceChannel ch[2];
    };

	void LoadSOFAs(UnityAudioEffectState* state)
	{
		if (state->samplerate < 8000) return; // terminate if sampling rate below 8 kHz
#if _DEBUG
		fprintf(sharedData.pConsole, "System sampling rate: %u\n", state->samplerate);
#endif
		fprintf(sharedData.Log, "System sampling rate: %u\n", state->samplerate);
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
					    mysofa_getAttribute(sharedData.mysofa[i]->attributes,"DatabaseName"),
						mysofa_getAttribute(sharedData.mysofa[i]->attributes, "Title"),
						mysofa_getAttribute(sharedData.mysofa[i]->attributes, "ListenerShortName")
						);
#endif
					fprintf(sharedData.Log, "%s: %u HRTFs, %u samples @ %u Hz. %s, %s, %s.\n",
						filename, sharedData.mysofa[i]->M, sharedData.mysofa[i]->N,
						(unsigned int)(sharedData.mysofa[i]->DataSamplingRate.values[0]),
						mysofa_getAttribute(sharedData.mysofa[i]->attributes, "DatabaseName"),
						mysofa_getAttribute(sharedData.mysofa[i]->attributes, "Title"),
						mysofa_getAttribute(sharedData.mysofa[i]->attributes, "ListenerShortName")
					);
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
						mysofa_loudness(sharedData.mysofa[i]);

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
				fprintf(sharedData.pConsole, "%s: Can't load file\n", filename);
#endif
				fprintf(sharedData.Log, "%s: Can't load file\n", filename);
			}
		}//iterate through all sofa files

	}

	void UnloadSOFAs()
	{
#if _DEBUG
		fprintf(sharedData.pConsole, "Unloaded!");
#endif
		fprintf(sharedData.Log, "Unloaded! **************************************************************************\n\n");
		fclose(sharedData.Log);

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
		LoadSOFAs(state);
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ReleaseCallback(UnityAudioEffectState* state)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        delete data;
		// Free the memory from all HRTFs 
		UnloadSOFAs();
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
            memcpy(outbuffer, inbuffer, length * outchannels * sizeof(float));
            return UNITY_AUDIODSP_OK; 
        }

        EffectData* data = state->GetEffectData<EffectData>();
		// in data->ch[0|1].h is space for the HRIRs

        static const float kRad2Deg = 180.0f / kPI;

        float* m = state->spatializerdata->listenermatrix;
        float* s = state->spatializerdata->sourcematrix;

        // Currently we ignore source orientation and only use the position
        float px = s[12];
        float py = s[13];
        float pz = s[14];

        float dir_x = m[0] * px + m[4] * py + m[8] * pz + m[12];
        float dir_y = m[1] * px + m[5] * py + m[9] * pz + m[13];
        float dir_z = m[2] * px + m[6] * py + m[10] * pz + m[14];

		unsigned int Selper = (unsigned int)(data->p[P_SOFASELECTOR]);
#if _DEBUG
		fprintf(sharedData.pConsole, "Set: #%d; ", (int)Selper);
#endif
		if(sharedData.mysofaflag[Selper] == 0)
		{
#if _DEBUG
			fprintf(sharedData.pConsole, "(not loaded)\n");
#endif
			// Set filters to zero (mute)
			for (int chidx = 1; chidx >= 0; chidx--)		// copy left-ear to ch=1 and right-ear to ch=0
				for (int n = 0; n < 2 * HRTFLEN; n++)
				{
					data->ch[chidx].h[n].re = 0;
					data->ch[chidx].h[n].im = 0;
				}
		}
		else
		{
#if _DEBUG
			// Calculate the source direction in spherical coordinates
			float azimuth = (fabsf(dir_z) < 0.001f) ? 0.0f : atan2f(dir_x, dir_z);
			if (azimuth < 0.0f)
				azimuth += 2.0f * kPI;
			azimuth = FastClip(azimuth * kRad2Deg, 0.0f, 360.0f);
			float elevation = atan2f(dir_y, sqrtf(dir_x * dir_x + dir_z * dir_z) + 0.001f) * kRad2Deg;
			fprintf(sharedData.pConsole, "Requested: (%d, %d); ", (int)azimuth, (int)elevation);
#endif
			// Calculate the source direction in cartesian coordinates for the look-up
			float t[3];
			t[0] = dir_z; // Z in Unity is X in SOFA
			t[1] = dir_x; // X in Unity is Y in SOFA
			t[2] = dir_y; // Y in Unity is Z in SOFA

			// Get the index to the nearest HRTF direction
			int nearest = mysofa_lookup(sharedData.mylookup[Selper], t);
#if _DEBUG
			fprintf(sharedData.pConsole, "Found: #%d ", nearest); 
			t[0] = sharedData.mysofa[Selper]->SourcePosition.values[3*nearest];
			t[1] = sharedData.mysofa[Selper]->SourcePosition.values[3*nearest+1];
			t[2] = sharedData.mysofa[Selper]->SourcePosition.values[3*nearest+2];
			mysofa_c2s(t);
			fprintf(sharedData.pConsole, "(%5.1f, %5.1f)\n", t[0], t[1]);
#endif

			// Create a pointer to the left-ear HRTF (the right-ear HRTF is right behind the left-ear)
			UnityComplexNumber *IRL;
			IRL = sharedData.hrtf[Selper] + nearest * (2 * HRTFLEN) * 2;             // nearest * N * R

			// Copy the HRTFs for both ears to the data array
			for (int chidx = 1; chidx >= 0; chidx--)		// copy left-ear to ch=1 and right-ear to ch=0
				for (int n = 0; n < 2 * HRTFLEN; n++)
				{
					data->ch[chidx].h[n].re = IRL->re;
					data->ch[chidx].h[n].im = IRL->im;
					IRL++;
				}
		}

        // From the FMOD documentation:
        //   A spread angle of 0 makes the stereo sound mono at the point of the 3D emitter.
        //   A spread angle of 90 makes the left part of the stereo sound place itself at 45 degrees to the left and the right part 45 degrees to the right.
        //   A spread angle of 180 makes the left part of the stero sound place itself at 90 degrees to the left and the right part 90 degrees to the right.
        //   A spread angle of 360 makes the stereo sound mono at the opposite speaker location to where the 3D emitter should be located (by moving the left part 180 degrees left and the right part 180 degrees right). So in this case, behind you when the sound should be in front of you!
        // Note that FMOD performs the spreading and panning in one go. We can't do this here due to the way that impulse-based spatialization works, so we perform the spread calculations on the left/right source signals before they enter the convolution processing.
        // That way we can still use it to control how the source signal downmixing takes place.
        float spread = cosf(state->spatializerdata->spread * kPI / 360.0f);
        float spreadmatrix[2] = { 2.0f - spread, spread };

		float spatialblend = state->spatializerdata->spatialblend;
		float reverbmix = state->spatializerdata->reverbzonemix;
		
		float* reverb = reverbmixbuffer;
        for (unsigned int sampleOffset = 0; sampleOffset < length; sampleOffset += HRTFLEN)
        {
            for (int c = 0; c < 2; c++)
            {
                // stereopan is in the [-1; 1] range, this acts the way fmod does it for stereo
                float stereopan = 1.0f - ((c == 0) ? FastMax(0.0f, state->spatializerdata->stereopan) : FastMax(0.0f, -state->spatializerdata->stereopan));

                InstanceChannel& ch = data->ch[c];

                for (int n = 0; n < HRTFLEN; n++)
                {
                    float left  = inbuffer[n * 2];
                    float right = inbuffer[n * 2 + 1];
                    ch.buffer[n] = ch.buffer[n + HRTFLEN];
                    ch.buffer[n + HRTFLEN] = left * spreadmatrix[c] + right * spreadmatrix[1 - c];
                }

                for (int n = 0; n < HRTFLEN * 2; n++)
                {
                    ch.x[n].re = ch.buffer[n];
                    ch.x[n].im = 0.0f;
                }

                FFT::Forward(ch.x, HRTFLEN * 2, false);

                for (int n = 0; n < HRTFLEN * 2; n++)
                    UnityComplexNumber::Mul<float, float, float>(ch.x[n], ch.h[n], ch.y[n]);

                FFT::Backward(ch.y, HRTFLEN * 2, false);

                for (int n = 0; n < HRTFLEN; n++)
                {
                    float s = inbuffer[n * 2 + c] * stereopan;
                    float y = s + (ch.y[n].re * GAINCORRECTION - s) * spatialblend;
                    outbuffer[n * 2 + c] = y;
                    reverb[n * 2 + c] += y * reverbmix;
                }
            }

            inbuffer += HRTFLEN * 2;
            outbuffer += HRTFLEN * 2;
            reverb += HRTFLEN * 2;
        }

        return UNITY_AUDIODSP_OK;
    }
}
