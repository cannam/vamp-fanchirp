/*
  copyright (C) 2011 I. Irigaray, M. Rocamora

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Remember to use a different guard symbol in each header!
#ifndef _FCHTRANSFORMF0GRAM_H_
#define _FCHTRANSFORMF0GRAM_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <vamp-sdk/Plugin.h>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
using std::string;

class FChTransformF0gram : public Vamp::Plugin {
public:
    FChTransformF0gram(float inputSampleRate);
    virtual ~FChTransformF0gram();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    string getCopyright() const;
    int getPluginVersion() const;

    InputDomain getInputDomain() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;
    size_t getPreferredStepSize() const;
    size_t getPreferredBlockSize() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:

    string m_currentProgram;
    int m_stepSize;
    int m_blockSize;
    float m_fs; // input sampling rate (inputSampleRate)

    // plugin-specific data and methods go here

    // =============  WARPING PARAMETERS  =============

    double m_fmax; // maximum frequency of interest (Hz)
    int m_nfft; // number of fft points (controls zero-padding)
    int m_hop; // hop in samples in the upsampled signal
    int m_num_f0s; // number of f0 values in F0gram grid
    //vector<float> m_f0s;    // vector of f0 values
    double *m_f0s; // vector of f0 values

    typedef struct {
        int nsamps_twarp; // number of samples of the warped signal frame
        double alpha_max; // maximum value of normalized frequency deviation (alpha)
        int num_warps; // number of warpings
        int fact_over_samp; // oversampling factor
        int alpha_dist; // distribution of alpha values, 'lin' or 'log' (0 - 1)
    } warping_parameters;

    warping_parameters m_warp_params;

    // =============   F0-GRAM PARAMETERS  =============

    typedef struct {
        double f0min; // minimun fundamental frequency
        int num_octs; // number of octaves
        int num_f0s_per_oct; // number of f0s per octave
        int num_f0_hyps; // number of f0 hypotesis to extract
        bool prefer; // whether to use a f0 preference guassian function
        int prefer_mean; // mean of f0 preference function (MIDI number for C4)
        int prefer_stdev; // stdev of f0 preference function (stdev in MIDI numbers)
    } f0_parameters;

    f0_parameters m_f0_params;
    bool m_f0gram_mode;

    // ======== GATHERED LOG SPECTRUM PARAMETERS =======

    typedef struct {
        bool HP_logS; //high-pass logS
        int att_subharms; // whether to attenuate subharmonics
        // model parameter variables (default values)
        double median_poly_coefs[3];
        double sigma_poly_coefs[3];
    } glogs_parameters;

    glogs_parameters m_glogs_params;

    // =============  WARPING DESIGN  =============

    typedef struct {
        double fs_orig; // sampling frequency after oversampling
        double fs_warp; // sampling frequency of warped signal
        double *chirp_rates; // chirp rates
        int nsamps_torig; // number of samples of the original signal frame
        int *pos_int; // index of previous sample to do the warping by interpolation efficiently
        double *pos_frac; // fractional value to do the warping by interpolation efficiently
    } warping_design;

    warping_design m_warpings;
    // LPFWindow
    double *mp_LPFWindow;
    double *LPF_time;
    fftw_complex *LPF_frequency;
    fftw_plan plan_backward_LPF;
    fftw_plan plan_forward_LPF;
    // timeWindow
    double *m_timeWindow;
    // Warpings
    double *x_warping;
    // Hanning window
    double *mp_HanningWindow;
    // FChT plan & transformed data structs
    double *m_absFanChirpTransform;
    fftw_complex *m_auxFanChirpTransform;
    fftw_plan plan_forward_xwarping;
    // GLogS
    double *m_glogs_f0;
    double *m_glogs;
    int *m_glogs_n;
    int *m_glogs_index;
    int *m_glogs_posint;
    double *m_glogs_posfrac;
    double *m_glogs_interp;
    int m_glogs_harmonic_count;
    int m_glogs_num_f0s;
    int m_glogs_init_f0s;
    int *m_glogs_third_harmonic_posint;
    double *m_glogs_third_harmonic_posfrac;
    double *m_glogs_third_harmonic;
    int *m_glogs_fifth_harmonic_posint;
    double *m_glogs_fifth_harmonic_posfrac;
    double *m_glogs_fifth_harmonic;
    double *m_glogs_f0_preference_weights;
    double *m_glogs_median_correction;
    double *m_glogs_sigma_correction;
    double *m_glogs_hf_smoothing_window;
    // auxiliar methods
    void design_GLogS();
    void design_FChT();
    void define_warps_linear_chirps(double *, double *);
    void design_warps(double *, double *, double *);
    void design_LPF();
    void clean_LPF();
    void apply_LPF();
    void design_FFT();
    void design_time_window();

    // FFT variables
    fftw_complex *in, *out;
    //TODO verificar que el tipo de datos de in_window es del tipo double, era del tipo float.
    double *in_window;
    fftw_plan planFFT;
};


#endif
