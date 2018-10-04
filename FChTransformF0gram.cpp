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

#include "FChTransformF0gram.h"
#include "FChTransformUtils.h"
#include <math.h>
#include <float.h>

#include <set>

#include "bqvec/Allocators.h"

using namespace breakfastquay;

#define DEBUG

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

FChTransformF0gram::FChTransformF0gram(ProcessingMode mode,
                                       float inputSampleRate) :
    Plugin(inputSampleRate),
    m_processingMode(mode),
    m_initialised(false),
    m_stepSize(256),
    m_blockSize(8192) {

    nsamp_options.push_back(256);
    nsamp_options.push_back(512);
    nsamp_options.push_back(1024);
    nsamp_options.push_back(2048);
    nsamp_options.push_back(4096);
    nsamp_options.push_back(8192);
    
    m_fs = inputSampleRate;
    // max frequency of interest (Hz)
    m_fmax = 10000.f;
    // warping parameters
    m_warp_params.nsamps_twarp = 2048;
    m_warp_params.alpha_max = 4;
    m_warp_params.num_warps = 21;
    m_warp_params.fact_over_samp = 2;
    m_warp_params.alpha_dist = 0;
    // f0 parameters
    m_f0_params.f0min = 80.0;
    m_f0_params.num_octs = 4;
    m_f0_params.num_f0s_per_oct = 192;
    m_f0_params.prefer = true;
    m_f0_params.prefer_mean = 60;
    m_f0_params.prefer_stdev = 18;
    // glogs parameters
    m_glogs_params.HP_logS = true;
    m_glogs_params.att_subharms = 1;
    // display parameters
    m_f0gram_mode = BestBinOfAllDirections;

    m_glogs_params.median_poly_coefs[0] = -0.000000058551680;
    m_glogs_params.median_poly_coefs[1] = -0.000006945207775;
    m_glogs_params.median_poly_coefs[2] = 0.002357223226588;

    m_glogs_params.sigma_poly_coefs[0] = 0.000000092782308;
    m_glogs_params.sigma_poly_coefs[1] = 0.000057283574898;
    m_glogs_params.sigma_poly_coefs[2] = 0.022199903714288;

    m_num_f0s = 0;
    m_f0s = 0;
}

FChTransformF0gram::~FChTransformF0gram()
{
    if (!m_initialised) {
        return; // nothing was allocated
    }

    deallocate(m_inputBuffer);
    
    deallocate(m_warpings.pos_int);
    deallocate(m_warpings.pos_frac);
    deallocate(m_warpings.chirp_rates);

    clean_LPF();

    deallocate(m_timeWindow);

    deallocate(mp_HanningWindow);
    
    // Warping
    deallocate(x_warping);
    delete fft_xwarping;
    deallocate(m_absFanChirpTransform); 
    deallocate(m_auxFanChirpTransform);

    // design_GLogS
    deallocate(m_glogs_f0);
    deallocate(m_glogs);
    deallocate(m_glogs_n);
    deallocate(m_glogs_index);
    deallocate(m_glogs_posint);
    deallocate(m_glogs_posfrac);
    deallocate(m_glogs_interp);
    deallocate(m_glogs_third_harmonic_posint);
    deallocate(m_glogs_third_harmonic_posfrac);
    deallocate(m_glogs_third_harmonic);
    deallocate(m_glogs_fifth_harmonic_posint);
    deallocate(m_glogs_fifth_harmonic_posfrac);
    deallocate(m_glogs_fifth_harmonic);
    deallocate(m_glogs_f0_preference_weights);
    deallocate(m_glogs_median_correction);
    deallocate(m_glogs_sigma_correction);

    deallocate(m_f0s);
}

string
FChTransformF0gram::getIdentifier() const {
    switch (m_processingMode) {
    case ModeF0Gram: return "fchtransformf0gram";
    case ModeSpectrogram: return "fchtransformspectrogram";
    case ModeRoughSpectrogram: return "fchtransformrough";
    }
    throw std::logic_error("unknown mode");
}

string
FChTransformF0gram::getName() const {
    switch (m_processingMode) {
    case ModeF0Gram: return "Fan Chirp Transform F0gram";
    case ModeSpectrogram: return "Fan Chirp Transform Spectrogram";
    case ModeRoughSpectrogram: return "Fan Chirp Transform Rough Spectrogram";
    }
    throw std::logic_error("unknown mode");
}

string
FChTransformF0gram::getDescription() const {
    switch (m_processingMode) {
    case ModeF0Gram: 
        return "This plug-in produces a representation, called F0gram, which exhibits the salience of the fundamental frequency of the sound sources in the audio file. The computation of the F0gram makes use of the Fan Chirp Transform analysis. It is based on the article \"Fan chirp transform for music representation\"  P. Cancela, E. Lopez, M. Rocamora, International Conference on Digital Audio Effects, 13th. DAFx-10. Graz, Austria - 6-10 Sep 2010.";
    case ModeSpectrogram:
        return "This plug-in produces a spectral representation of the audio using Fan Chirp Transform analysis.";
    case ModeRoughSpectrogram:
        return "This plug-in produces a more approximate spectral representation of the audio using Fan Chirp Transform analysis.";
    }
    throw std::logic_error("unknown mode");
}

string
FChTransformF0gram::getMaker() const {
    // Your name here
    return "Audio Processing Group \n Universidad de la Republica";
}

int
FChTransformF0gram::getPluginVersion() const {
    // Increment this each time you release a version that behaves
    // differently from the previous one
    //
    // 0 - initial version from scratch
    return 1;
}

string
FChTransformF0gram::getCopyright() const {
    // This function is not ideally named.  It does not necessarily
    // need to say who made the plugin -- getMaker does that -- but it
    // should indicate the terms under which it is distributed.  For
    // example, "Copyright (year). All Rights Reserved", or "GPL"
    return "copyright (C) 2011 GPL - Audio Processing Group, UdelaR";
}

FChTransformF0gram::InputDomain
FChTransformF0gram::getInputDomain() const {
    return TimeDomain;
}

size_t FChTransformF0gram::getPreferredBlockSize() const {
    // We do our own accumulating into blocks within process()
    return m_blockSize/2;
}

size_t
FChTransformF0gram::getPreferredStepSize() const {
    return m_stepSize;
}

size_t
FChTransformF0gram::getMinChannelCount() const {
    return 1;
}

size_t
FChTransformF0gram::getMaxChannelCount() const {
    return 1;
}

FChTransformF0gram::ParameterList
FChTransformF0gram::getParameterDescriptors() const {
    ParameterList list;

    // If the plugin has no adjustable parameters, return an empty
    // list here (and there's no need to provide implementations of
    // getParameter and setParameter in that case either).

    // Note that it is your responsibility to make sure the parameters
    // start off having their default values (e.g. in the constructor
    // above).  The host needs to know the default value so it can do
    // things like provide a "reset to default" function, but it will
    // not explicitly set your parameters to their defaults for you if
    // they have not changed in the mean time.

    // =============  WARPING PARAMETERS  =============

    ParameterDescriptor fmax;
    fmax.identifier = "fmax";
    fmax.name = "Maximum frequency";
    fmax.description = "Maximum frequency of interest for the analysis.";
    fmax.unit = "Hz";
    fmax.minValue = 2000;
    fmax.maxValue = 22050;
    fmax.defaultValue = 10000;
    fmax.isQuantized = true;
    fmax.quantizeStep = 1.0;
    list.push_back(fmax);

    ParameterDescriptor nsamp;
    nsamp.identifier = "nsamp_ix";
    nsamp.name = "Number of samples";
    nsamp.description = "Number of samples of the time warped frame";
    nsamp.minValue = 0;
    nsamp.maxValue = nsamp_options.size()-1;
    nsamp.defaultValue = 3;
    nsamp.isQuantized = true;
    nsamp.quantizeStep = 1.0;
    char label[100];
    for (int i = 0; i < int(nsamp_options.size()); ++i) {
        sprintf(label, "%d", nsamp_options[i]);
        nsamp.valueNames.push_back(label);
    }
    nsamp.isQuantized = true;
    nsamp.quantizeStep = 1.0;
    list.push_back(nsamp);

    ParameterDescriptor alpha_max;
    alpha_max.identifier = "alpha_max";
    alpha_max.name = "Maximum alpha value";
    alpha_max.description = "Maximum value for the alpha parameter of the transform.";
    alpha_max.unit = "Hz/s";
    alpha_max.minValue = -10;
    alpha_max.maxValue = 10;
    alpha_max.defaultValue = 5;
    alpha_max.isQuantized = true;
    alpha_max.quantizeStep = 1.0;
    list.push_back(alpha_max);

    ParameterDescriptor num_warps;
    num_warps.identifier = "num_warps";
    num_warps.name = "Number of warpings";
    num_warps.description = "Number of different warpings in the specified range (must be odd).";
    num_warps.unit = "";
    num_warps.minValue = 1;
    num_warps.maxValue = 101;
    num_warps.defaultValue = 21;
    num_warps.isQuantized = true;
    num_warps.quantizeStep = 2.0;
    list.push_back(num_warps);

    ParameterDescriptor alpha_dist;
    alpha_dist.identifier = "alpha_dist";
    alpha_dist.name = "alpha distribution";
    alpha_dist.description = "Type of distribution of alpha values (linear or log).";
    alpha_dist.unit = "";
    alpha_dist.minValue = 0;
    alpha_dist.maxValue = 1;
    alpha_dist.defaultValue = 1;
    alpha_dist.isQuantized = true;
    alpha_dist.quantizeStep = 1.0;
    // lin (0), log (1)
    alpha_dist.valueNames.push_back("lin");
    alpha_dist.valueNames.push_back("log");
    list.push_back(alpha_dist);

    // =============   F0-GRAM PARAMETERS  =============

    ParameterDescriptor f0min;
    f0min.identifier = "f0min";
    f0min.name = "min f0";
    f0min.description = "Minimum fundamental frequency (f0) value.";
    f0min.unit = "Hz";
    f0min.minValue = 1;
    f0min.maxValue = 500;
    f0min.defaultValue = 80;
    f0min.isQuantized = true;
    f0min.quantizeStep = 1.0;
    list.push_back(f0min);

    ParameterDescriptor num_octs;
    num_octs.identifier = "num_octs";
    num_octs.name = "number of octaves";
    num_octs.description = "Number of octaves for F0gram computation.";
    num_octs.unit = "";
    num_octs.minValue = 1;
    num_octs.maxValue = 10;
    num_octs.defaultValue = 4;
    num_octs.isQuantized = true;
    num_octs.quantizeStep = 1.0;
    list.push_back(num_octs);

    ParameterDescriptor f0s_per_oct;
    f0s_per_oct.identifier = "f0s_per_oct";
    f0s_per_oct.name = "f0 values per octave";
    f0s_per_oct.description = "Number of f0 values per octave.";
    f0s_per_oct.unit = "";
    f0s_per_oct.minValue = 12;
    f0s_per_oct.maxValue = 768;
    f0s_per_oct.defaultValue = 192;
    f0s_per_oct.isQuantized = true;
    f0s_per_oct.quantizeStep = 1.0;
    list.push_back(f0s_per_oct);

    ParameterDescriptor f0_prefer_fun;
    f0_prefer_fun.identifier = "f0_prefer_fun";
    f0_prefer_fun.name = "Use f0 weighting";
    f0_prefer_fun.description = "Whether to use a f0 weighting function to prefer frequencies nearer a mean value.";
    f0_prefer_fun.unit = "";
    f0_prefer_fun.minValue = 0;
    f0_prefer_fun.maxValue = 1;
    f0_prefer_fun.defaultValue = 1;
    f0_prefer_fun.isQuantized = true;
    f0_prefer_fun.quantizeStep = 1.0;
    list.push_back(f0_prefer_fun);

    ParameterDescriptor f0_prefer_mean;
    f0_prefer_mean.identifier = "f0_prefer_mean";
    f0_prefer_mean.name = "Mean pitch for f0 weighting";
    f0_prefer_mean.description = "Mean value for f0 weighting function (MIDI number).";
    f0_prefer_mean.unit = "";
    f0_prefer_mean.minValue = 1;
    f0_prefer_mean.maxValue = 127;
    f0_prefer_mean.defaultValue = 60;
    f0_prefer_mean.isQuantized = true;
    f0_prefer_mean.quantizeStep = 1.0;
    list.push_back(f0_prefer_mean);

    ParameterDescriptor f0_prefer_stdev;
    f0_prefer_stdev.identifier = "f0_prefer_stdev";
    f0_prefer_stdev.name = "Stdev for f0 weighting";
    f0_prefer_stdev.description = "Standard deviation for f0 weighting function (MIDI number).";
    f0_prefer_stdev.unit = "";
    f0_prefer_stdev.minValue = 1;
    f0_prefer_stdev.maxValue = 127;
    f0_prefer_stdev.defaultValue = 18;
    f0_prefer_stdev.isQuantized = true;
    f0_prefer_stdev.quantizeStep = 1.0;
    list.push_back(f0_prefer_stdev);

    ParameterDescriptor f0gram_mode;
    f0gram_mode.identifier = "f0gram_mode";
    f0gram_mode.name = "display mode of f0gram";
    f0gram_mode.description = "Display all bins of the best direction, or the best bin for each direction.";
    f0gram_mode.unit = "";
    f0gram_mode.minValue = 0;
    f0gram_mode.maxValue = 1;
    f0gram_mode.defaultValue = 1;
    f0gram_mode.isQuantized = true;
    f0gram_mode.quantizeStep = 1.0;
    list.push_back(f0gram_mode);

    return list;
}

float
FChTransformF0gram::getParameter(string identifier) const {

    if (identifier == "fmax") {
        return m_fmax;
    } else if (identifier == "nsamp_ix") {
        for (int i = 0; i < int(nsamp_options.size()); ++i) {
            if (m_warp_params.nsamps_twarp == nsamp_options[i]) {
                return i;
            }
        }
        throw std::logic_error("internal error: nsamps_twarp not in nsamp_options");
    } else if (identifier == "alpha_max") {
        return m_warp_params.alpha_max;
    } else if (identifier == "num_warps") {
        return m_warp_params.num_warps;
    } else if (identifier == "alpha_dist") {
        return m_warp_params.alpha_dist;
    } else if (identifier == "f0min") {
        return m_f0_params.f0min;
    } else if (identifier == "num_octs") {
        return m_f0_params.num_octs;
    } else if (identifier == "f0s_per_oct") {
        return m_f0_params.num_f0s_per_oct;
    } else if (identifier == "f0_prefer_fun") {
        return m_f0_params.prefer ? 1.0 : 0.0;
    } else if (identifier == "f0_prefer_mean") {
        return m_f0_params.prefer_mean;
    } else if (identifier == "f0_prefer_stdev") {
        return m_f0_params.prefer_stdev;
    } else if (identifier == "f0gram_mode") {
        return m_f0gram_mode == BestBinOfAllDirections ? 1.0 : 0.0;
    } else {
        return 0.f;
    }

}

void FChTransformF0gram::setParameter(string identifier, float value)
{
    if (identifier == "fmax") {
        m_fmax = value;
    } else if (identifier == "nsamp_ix") {
        int n = int(roundf(value));
        for (int i = 0; i < int(nsamp_options.size()); ++i) {
            if (i == n) {
                m_warp_params.nsamps_twarp = nsamp_options[i];
                m_blockSize = m_warp_params.nsamps_twarp * 4;
            }
        }
    } else if (identifier == "alpha_max") {
        m_warp_params.alpha_max = value;
    } else if (identifier == "num_warps") {
        m_warp_params.num_warps = value;
    } else if (identifier == "alpha_dist") {
        m_warp_params.alpha_dist = value;
    } else if (identifier == "f0min") {
        m_f0_params.f0min = value;
    } else if (identifier == "num_octs") {
        m_f0_params.num_octs = value;
    } else if (identifier == "f0s_per_oct") {
        m_f0_params.num_f0s_per_oct = value;
    } else if (identifier == "f0_prefer_fun") {
        m_f0_params.prefer = (value > 0.5);
    } else if (identifier == "f0_prefer_mean") {
        m_f0_params.prefer_mean = value;
    } else if (identifier == "f0_prefer_stdev") {
        m_f0_params.prefer_stdev = value;
    } else if (identifier == "f0gram_mode") {
        m_f0gram_mode = (value > 0.5 ?
                         BestBinOfAllDirections :
                         AllBinsOfBestDirection);
    } else {
        cerr << "WARNING: Unknown parameter id \""
             << identifier << "\"" << endl;
    }
}

FChTransformF0gram::ProgramList
FChTransformF0gram::getPrograms() const {
    ProgramList list;
    return list;
}

FChTransformF0gram::OutputList
FChTransformF0gram::getOutputDescriptors() const {

    OutputList list;

    vector<string> labels;
    char label[100];

    if (m_processingMode == ModeF0Gram) {

        /* f0 values of F0gram grid as string values */
        for (int i = 0; i < m_num_f0s; ++i) {
            sprintf(label, "%4.2f Hz", m_f0s[i]);
            labels.push_back(label);
        }
    
        /* The F0gram */
        OutputDescriptor d;
        d.identifier = "f0gram";
        d.name = "F0gram";
        d.description = "The salience of the different f0s in the signal.";
        d.hasFixedBinCount = true;
        d.binCount = m_f0_params.num_octs * m_f0_params.num_f0s_per_oct;
        d.binNames = labels;
        d.hasKnownExtents = false;
        d.isQuantized = false;
        d.sampleType = OutputDescriptor::OneSamplePerStep;
        d.hasDuration = false;
        list.push_back(d);

        d.identifier = "pitch";
        d.name = "Most salient pitch";
        d.description = "The most salient f0 in the signal for each time step.";
        d.unit = "Hz";
        d.hasFixedBinCount = true;
        d.binCount = 1;
        d.binNames.clear();
        d.hasKnownExtents = false;
        d.isQuantized = false;
        d.sampleType = OutputDescriptor::OneSamplePerStep;
        d.hasDuration = false;
        list.push_back(d);

    } else {

        for (int i = 0; i < m_warp_params.nsamps_twarp/2+1; ++i) {
            double freq = i * (m_warpings.fs_warp / m_warp_params.nsamps_twarp);
            sprintf(label, "%4.2f Hz", freq);
            labels.push_back(label);
        }
    
        OutputDescriptor d;
        d.identifier = "spectrogram";
        d.name = "Spectrogram";
        d.description = "Time/frequency spectrogram derived from the Fan Chirp Transform output";
        d.hasFixedBinCount = true;
        d.binCount = m_warp_params.nsamps_twarp/2+1;
        d.binNames = labels;
        d.hasKnownExtents = false;
        d.isQuantized = false;
        d.sampleType = OutputDescriptor::OneSamplePerStep;
        d.hasDuration = false;
        list.push_back(d);
    }
    
    return list;
}

bool
FChTransformF0gram::initialise(size_t channels, size_t stepSize, size_t blockSize) {
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount() ||
        blockSize != m_blockSize/2 ||
        stepSize != m_stepSize) {
        return false;
    }

    m_inputBuffer = allocate_and_zero<float>(m_blockSize);
    
    // WARNING !!!
    // these values in fact are determined by the sampling frequency m_fs
    // the parameters used below correspond to default values i.e. m_fs = 44.100 Hz
    //m_blockSize = 4 * m_warp_params.nsamps_twarp;
//    m_stepSize = floor(m_hop / m_warp_params.fact_over_samp);

    /* design of FChT */
    design_FChT();

    /* initialise m_glogs_params */
    design_GLogS();

    design_LPF();

    design_time_window();

    // Create Hanning window for warped signals
    mp_HanningWindow = allocate<double>(m_warp_params.nsamps_twarp);
    bool normalize = false;
    Utils::hanning_window(mp_HanningWindow, m_warp_params.nsamps_twarp, normalize);

    m_num_f0s = m_f0_params.num_octs * m_f0_params.num_f0s_per_oct;
    m_f0s = allocate<double>(m_num_f0s);
    for (int i = 0; i < m_num_f0s; ++i) {
        m_f0s[i] = m_glogs_f0[m_glogs_init_f0s + i];
    }

    m_initialised = true;
    return true;
}

void
FChTransformF0gram::design_GLogS() {

    // total number & initial quantity of f0s

    m_glogs_init_f0s = (int)(((double)m_f0_params.num_f0s_per_oct)*log2(5.0))+1;
    m_glogs_num_f0s = (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct + m_glogs_init_f0s;

    // Initialize arrays
    m_glogs_f0 = allocate<double>(m_glogs_num_f0s);
    m_glogs = allocate<double>(m_glogs_num_f0s*m_warp_params.num_warps);
    m_glogs_n = allocate<int>(m_glogs_num_f0s);
    m_glogs_index = allocate<int>(m_glogs_num_f0s);

    // Compute f0 values
    m_glogs_harmonic_count = 0;
    double factor = (double)(m_warp_params.nsamps_twarp/2)/(double)(m_warp_params.nsamps_twarp/2+1);
    for (int i = 0; i < m_glogs_num_f0s; i++) {
        m_glogs_f0[i] = (m_f0_params.f0min/5.0)*pow(2.0,(double)i/(double)m_f0_params.num_f0s_per_oct);
        // for every f0 compute number of partials less or equal than m_fmax.
        m_glogs_n[i] = m_fmax*factor/m_glogs_f0[i];
        m_glogs_index[i] = m_glogs_harmonic_count; 
        m_glogs_harmonic_count += m_glogs_n[i];
    }

    // Initialize arrays for interpolation
    m_glogs_posint = allocate<int>(m_glogs_harmonic_count);
    m_glogs_posfrac = allocate<double>(m_glogs_harmonic_count);
    m_glogs_interp = allocate<double>(m_glogs_harmonic_count);

    // Compute int & frac of interpolation positions
    int aux_index = 0;
    double aux_pos;
    for (int i = 0; i < m_glogs_num_f0s; i++) {
        for (int j = 1; j <= m_glogs_n[i]; j++) {
            aux_pos = ((double)j * m_glogs_f0[i]) * ((double)(m_warp_params.nsamps_twarp))/m_warpings.fs_warp;
            m_glogs_posint[aux_index] = (int)aux_pos;
            m_glogs_posfrac[aux_index] = aux_pos - (double)m_glogs_posint[aux_index];
            aux_index++;
        }
    }

    // Third harmonic attenuation
    double aux_third_harmonic;
    m_glogs_third_harmonic_posint = allocate<int>((m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
    m_glogs_third_harmonic_posfrac = allocate<double>((m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
    for (int i = 0; i < (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct; i++) {
        aux_third_harmonic = (double)i + (double)m_glogs_init_f0s - ((double)m_f0_params.num_f0s_per_oct)*log2(3.0);
        m_glogs_third_harmonic_posint[i] = (int)aux_third_harmonic;
        m_glogs_third_harmonic_posfrac[i] = aux_third_harmonic - (double)(m_glogs_third_harmonic_posint[i]);
    }
    m_glogs_third_harmonic = allocate<double>((m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);

    // Fifth harmonic attenuation
    double aux_fifth_harmonic;
    m_glogs_fifth_harmonic_posint = allocate<int>((m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
    m_glogs_fifth_harmonic_posfrac = allocate<double>((m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
    for (int i = 0; i < (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct; i++) {
        aux_fifth_harmonic = (double)i + (double)m_glogs_init_f0s - ((double)m_f0_params.num_f0s_per_oct)*log2(5.0);
        m_glogs_fifth_harmonic_posint[i] = (int)aux_fifth_harmonic;
        m_glogs_fifth_harmonic_posfrac[i] = aux_fifth_harmonic - (double)(m_glogs_fifth_harmonic_posint[i]);
    }
    m_glogs_fifth_harmonic = allocate<double>((m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);

    // Normalization & attenuation windows
    m_glogs_f0_preference_weights = allocate<double>(m_f0_params.num_octs*m_f0_params.num_f0s_per_oct);
    m_glogs_median_correction = allocate<double>(m_f0_params.num_octs*m_f0_params.num_f0s_per_oct);
    m_glogs_sigma_correction = allocate<double>(m_f0_params.num_octs*m_f0_params.num_f0s_per_oct);
    double MIDI_value;
    for (int i = 0; i < m_f0_params.num_octs*m_f0_params.num_f0s_per_oct; i++) {
        if (m_f0_params.prefer) {
            MIDI_value = 69.0 + 12.0 * log2(m_glogs_f0[i + m_glogs_init_f0s]/440.0);
            m_glogs_f0_preference_weights[i] = 1.0/sqrt(2.0*M_PI*m_f0_params.prefer_stdev*m_f0_params.prefer_stdev)*exp(-(MIDI_value-m_f0_params.prefer_mean)*(MIDI_value-m_f0_params.prefer_mean)/(2.0*m_f0_params.prefer_stdev*m_f0_params.prefer_stdev));
            m_glogs_f0_preference_weights[i] = (0.01 + m_glogs_f0_preference_weights[i]) / (1.01);
        } else {
            m_glogs_f0_preference_weights[i] = 1.0;
        }
		
        m_glogs_median_correction[i] = m_glogs_params.median_poly_coefs[0]*(i+1.0)*(i+1.0) + m_glogs_params.median_poly_coefs[1]*(i+1.0) + m_glogs_params.median_poly_coefs[2];
        m_glogs_sigma_correction[i] = 1.0 / (m_glogs_params.sigma_poly_coefs[0]*(i+1.0)*(i+1.0) + m_glogs_params.sigma_poly_coefs[1]*(i+1.0) + m_glogs_params.sigma_poly_coefs[2]);
    }
}

void
FChTransformF0gram::design_FChT() {

    /*  =============  WARPING DESIGN   ============= */

    // sampling frequency after oversampling
    m_warpings.fs_orig = m_warp_params.fact_over_samp * m_fs;

    // number of samples of the original signal frame
    m_warpings.nsamps_torig = 4 * m_warp_params.fact_over_samp * m_warp_params.nsamps_twarp;
    // equivalent to: m_warpings.nsamps_torig = m_warp_params.fact_over_samp * m_blockSize;

    // time instants of the original signal frame
    double *t_orig = allocate<double>(m_warpings.nsamps_torig);
    for (int ind = 0; ind < m_warpings.nsamps_torig; ind++) {
        t_orig[ind] = ((double)(ind + 1) - (double)m_warpings.nsamps_torig / 2.0) / m_warpings.fs_orig;
    }

    // linear chirps warping definition as relative frequency deviation
    //TODO
    double *freq_relative = allocate<double>(m_warpings.nsamps_torig * m_warp_params.num_warps);
    define_warps_linear_chirps(freq_relative, t_orig);

    // maximum relative frequency deviation
    double freq_relative_max = 0;
    for (int i = 0; i < m_warpings.nsamps_torig; i++) {
        for (int j = 0; j < m_warp_params.num_warps; j++) {
            if (freq_relative_max < freq_relative[j * m_warpings.nsamps_torig + i]) {
                freq_relative_max = freq_relative[j * m_warpings.nsamps_torig + i];
            }
        }
    }

    // sampling frequency of warped signal to be free of aliasing up to fmax
    m_warpings.fs_warp = 2 * m_fmax * freq_relative_max;

    // time instants of the warped signal frame
    double *t_warp = allocate<double>(m_warp_params.nsamps_twarp);
    for (int ind = 0; ind < m_warp_params.nsamps_twarp; ind++) {
        t_warp[ind] = ((double)((int)(ind + 1)- (int)m_warp_params.nsamps_twarp / 2)) / (double)m_warpings.fs_warp;
    }

    // design of warpings for efficient interpolation
    design_warps(freq_relative, t_orig, t_warp);

    deallocate(freq_relative);
    deallocate(t_orig);
    deallocate(t_warp);

    /*  =============  FFTW PLAN DESIGN   ============= */
    // Initialize 2-d array for warped signals
    x_warping = allocate<double>(m_warp_params.nsamps_twarp);
    m_absFanChirpTransform = allocate<double>(m_warp_params.num_warps * (m_warp_params.nsamps_twarp/2 + 1));
    m_auxFanChirpTransform = allocate<double>(2 * (m_warp_params.nsamps_twarp/2 + 1));
    fft_xwarping = new FFTReal(m_warp_params.nsamps_twarp);
}

void
FChTransformF0gram::design_warps(double * freq_relative, double * t_orig, double * t_warp) {
    /* the warping is done by interpolating the original signal in time instants
       given by the desired frequency deviation, to do this, the interpolation
       instants are stored in a structure as an integer index and a fractional value
       hypothesis: sampling frequency at the central point equals the original
    */

    m_warpings.pos_int = allocate<int>(m_warp_params.num_warps * m_warp_params.nsamps_twarp);
    m_warpings.pos_frac = allocate<double>(m_warp_params.num_warps * m_warp_params.nsamps_twarp);

    // vector of phase values
    double *phi = allocate<double>(m_warpings.nsamps_torig);
    double aux;

    // warped positions
    double *pos1 = allocate<double>(m_warp_params.nsamps_twarp*m_warp_params.num_warps);
	
    for (int i = 0; i < m_warp_params.num_warps; i++) {
		
        // integration of relative frequency to obtain phase values
        Utils::cumtrapz(t_orig, freq_relative + i*(m_warpings.nsamps_torig), m_warpings.nsamps_torig, phi);

        // centering of phase values to force original frequency in the middle
        aux = phi[m_warpings.nsamps_torig/2];
        for (int j = 0; j < m_warpings.nsamps_torig; j++) {
            phi[j] -= aux;
        } //for

        // interpolation of phase values to obtain warped positions
        Utils::interp1(phi, t_orig, m_warpings.nsamps_torig, t_warp, pos1 + i*m_warp_params.nsamps_twarp, m_warp_params.nsamps_twarp);
    }

    // % previous sample index
    // pos1_int = uint32(floor(pos1))';
    // % integer corresponding to previous sample index in "c"
    // warps.pos1_int = (pos1_int - uint32(1));
    // % fractional value that defines the warped position
    // warps.pos1_frac = (double(pos1)' - double(pos1_int));

    for (int j = 0; j < m_warp_params.nsamps_twarp*m_warp_params.num_warps; j++) {
        // previous sample index
        pos1[j] = pos1[j]*m_warpings.fs_orig + m_warpings.nsamps_torig/2 + 1;
        m_warpings.pos_int[j] = (int) pos1[j];
        m_warpings.pos_frac[j] = pos1[j] - (double)(m_warpings.pos_int[j]);
    } //for

    deallocate(phi);
    deallocate(pos1);
}

void
FChTransformF0gram::define_warps_linear_chirps(double * freq_relative, double * t_orig) {
    /**  define warps as relative frequency deviation from original frequency
         t_orig : time vector
         freq_relative : relative frequency deviations
    */
    if (m_warp_params.alpha_dist == 0) {

        // linear alpha values spacing
        m_warpings.chirp_rates = allocate<double>(m_warp_params.num_warps);
        // WARNING m_warp_params.num_warps must be odd
        m_warpings.chirp_rates[0] = -m_warp_params.alpha_max;
        double increment = (double) m_warp_params.alpha_max / ((m_warp_params.num_warps - 1) / 2);

        for (int ind = 1; ind < m_warp_params.num_warps; ind++) {
            m_warpings.chirp_rates[ind] = m_warpings.chirp_rates[ind - 1] + increment;
        }
        // force zero value
        m_warpings.chirp_rates[(int) ((m_warp_params.num_warps - 1) / 2)] = 0;

    } else {
        // log alpha values spacing
        m_warpings.chirp_rates = allocate<double>(m_warp_params.num_warps);

        // force zero value
        int middle_point = (int) ((m_warp_params.num_warps - 1) / 2);
        m_warpings.chirp_rates[middle_point] = 0;

        double logMax = log10(m_warp_params.alpha_max + 1);
        double increment = logMax / ((m_warp_params.num_warps - 1) / 2.0f);
        double exponent = 0;

        // fill positive values
        int ind_log = middle_point;
        for (int ind = 0; ind < (m_warp_params.num_warps + 1) / 2; ind++) {
            m_warpings.chirp_rates[ind_log] = pow(10, exponent) - 1;
            exponent += increment;
            ind_log++;
        }
        // fill negative values
        for (int ind = 0; ind < (m_warp_params.num_warps - 1) / 2; ind++) {
            m_warpings.chirp_rates[ind] = -m_warpings.chirp_rates[m_warp_params.num_warps - 1 - ind];
        }
    }

    // compute relative frequency deviation
    for (int i = 0; i < m_warpings.nsamps_torig; i++) {
        for (int j = 0; j < m_warp_params.num_warps; j++) {
            freq_relative[j * m_warpings.nsamps_torig + i] = 1.0 + t_orig[i] * m_warpings.chirp_rates[j];
        }
    }
}

void
FChTransformF0gram::design_LPF()
{
    double *lp_LPFWindow_aux = allocate<double>(m_blockSize/2+1);
    mp_LPFWindow = allocate<double>(m_blockSize/2+1);
    
    int i_max = (int) ((2.0*m_fmax/m_fs) * ( (double)m_blockSize / 2.0 + 1.0 ));
    for (int i = 0; i < m_blockSize/2+1; i++) {
        if (i >= i_max) {
            lp_LPFWindow_aux[i] = 0.0;
        } else {
            lp_LPFWindow_aux[i] = 1.0;
        }
    }

    LPF_time = allocate_and_zero<double>(m_warpings.nsamps_torig);
    LPF_frequency = allocate_and_zero<double>(2 * (m_warpings.nsamps_torig/2 + 1));

    fft_forward_LPF = new FFTReal(m_blockSize);
    fft_inverse_LPF = new FFTReal(m_warpings.nsamps_torig);
    
    int winWidth = 11;
    double *lp_hanningWindow = allocate<double>(winWidth); 
    double accum=0;
    for (int i = 0; i < winWidth; i++) {
        lp_hanningWindow[i]=0.5*(1.0-cos(2*M_PI*(double)(i+1)/((double)winWidth+1.0)));
        accum+=lp_hanningWindow[i];
        
    }
    for (int i = 0; i < winWidth; i++) { //window normalization
        lp_hanningWindow[i]=lp_hanningWindow[i]/accum;
    }
    for (int i = 0; i < m_blockSize/2+1; i++) {
        //if (((i-(winWidth-1)/2)<0)||(i+(winWidth-1))/2>m_blockSize/2-1) {//consideramos winWidth impar, si la ventana sale del arreglo se rellena con el valor origianl
        if ( (i > (i_max + (winWidth-1)/2)) ||  (i <= (i_max - (winWidth-1)/2)) ) {
            mp_LPFWindow[i]=lp_LPFWindow_aux[i];
        } else {
            accum=0;
            for (int j = -((winWidth-1)/2); j <= (winWidth-1)/2; j++) {
            	accum+=lp_LPFWindow_aux[i-j]*lp_hanningWindow[j+(winWidth-1)/2];
            }
            mp_LPFWindow[i]=accum;
        }
    }

    deallocate(lp_LPFWindow_aux);
    deallocate(lp_hanningWindow);
}

void FChTransformF0gram::apply_LPF()
{
    fft_forward_LPF->forward(LPF_time, LPF_frequency);

    for (int i = 0; i < m_blockSize/2+1; i++) {
        LPF_frequency[i*2]     *= mp_LPFWindow[i];
        LPF_frequency[i*2 + 1] *= mp_LPFWindow[i];
    }

    fft_inverse_LPF->inverse(LPF_frequency, LPF_time);
    
    // TODO ver si hay que hacer fftshift para corregir la fase respecto al centro del frame.
    // nota: además de aplicar el LPF, esta función resamplea la señal original.
}

void FChTransformF0gram::clean_LPF()
{
    delete fft_forward_LPF;
    delete fft_inverse_LPF;
    deallocate(LPF_time);
    deallocate(LPF_frequency);
    deallocate(mp_LPFWindow);
}

void FChTransformF0gram::reset()
{
}

FChTransformF0gram::FeatureSet
FChTransformF0gram::process(const float *const *inputBuffers, Vamp::RealTime) {

    if (!m_initialised) return FeatureSet();
    
    /* PSEUDOCÓDIGO:
       - Aplicar FFT al frame entero.	
       - Filtro pasabajos en frecuencia.
       - FFT inversa al frame entero.
       -----------------------------------------------------------------------------
       - Para cada warp: *Si es un espectrograma direccional (un solo warp 
       => no es para cada warp sino para el elegido)
       - Hacer la interpolación con interp1q.
       - Aplicar la FFT al frame warpeado.
       - (Opcional) GLogS.
       - ...
    */

//---------------------------------------------------------------------------
    FeatureSet fs;
	
#ifdef DEBUG
    fprintf(stderr, "\n	----- DEBUG INFORMATION ----- \n");
    fprintf(stderr, "	m_fs = %f Hz.\n",m_fs);
    fprintf(stderr, "	fs_orig = %f Hz.\n",m_warpings.fs_orig);
    fprintf(stderr, "	fs_warp = %f Hz.\n",m_warpings.fs_warp);
    fprintf(stderr, "	m_blockSize = %d.\n",m_blockSize);
    fprintf(stderr, "	m_warpings.nsamps_torig = %d.\n",m_warpings.nsamps_torig);
    fprintf(stderr, "	m_warp_params.nsamps_twarp = %d.\n",m_warp_params.nsamps_twarp);
    fprintf(stderr, "	m_warp_params.num_warps = %d.\n",m_warp_params.num_warps);
    fprintf(stderr, "	m_glogs_harmonic_count = %d.\n",m_glogs_harmonic_count);
#endif

    for (int i = 0; i < m_blockSize - m_stepSize; ++i) {
        m_inputBuffer[i] = m_inputBuffer[i + m_stepSize];
    }
    for (int i = 0; i < m_blockSize/2; ++i) {
        m_inputBuffer[m_blockSize/2 + i] = inputBuffers[0][i];
    }
    for (int i = 0; i < m_blockSize; ++i) {
        LPF_time[i] = m_inputBuffer[i] * m_timeWindow[i];
    }
    for (int i = 0; i < m_blockSize; ++i) {
        LPF_time[m_blockSize + i] = 0.0;
    }
        
    apply_LPF();
    // Señal filtrada queda en LPF_time

    Feature feature;
    feature.hasTimestamp = false;

    if (m_processingMode == ModeRoughSpectrogram) {
        feature.values = vector<float>(m_warp_params.nsamps_twarp/2+1, 0.f);
    }

// ----------------------------------------------------------------------------------------------
// 		Hanning window & FFT for all warp directions

    double max_glogs = -DBL_MAX;
    int ind_max_glogs = 0;

    for (int i_warp = 0; i_warp < m_warp_params.num_warps; i_warp++) {
        
        // Interpolate
        Utils::interp1q(LPF_time, (m_warpings.pos_int) + i_warp*m_warp_params.nsamps_twarp, m_warpings.pos_frac + i_warp*m_warp_params.nsamps_twarp, x_warping, m_warp_params.nsamps_twarp);

        // Apply window
        for (int i = 0; i < m_warp_params.nsamps_twarp; i++) {
            x_warping[i] *= mp_HanningWindow[i];
        }

        // Transform
        fft_xwarping->forward(x_warping, m_auxFanChirpTransform);

        if (m_processingMode == ModeRoughSpectrogram) {
            for (int i = 0; i < (m_warp_params.nsamps_twarp/2+1); i++) {
                double abs = sqrt(m_auxFanChirpTransform[i*2]*m_auxFanChirpTransform[i*2]+m_auxFanChirpTransform[i*2+1]*m_auxFanChirpTransform[i*2+1]);
                if (abs > feature.values[i]) {
                    feature.values[i] = abs;
                }
            }
            continue;
        }

        // Copy result
        double *aux_abs_fcht = m_absFanChirpTransform + i_warp*(m_warp_params.nsamps_twarp/2+1);
        for (int i = 0; i < (m_warp_params.nsamps_twarp/2+1); i++) {
            aux_abs_fcht[i] = log10(1.0 + 10.0*sqrt(m_auxFanChirpTransform[i*2]*m_auxFanChirpTransform[i*2]+m_auxFanChirpTransform[i*2+1]*m_auxFanChirpTransform[i*2+1]));
        }
		
//      -----------------------------------------------------------------------------------------
// 		GLogS
        Utils::interp1q(aux_abs_fcht, m_glogs_posint, m_glogs_posfrac, m_glogs_interp, m_glogs_harmonic_count);
        int glogs_ind = 0;
        for (int i = 0; i < m_glogs_num_f0s; i++) {
            double glogs_accum = 0;
            for (int j = 1; j <= m_glogs_n[i]; j++) {
                glogs_accum += m_glogs_interp[glogs_ind++];
            }
            m_glogs[i + i_warp*m_glogs_num_f0s] = glogs_accum/(double)m_glogs_n[i];
        }

//		Sub/super harmonic correction
        Utils::interp1q(m_glogs + i_warp*m_glogs_num_f0s, m_glogs_third_harmonic_posint, m_glogs_third_harmonic_posfrac, m_glogs_third_harmonic, (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
        Utils::interp1q(m_glogs + i_warp*m_glogs_num_f0s, m_glogs_fifth_harmonic_posint, m_glogs_fifth_harmonic_posfrac, m_glogs_fifth_harmonic, (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
        for (int i = m_glogs_num_f0s-1; i >= m_glogs_init_f0s; i--) {
            m_glogs[i + i_warp*m_glogs_num_f0s] -= MAX(MAX(m_glogs[i-m_f0_params.num_f0s_per_oct + i_warp*m_glogs_num_f0s],m_glogs_third_harmonic[i-m_glogs_init_f0s]),m_glogs_fifth_harmonic[i-m_glogs_init_f0s]);
        }
        for (int i = m_glogs_init_f0s; i < m_glogs_num_f0s-m_f0_params.num_f0s_per_oct; i++) {
            m_glogs[i + i_warp*m_glogs_num_f0s] -= 0.3*m_glogs[i+m_f0_params.num_f0s_per_oct + i_warp*m_glogs_num_f0s];
            // Median, sigma $ weights correction
            m_glogs[i + i_warp*m_glogs_num_f0s] = (m_glogs[i + i_warp*m_glogs_num_f0s]-m_glogs_median_correction[i-m_glogs_init_f0s])*m_glogs_sigma_correction[i-m_glogs_init_f0s]*m_glogs_f0_preference_weights[i-m_glogs_init_f0s];
        }
	
        // Look for maximum value to determine best direction
        for (int i = m_glogs_init_f0s; i < m_glogs_num_f0s-m_f0_params.num_f0s_per_oct; i++) {
            if (m_glogs[i + i_warp*m_glogs_num_f0s] > max_glogs) {
                max_glogs = m_glogs[i + i_warp*m_glogs_num_f0s];
                ind_max_glogs = i_warp;
            }
        }
    }	

    if (m_processingMode == ModeRoughSpectrogram) {

        // already accumulated our return values in feature
        fs[0].push_back(feature);

    } else if (m_processingMode == ModeSpectrogram) {

        for (int i = 0; i < m_warp_params.nsamps_twarp/2+1; i++) {
            feature.values.push_back(pow(10.0, m_absFanChirpTransform[ind_max_glogs * (m_warp_params.nsamps_twarp/2+1) + i]) - 1.0);
        }
        fs[0].push_back(feature);

    } else { // f0gram

        int bestIndex = -1;
        
        for (int i=m_glogs_init_f0s; i< m_glogs_num_f0s - m_f0_params.num_f0s_per_oct; i++) {
            double value = 0.0;
            switch (m_f0gram_mode) {
            case AllBinsOfBestDirection:
                value = m_glogs[i+(int)ind_max_glogs*(int)m_glogs_num_f0s];
                break;
            case BestBinOfAllDirections:
                max_glogs = -DBL_MAX;
                for (int i_warp = 0; i_warp < m_warp_params.num_warps; i_warp++) {
                    if (m_glogs[i + i_warp*m_glogs_num_f0s] > max_glogs) {
                        max_glogs = m_glogs[i + i_warp*m_glogs_num_f0s];
                        ind_max_glogs = i_warp;
                    }
                }
                value = max_glogs;
                break;
            }
            if (bestIndex < 0 || float(value) > feature.values[bestIndex]) {
                bestIndex = int(feature.values.size());
            }
            feature.values.push_back(float(value));
        }
        
        fs[0].push_back(feature);

        if (bestIndex >= 0) {

            double bestValue = feature.values[bestIndex];
            set<double> ordered(feature.values.begin(), feature.values.end());
            vector<double> flattened(ordered.begin(), ordered.end());
            double median = flattened[flattened.size()/2];
            if (bestValue > median * 8.0) {
                Feature pfeature;
                pfeature.hasTimestamp = false;
                pfeature.values.push_back(m_f0s[bestIndex]);
                fs[1].push_back(pfeature);
            }
        }
    }

    return fs;
}

FChTransformF0gram::FeatureSet
FChTransformF0gram::getRemainingFeatures() {
    return FeatureSet();
}

void
FChTransformF0gram::design_time_window() {

    int transitionWidth = (int)m_blockSize/128 + 128;
    m_timeWindow = allocate<double>(m_blockSize);
    double *lp_transitionWindow = allocate<double>(transitionWidth);

    for (int i = 0; i < m_blockSize; i++) {
        m_timeWindow[i] = 1.0;
    }

    for (int i = 0; i < transitionWidth; i++) {
        lp_transitionWindow[i]=0.5*(1.0-cos(2*M_PI*(double)(i+1)/((double)transitionWidth+1.0)));
    }

    for (int i = 0; i < transitionWidth/2; i++) {
        m_timeWindow[i] = lp_transitionWindow[i];
        m_timeWindow[m_blockSize-1-i] = lp_transitionWindow[transitionWidth-1-i];
    }

#ifdef DEBUG
    for (int i = 0; i < m_blockSize; i++) {
        if ((i<transitionWidth)) {
            fprintf(stderr, "	m_timeWindow[%d] = %f.\n",i,m_timeWindow[i]);
        }
    }
#endif

    deallocate(lp_transitionWindow);
}

