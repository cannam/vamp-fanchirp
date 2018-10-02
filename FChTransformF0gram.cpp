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
//#define DEBUG
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

FChTransformF0gram::FChTransformF0gram(float inputSampleRate) :
Plugin(inputSampleRate),
m_currentProgram("default"),
m_stepSize(0), // We are using 0 for step and block size to indicate "not yet set".
m_blockSize(0) {

    m_fs = inputSampleRate;
    // max frequency of interest (Hz)
    m_fmax = 10000.f;
    // warping parameters
    m_warp_params.nsamps_twarp   = 2048;
    //m_warp_params.nsamps_twarp = 8;
    m_warp_params.alpha_max = 4;
    m_warp_params.num_warps = 21;
    //m_warp_params.num_warps = 11;
    m_warp_params.fact_over_samp = 2;
    m_warp_params.alpha_dist = 0;
    // f0 parameters
    m_f0_params.f0min = 80.0;
    m_f0_params.num_octs = 4;
    m_f0_params.num_f0s_per_oct = 192;
    m_f0_params.num_f0_hyps = 5;
    m_f0_params.prefer = true;
    m_f0_params.prefer_mean = 60;
    m_f0_params.prefer_stdev = 18;
    // glogs parameters
    m_glogs_params.HP_logS = true;
    m_glogs_params.att_subharms = 1;
	// display parameters
	m_f0gram_mode = true;

    m_glogs_params.median_poly_coefs[0] = -0.000000058551680;
    m_glogs_params.median_poly_coefs[1] = -0.000006945207775;
    m_glogs_params.median_poly_coefs[2] = 0.002357223226588;

    m_glogs_params.sigma_poly_coefs[0] = 0.000000092782308;
    m_glogs_params.sigma_poly_coefs[1] = 0.000057283574898;
    m_glogs_params.sigma_poly_coefs[2] = 0.022199903714288;

    // number of fft points (controls zero-padding)
    m_nfft = m_warp_params.nsamps_twarp;
    // hop in samples
    m_hop = m_warp_params.fact_over_samp * 256;

    m_num_f0s = 0;

}

FChTransformF0gram::~FChTransformF0gram() {
    // remeber to delete everything that deserves to
}

string
FChTransformF0gram::getIdentifier() const {
    return "fchtransformf0gram";
}

string
FChTransformF0gram::getName() const {
    return "Fan Chirp Transform F0gram";
}

string
FChTransformF0gram::getDescription() const {
    // Return something helpful here!
    return "This plug-in produces a representation, called F0gram, which exhibits the salience of the fundamental frequency of the sound sources in the audio file. The computation of the F0gram makes use of the Fan Chirp Transform analysis. It is based on the article \"Fan chirp transform for music representation\"  P. Cancela, E. Lopez, M. Rocamora, International Conference on Digital Audio Effects, 13th. DAFx-10. Graz, Austria - 6-10 Sep 2010.";
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
    return 0;
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
    return 8192; // 0 means "I can handle any block size"
}

size_t
FChTransformF0gram::getPreferredStepSize() const {
    return 256; // 0 means "anything sensible"; in practice this
    // means the same as the block size for TimeDomain
    // plugins, or half of it for FrequencyDomain plugins
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
    nsamp.identifier = "nsamp";
    nsamp.name = "Number of samples";
    nsamp.description = "Number of samples of the time warped frame";
    nsamp.unit = "samples";
    nsamp.minValue = 128;
    nsamp.maxValue = 4096;
    nsamp.defaultValue = 2048;
    nsamp.isQuantized = true;
    nsamp.quantizeStep = 1.0;
    list.push_back(nsamp);

    ParameterDescriptor nfft;
    nfft.identifier = "nfft";
    nfft.name = "FFT number of points";
    nfft.description = "Number of FFT points (controls zero-padding)";
    nfft.unit = "samples";
    nfft.minValue = 0;
    nfft.maxValue = 4;
    nfft.defaultValue = 3;
    nfft.isQuantized = true;
    nfft.quantizeStep = 1.0;
    nfft.valueNames.push_back("256");
    nfft.valueNames.push_back("512");
    nfft.valueNames.push_back("1024");
    nfft.valueNames.push_back("2048");
    nfft.valueNames.push_back("4096");
    nfft.valueNames.push_back("8192");
    list.push_back(nfft);

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

    ParameterDescriptor num_f0_hyps;
    num_f0_hyps.identifier = "num_f0_hyps";
    num_f0_hyps.name = "number of f0 hypotesis";
    num_f0_hyps.description = "Number of f0 hypotesis to extract.";
    num_f0_hyps.unit = "";
    num_f0_hyps.minValue = 1;
    num_f0_hyps.maxValue = 100;
    num_f0_hyps.defaultValue = 10;
    num_f0_hyps.isQuantized = true;
    num_f0_hyps.quantizeStep = 1.0;
    list.push_back(num_f0_hyps);

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
    f0_prefer_fun.name = "f0 preference function";
    f0_prefer_fun.description = "Whether to use a f0 weighting function.";
    f0_prefer_fun.unit = "";
    f0_prefer_fun.minValue = 0;
    f0_prefer_fun.maxValue = 1;
    f0_prefer_fun.defaultValue = 1;
    f0_prefer_fun.isQuantized = true;
    f0_prefer_fun.quantizeStep = 1.0;
    list.push_back(f0_prefer_fun);

    ParameterDescriptor f0_prefer_mean;
    f0_prefer_mean.identifier = "f0_prefer_mean";
    f0_prefer_mean.name = "mean f0 preference function";
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
    f0_prefer_stdev.name = "stdev of f0 preference function";
    f0_prefer_stdev.description = "Stdev for f0 weighting function (MIDI number).";
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
    } else if (identifier == "nsamp") {
        return m_warp_params.nsamps_twarp;
    } else if (identifier == "alpha_max") {
        return m_warp_params.alpha_max;
    } else if (identifier == "num_warps") {
        return m_warp_params.num_warps;
    } else if (identifier == "alpha_dist") {
        return m_warp_params.alpha_dist;
    } else if (identifier == "nfft") {
        return m_nfft;
    } else if (identifier == "f0min") {
        return m_f0_params.f0min;
    } else if (identifier == "num_octs") {
        return m_f0_params.num_octs;
    } else if (identifier == "f0s_per_oct") {
        return m_f0_params.num_f0s_per_oct;
    } else if (identifier == "num_f0_hyps") {
        return m_f0_params.num_f0_hyps;
    } else if (identifier == "f0_prefer_fun") {
        return m_f0_params.prefer;
    } else if (identifier == "f0_prefer_mean") {
        return m_f0_params.prefer_mean;
    } else if (identifier == "f0_prefer_stdev") {
        return m_f0_params.prefer_stdev;
	} else if (identifier == "f0gram_mode") {
        return m_f0gram_mode;
    } else {
        return 0.f;
    }

}

void FChTransformF0gram::setParameter(string identifier, float value) {

    if (identifier == "fmax") {
        m_fmax = value;
    } else if (identifier == "nsamp") {
        m_warp_params.nsamps_twarp = value;
    } else if (identifier == "alpha_max") {
        m_warp_params.alpha_max = value;
    } else if (identifier == "num_warps") {
        m_warp_params.num_warps = value;
    } else if (identifier == "alpha_dist") {
        m_warp_params.alpha_dist = value;
    } else if (identifier == "nfft") {
        m_nfft = value;
    } else if (identifier == "f0min") {
        m_f0_params.f0min = value;
    } else if (identifier == "num_octs") {
        m_f0_params.num_octs = value;
    } else if (identifier == "f0s_per_oct") {
        m_f0_params.num_f0s_per_oct = value;
    } else if (identifier == "num_f0_hyps") {
        m_f0_params.num_f0_hyps = value;
    } else if (identifier == "f0_prefer_fun") {
        m_f0_params.prefer = value;
    } else if (identifier == "f0_prefer_mean") {
        m_f0_params.prefer_mean = value;
    } else if (identifier == "f0_prefer_stdev") {
        m_f0_params.prefer_stdev = value;
    } else if (identifier == "f0gram_mode") {
        m_f0gram_mode = value;
    }

}

FChTransformF0gram::ProgramList
FChTransformF0gram::getPrograms() const {
    ProgramList list;

    list.push_back("default");

    return list;
}

string
FChTransformF0gram::getCurrentProgram() const {
    return m_currentProgram;
}

void
FChTransformF0gram::selectProgram(string name) {

    m_currentProgram = name;

    if (name == "default") {
        m_fmax = 10000.f;

        m_warp_params.nsamps_twarp = 2048;
        m_warp_params.alpha_max = 4;
        m_warp_params.num_warps = 21;
        m_warp_params.fact_over_samp = 2;
        m_warp_params.alpha_dist = 0;

        m_f0_params.f0min = 80.0;
        m_f0_params.num_octs = 4;
        m_f0_params.num_f0s_per_oct = 192;
        m_f0_params.num_f0_hyps = 5;
        m_f0_params.prefer = true;
        m_f0_params.prefer_mean = 60;
        m_f0_params.prefer_stdev = 18;

        m_glogs_params.HP_logS = true;
        m_glogs_params.att_subharms = 1;

        m_glogs_params.median_poly_coefs[0] = -0.000000058551680;
        m_glogs_params.median_poly_coefs[1] = -0.000006945207775;
        m_glogs_params.median_poly_coefs[2] = 0.002357223226588;

        m_glogs_params.sigma_poly_coefs[0] = 0.000000092782308;
        m_glogs_params.sigma_poly_coefs[1] = 0.000057283574898;
        m_glogs_params.sigma_poly_coefs[2] = 0.022199903714288;

        m_nfft = m_warp_params.nsamps_twarp;
        m_hop = m_warp_params.fact_over_samp * 256;

        m_num_f0s = 0;

		m_f0gram_mode = 1;

    }
}

FChTransformF0gram::OutputList
FChTransformF0gram::getOutputDescriptors() const {

    OutputList list;

    // See OutputDescriptor documentation for the possibilities here.
    // Every plugin must have at least one output.

    /* f0 values of F0gram grid as string values */
    vector<string> f0values;
    size_t ind = 0;
    char f0String[10];
    while (ind < m_num_f0s) {
        sprintf(f0String, "%4.2f", m_f0s[ind]);
        f0values.push_back(f0String);
        ind++;
    }

    /* The F0gram */
    OutputDescriptor d;
    d.identifier = "f0gram";
    d.name = "F0gram: salience of f0s";
    d.description = "This representation show the salience of the different f0s in the signal.";
    d.unit = "Hertz";
    d.hasFixedBinCount = true;
    //d.binCount = m_num_f0s;
	//d.binCount = m_blockSize/2+1;
	//d.binCount = m_warp_params.nsamps_twarp/2+1;
	//d.binCount = m_warpings.nsamps_torig;
	d.binCount = m_f0_params.num_octs*m_f0_params.num_f0s_per_oct;
    d.binNames = f0values;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::OneSamplePerStep;
    d.hasDuration = false;
    list.push_back(d);

    return list;
}

bool
FChTransformF0gram::initialise(size_t channels, size_t stepSize, size_t blockSize) {
    if (channels < getMinChannelCount() ||
            channels > getMaxChannelCount()) return false;

    // set blockSize and stepSize (but changed below)
    m_blockSize = blockSize;
    m_stepSize = stepSize;

    // WARNING !!!
    // these values in fact are determined by the sampling frequency m_fs
    // the parameters used below correspond to default values i.e. m_fs = 44.100 Hz
    //m_blockSize = 4 * m_warp_params.nsamps_twarp;
    m_stepSize = floor(m_hop / m_warp_params.fact_over_samp);

    /* initialise m_warp_params  */
    //    FChTF0gram:warping_design m_warpings = new warping_design;
    /* initialise m_f0_params    */

    /* initialise m_glogs_params */
	design_GLogS();

    /* design of FChT */
    // design_fcht(m_warps, m_accums, m_f0s)
    design_FChT();

	design_FFT();

	design_LPF();

	design_time_window();

	// Create Hanning window for warped signals
	mp_HanningWindow = new double[m_warp_params.nsamps_twarp];
	bool normalize = false;
	hanning_window(mp_HanningWindow, m_warp_params.nsamps_twarp, normalize);

    return true;
}

void
FChTransformF0gram::design_GLogS() {

	// total number & initial quantity of f0s
	m_glogs_init_f0s = (size_t)(((double)m_f0_params.num_f0s_per_oct)*log2(5.0))+1;
	m_glogs_num_f0s = (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct + m_glogs_init_f0s;

	// Initialize arrays
	m_glogs_f0 = new double[m_glogs_num_f0s];
	m_glogs = new double[m_glogs_num_f0s*m_warp_params.num_warps];
	m_glogs_n = new size_t[m_glogs_num_f0s];
	m_glogs_index = new size_t[m_glogs_num_f0s];

	// Compute f0 values
	m_glogs_harmonic_count = 0;
	double factor = (double)(m_warp_params.nsamps_twarp/2)/(double)(m_warp_params.nsamps_twarp/2+1);
	for (size_t i = 0; i < m_glogs_num_f0s; i++) {
		m_glogs_f0[i] = (m_f0_params.f0min/5.0)*pow(2.0,(double)i/(double)m_f0_params.num_f0s_per_oct);
		// for every f0 compute number of partials less or equal than m_fmax.
		m_glogs_n[i] = m_fmax*factor/m_glogs_f0[i];
		m_glogs_index[i] = m_glogs_harmonic_count; 
		m_glogs_harmonic_count += m_glogs_n[i];
	}

	// Initialize arrays for interpolation
	m_glogs_posint = new size_t[m_glogs_harmonic_count];
	m_glogs_posfrac = new double[m_glogs_harmonic_count];
	m_glogs_interp = new double[m_glogs_harmonic_count];

	// Compute int & frac of interpolation positions
	size_t aux_index = 0;
	double aux_pos;
	for (size_t i = 0; i < m_glogs_num_f0s; i++) {
		for (size_t j = 1; j <= m_glogs_n[i]; j++) {
			// indice en el vector de largo t_warp/2+1 donde el ultimo valor corresponde a f=m_fmax
			aux_pos = ((double)j*m_glogs_f0[i])*((double)(m_warp_params.nsamps_twarp/2+1))/m_fmax;
			m_glogs_posint[aux_index] = (size_t)aux_pos;
			m_glogs_posfrac[aux_index] = aux_pos - (double)m_glogs_posint[aux_index];
			aux_index++;
		}
	}

	// Third harmonic attenuation
	double aux_third_harmonic;
	m_glogs_third_harmonic_posint = new size_t[(m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct];
	m_glogs_third_harmonic_posfrac = new double[(m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct];
	for (size_t i = 0; i < (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct; i++) {
		aux_third_harmonic = (double)i + (double)m_glogs_init_f0s - ((double)m_f0_params.num_f0s_per_oct)*log2(3.0);
		m_glogs_third_harmonic_posint[i] = (size_t)aux_third_harmonic;
		m_glogs_third_harmonic_posfrac[i] = aux_third_harmonic - (double)(m_glogs_third_harmonic_posint[i]);
	}
	m_glogs_third_harmonic = new double[(m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct];

	// Fifth harmonic attenuation
	double aux_fifth_harmonic;
	m_glogs_fifth_harmonic_posint = new size_t[(m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct];
	m_glogs_fifth_harmonic_posfrac = new double[(m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct];
	for (size_t i = 0; i < (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct; i++) {
		aux_fifth_harmonic = (double)i + (double)m_glogs_init_f0s - ((double)m_f0_params.num_f0s_per_oct)*log2(5.0);
		m_glogs_fifth_harmonic_posint[i] = (size_t)aux_fifth_harmonic;
		m_glogs_fifth_harmonic_posfrac[i] = aux_fifth_harmonic - (double)(m_glogs_fifth_harmonic_posint[i]);
	}
	m_glogs_fifth_harmonic = new double[(m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct];

	// Normalization & attenuation windows
	m_glogs_f0_preference_weights = new double[m_f0_params.num_octs*m_f0_params.num_f0s_per_oct];
	m_glogs_median_correction = new double[m_f0_params.num_octs*m_f0_params.num_f0s_per_oct];
	m_glogs_sigma_correction = new double[m_f0_params.num_octs*m_f0_params.num_f0s_per_oct];
	m_glogs_hf_smoothing_window = new double[m_warp_params.nsamps_twarp/2+1];
	double MIDI_value;
	for (size_t i = 0; i < m_f0_params.num_octs*m_f0_params.num_f0s_per_oct; i++) {
		MIDI_value = 69.0 + 12.0 * log2(m_glogs_f0[i + m_glogs_init_f0s]/440.0);
		m_glogs_f0_preference_weights[i] = 1.0/sqrt(2.0*M_PI*m_f0_params.prefer_stdev*m_f0_params.prefer_stdev)*exp(-(MIDI_value-m_f0_params.prefer_mean)*(MIDI_value-m_f0_params.prefer_mean)/(2.0*m_f0_params.prefer_stdev*m_f0_params.prefer_stdev));
		m_glogs_f0_preference_weights[i] = (0.01 + m_glogs_f0_preference_weights[i]) / (1.01);
		
		m_glogs_median_correction[i] = m_glogs_params.median_poly_coefs[0]*(i+1.0)*(i+1.0) + m_glogs_params.median_poly_coefs[1]*(i+1.0) + m_glogs_params.median_poly_coefs[2];
		m_glogs_sigma_correction[i] = 1.0 / (m_glogs_params.sigma_poly_coefs[0]*(i+1.0)*(i+1.0) + m_glogs_params.sigma_poly_coefs[1]*(i+1.0) + m_glogs_params.sigma_poly_coefs[2]);
	}
	
	double smooth_width = 1000.0; // hertz.
	double smooth_aux = (double)(m_warp_params.nsamps_twarp/2+1)*(m_fmax-smooth_width)/m_fmax;
	for (size_t i = 0; i < m_warp_params.nsamps_twarp/2+1; i++) {
		if (i <  smooth_aux) {
			m_glogs_hf_smoothing_window[i] = 1.0;
		} else {
			m_glogs_hf_smoothing_window[i] = ((double)i - (double)m_warp_params.nsamps_twarp/2.0)*(-1.0/((double)(m_warp_params.nsamps_twarp/2+1)-smooth_aux));
		}
	}
}

void
FChTransformF0gram::design_FFT() {
    in = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * m_nfft);
    out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * m_nfft);
	//TODO verificar que el tipo de datos de in_window es del tipo double, era float.
    in_window = (double*) fftw_malloc(sizeof (double) * m_nfft);
    planFFT = fftw_plan_dft_1d(m_nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	//TODO hacer diseño del FFT para el filtrado pasabajos.

}

void
FChTransformF0gram::design_FChT() {

    /*
     * FILES FOR DEBUGGING
     */

    //ofstream output("output.txt");


    /*  =============  WARPING DESIGN   ============= */

    // sampling frequency after oversampling
    m_warpings.fs_orig = m_warp_params.fact_over_samp * m_fs;

    // number of samples of the original signal frame
    m_warpings.nsamps_torig = 4 * m_warp_params.fact_over_samp * m_warp_params.nsamps_twarp;
    // equivalent to: m_warpings.nsamps_torig = m_warp_params.fact_over_samp * m_blockSize;

    // time instants of the original signal frame
    double t_orig[m_warpings.nsamps_torig];
    //float * t_orig = new float [m_warpings.nsamps_torig];
    for (size_t ind = 0; ind < m_warpings.nsamps_torig; ind++) {
        t_orig[ind] = ((double)(ind + 1) - (double)m_warpings.nsamps_torig / 2.0) / m_warpings.fs_orig;
    }

    // linear chirps warping definition as relative frequency deviation
    //double * freq_relative = new double [m_warpings.nsamps_torig * m_warp_params.num_warps];
	//TODO
	double *freq_relative = new double [m_warpings.nsamps_torig * m_warp_params.num_warps];
    define_warps_linear_chirps(freq_relative, t_orig);

    // maximum relative frequency deviation
    double freq_relative_max = 0;
    for (size_t i = 0; i < m_warpings.nsamps_torig; i++)
        for (size_t j = 0; j < m_warp_params.num_warps; j++)
            if (freq_relative_max < freq_relative[j * m_warpings.nsamps_torig + i])
                freq_relative_max = freq_relative[j * m_warpings.nsamps_torig + i];

    // sampling frequency of warped signal to be free of aliasing up to fmax
    m_warpings.fs_warp = 2 * m_fmax * freq_relative_max;

    // time instants of the warped signal frame
    double t_warp[m_warp_params.nsamps_twarp];
    for (size_t ind = 0; ind < m_warp_params.nsamps_twarp; ind++) {
        t_warp[ind] = ((double)((int)(ind + 1)- (int)m_warp_params.nsamps_twarp / 2)) / (double)m_warpings.fs_warp;
    }

    // design of warpings for efficient interpolation
    design_warps(freq_relative, t_orig, t_warp);


    /*
     * FILES FOR DEBUGGING
     */

    /*
    output << "chirp_rates" << endl;
    for (size_t j = 0; j < m_warp_params.num_warps; j++){
        output << m_warpings.chirp_rates[j];
        output << " ";
    }
    output << endl << "freq_relative" << endl;

    for (size_t i = 0; i < m_warpings.nsamps_torig; i++){
            for (size_t j = 0; j < m_warp_params.num_warps; j++){
                   output << freq_relative[j * m_warpings.nsamps_torig + i];
                   output << " ";
            }
            output << endl;
    }

    output << endl << "t_orig" << endl;

    for (size_t i = 0; i < m_warpings.nsamps_torig; i++){
                  output << t_orig[i] << endl ;
    }
     */

	delete [] freq_relative;
    //output.close();

    /*  =============  FFTW PLAN DESIGN   ============= */
	// Initialize 2-d array for warped signals
	x_warping = new double[m_warp_params.nsamps_twarp];
	m_absFanChirpTransform = (double*)fftw_malloc(sizeof (double) * m_warp_params.num_warps * (m_warp_params.nsamps_twarp/2 + 1));
	m_auxFanChirpTransform = (fftw_complex*)fftw_malloc(sizeof ( fftw_complex) * (m_warp_params.nsamps_twarp/2 + 1));
	plan_forward_xwarping = fftw_plan_dft_r2c_1d(m_warp_params.nsamps_twarp, x_warping, m_auxFanChirpTransform, FFTW_ESTIMATE);

}

void
FChTransformF0gram::design_warps(double * freq_relative, double * t_orig, double * t_warp) {
    /* the warping is done by interpolating the original signal in time instants
       given by the desired frequency deviation, to do this, the interpolation
       instants are stored in a structure as an integer index and a fractional value
       hypothesis: sampling frequency at the central point equals the original
     */

    m_warpings.pos_int = new size_t[m_warp_params.num_warps * m_warp_params.nsamps_twarp];
	m_warpings.pos_frac = new double[m_warp_params.num_warps * m_warp_params.nsamps_twarp];

	// vector of phase values
	double *phi = new double[m_warpings.nsamps_torig];
	double aux;

	// warped positions
	double *pos1 = new double[m_warp_params.nsamps_twarp*m_warp_params.num_warps];
	
    for (size_t i = 0; i < m_warp_params.num_warps; i++) {
        // vector of phase values
        //    float * phi;
        // integration of relative frequency to obtain phase values
        // phi = cumtrapz(t_orig,freq_relative(:,i)');
        // centering of phase values to force original frequency in the middle
        //phi = phi - phi(end/2);
        // interpolation of phase values to obtain warped positions
        //pos1(i,:) = interp1(phi,t_orig,t_warp)*fs_orig + length(t_orig)/2;
		
		// integration of relative frequency to obtain phase values
		cumtrapz(t_orig, freq_relative + i*(m_warpings.nsamps_torig), m_warpings.nsamps_torig, phi);

		// centering of phase values to force original frequency in the middle
		aux = phi[m_warpings.nsamps_torig/2];
		for (size_t j = 0; j < m_warpings.nsamps_torig; j++) {
			phi[j] -= aux;
		} //for

		// interpolation of phase values to obtain warped positions
		interp1(phi, t_orig, m_warpings.nsamps_torig, t_warp, pos1 + i*m_warp_params.nsamps_twarp, m_warp_params.nsamps_twarp);
		
    }

    // % previous sample index
    // pos1_int = uint32(floor(pos1))';
    // % integer corresponding to previous sample index in "c"
    // warps.pos1_int = (pos1_int - uint32(1));
    // % fractional value that defines the warped position
    // warps.pos1_frac = (double(pos1)' - double(pos1_int));

	// m_warpings.pos_int = new size_t[m_warp_params.num_warps * m_warp_params.nsamps_twarp];
	for (size_t j = 0; j < m_warp_params.nsamps_twarp*m_warp_params.num_warps; j++) {
		// previous sample index
		pos1[j] = pos1[j]*m_warpings.fs_orig + m_warpings.nsamps_torig/2 + 1;
		m_warpings.pos_int[j] = (size_t) pos1[j];
		m_warpings.pos_frac[j] = pos1[j] - (double)(m_warpings.pos_int[j]);
	} //for

	delete [] phi;
	delete [] pos1;
}

void
FChTransformF0gram::define_warps_linear_chirps(double * freq_relative, double * t_orig) {
    /**  define warps as relative frequency deviation from original frequency
                 t_orig : time vector
          freq_relative : relative frequency deviations
     */
    if (m_warp_params.alpha_dist == 0) {

        // linear alpha values spacing
        m_warpings.chirp_rates = new double [m_warp_params.num_warps];
        // WARNING m_warp_params.num_warps must be odd
        m_warpings.chirp_rates[0] = -m_warp_params.alpha_max;
        double increment = (double) m_warp_params.alpha_max / ((m_warp_params.num_warps - 1) / 2);

        for (size_t ind = 1; ind < m_warp_params.num_warps; ind++) {
            m_warpings.chirp_rates[ind] = m_warpings.chirp_rates[ind - 1] + increment;
        }
        // force zero value
        m_warpings.chirp_rates[(int) ((m_warp_params.num_warps - 1) / 2)] = 0;

    } else {
        // log alpha values spacing
        m_warpings.chirp_rates = new double [m_warp_params.num_warps];

        // force zero value
        int middle_point = (int) ((m_warp_params.num_warps - 1) / 2);
        m_warpings.chirp_rates[middle_point] = 0;

        double logMax = log10(m_warp_params.alpha_max + 1);
        double increment = logMax / ((m_warp_params.num_warps - 1) / 2.0f);
        double exponent = 0;

        // fill positive values
        int ind_log = middle_point;
        for (size_t ind = 0; ind < (m_warp_params.num_warps + 1) / 2; ind++) {
            m_warpings.chirp_rates[ind_log] = pow(10, exponent) - 1;
            exponent += increment;
            ind_log++;
        }
        // fill negative values
        for (size_t ind = 0; ind < (m_warp_params.num_warps - 1) / 2; ind++) {
            m_warpings.chirp_rates[ind] = -m_warpings.chirp_rates[m_warp_params.num_warps - 1 - ind];
        }
    }

    // compute relative frequency deviation
    for (size_t i = 0; i < m_warpings.nsamps_torig; i++)
        for (size_t j = 0; j < m_warp_params.num_warps; j++)
            freq_relative[j * m_warpings.nsamps_torig + i] = 1.0 + t_orig[i] * m_warpings.chirp_rates[j];
    //freq_relative[i * m_warpings.nsamps_torig + j] = 1.0 + t_orig[i] * m_warpings.chirp_rates[j];
    //freq_relative[i][j] = 1.0 + t_orig[i] * m_warpings.chirp_rates[j];
}

void
FChTransformF0gram::design_LPF() {

    //    in = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * tamanoVentana);
    //    out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * tamanoVentana);
    //    in_window = (float*) fftw_malloc(sizeof (float) * tamanoVentana);
    //    p = fftw_plan_dft_1d(tamanoVentana, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    double *lp_LPFWindow_aux = new double[m_blockSize/2+1];
    mp_LPFWindow = new double[m_blockSize/2+1];
    
    size_t i_max = (size_t) ((2.0*m_fmax/m_fs) * ( (double)m_blockSize / 2.0 + 1.0 ));
    for (size_t i = 0; i < m_blockSize/2+1; i++) {
        if (i >= i_max) {
            lp_LPFWindow_aux[i] = 0.0;
        } else {
            lp_LPFWindow_aux[i] = 1.0;
        }
    }
    LPF_time = (double*)fftw_malloc(sizeof ( double) * m_warpings.nsamps_torig);
		//memset((char*)LPF_time, 0, m_warpings.nsamps_torig * sizeof(double));	
			// sustituyo el memset por un for:
			for (size_t i = 0; i < m_warpings.nsamps_torig; i++) {
				LPF_time[i] = 0.0;
			}
		#ifdef DEBUG
			printf("	Corrio primer memset...\n");
		#endif
    LPF_frequency = (fftw_complex*)fftw_malloc(sizeof ( fftw_complex) * (m_warpings.nsamps_torig/2 + 1)); //tamaño de la fft cuando la entrada es real
		//memset((char*)LPF_frequency, 0, sizeof(fftw_complex) * (m_warpings.nsamps_torig/2 + 1));
			// sustituyo el memset por un for:
			for (size_t i = 0; i < (m_warpings.nsamps_torig/2 + 1); i++) {
				LPF_frequency[i][0] = 0.0;
				LPF_frequency[i][1] = 0.0;
			}
//	for (int i=0; i<(m_blockSize/2+1); i++) {
//		LPF_frequency[i] =  new fftw_complex;
//	}
    plan_forward_LPF = fftw_plan_dft_r2c_1d(m_blockSize, LPF_time, LPF_frequency, FFTW_ESTIMATE);
    plan_backward_LPF = fftw_plan_dft_c2r_1d(m_warpings.nsamps_torig, LPF_frequency, LPF_time, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
    
	size_t winWidth = 11;
    double *lp_hanningWindow = new double[winWidth]; 
    double accum=0;
    for (size_t i = 0; i < winWidth; i++) {
        lp_hanningWindow[i]=0.5*(1.0-cos(2*M_PI*(double)(i+1)/((double)winWidth+1.0)));
        accum+=lp_hanningWindow[i];
        
    }
    for (size_t i = 0; i < winWidth; i++) { //window normalization
        lp_hanningWindow[i]=lp_hanningWindow[i]/accum;
    }
    for (size_t i = 0; i < m_blockSize/2+1; i++) {
        //if (((i-(winWidth-1)/2)<0)||(i+(winWidth-1))/2>m_blockSize/2-1) {//consideramos winWidth impar, si la ventana sale del arreglo se rellena con el valor origianl
		if ( (i > (i_max + (winWidth-1)/2)) ||  (i <= (i_max - (winWidth-1)/2)) ) {
            mp_LPFWindow[i]=lp_LPFWindow_aux[i];
        } else {
            accum=0;
        	for (size_t j = -((winWidth-1)/2); j <= (winWidth-1)/2; j++) {
            	accum+=lp_LPFWindow_aux[i-j]*lp_hanningWindow[j+(winWidth-1)/2];
        	}
            mp_LPFWindow[i]=accum;
        }
    }

    delete[] lp_LPFWindow_aux;
    delete[] lp_hanningWindow;
}

void FChTransformF0gram::apply_LPF() {
    fftw_execute(plan_forward_LPF);
    for (size_t i = 0; i < m_blockSize/2+1; i++) {
        LPF_frequency[i][0]*=mp_LPFWindow[i];
        LPF_frequency[i][1]*=mp_LPFWindow[i];
    }
    fftw_execute(plan_backward_LPF);

	// TODO ver si hay que hacer fftshift para corregir la fase respecto al centro del frame.
	// nota: además de aplicar el LPF, esta función resamplea la señal original.
}

void FChTransformF0gram::clean_LPF() {
    delete[] mp_LPFWindow;

	fftw_destroy_plan(plan_forward_LPF);
	fftw_destroy_plan(plan_backward_LPF);
	fftw_free(LPF_time);
	fftw_free(LPF_frequency);
}

void FChTransformF0gram::reset() {

    // Clear buffers, reset stored values, etc

	delete [] m_warpings.pos_int;
	delete [] m_warpings.pos_frac;

    fftw_destroy_plan(planFFT);
	fftw_free(in); 
	fftw_free(out);

	clean_LPF();

	delete [] m_timeWindow;

	delete [] mp_HanningWindow;

	// Warping
	delete [] x_warping;
	fftw_destroy_plan(plan_forward_xwarping);
	fftw_free(m_absFanChirpTransform); 
	fftw_free(m_auxFanChirpTransform);

	// design_GLogS
	delete [] m_glogs_f0;
	delete [] m_glogs;
	delete [] m_glogs_n;
	delete [] m_glogs_index;
	delete [] m_glogs_posint;
	delete [] m_glogs_posfrac;
	delete [] m_glogs_third_harmonic_posint;
	delete [] m_glogs_third_harmonic_posfrac;
	delete [] m_glogs_third_harmonic;
	delete [] m_glogs_fifth_harmonic_posint;
	delete [] m_glogs_fifth_harmonic_posfrac;
	delete [] m_glogs_fifth_harmonic;
	delete [] m_glogs_f0_preference_weights;
	delete [] m_glogs_median_correction;
	delete [] m_glogs_sigma_correction;
	delete [] m_glogs_hf_smoothing_window;

}

FChTransformF0gram::FeatureSet
FChTransformF0gram::process(const float *const *inputBuffers, Vamp::RealTime) {

    //    // Do actual work!
    //

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
		printf("\n	----- DEBUG INFORMATION ----- \n");
		printf("	m_fs = %f Hz.\n",m_fs);
		printf("	fs_orig = %f Hz.\n",m_warpings.fs_orig);
		printf("	fs_warp = %f Hz.\n",m_warpings.fs_warp);
		printf("	m_nfft = %d.\n",m_nfft);
		printf("	m_blockSize = %d.\n",m_blockSize);
		printf("	m_warpings.nsamps_torig = %d.\n",m_warpings.nsamps_torig);
		printf("	m_warp_params.num_warps = %d.\n",m_warp_params.num_warps);
		printf("	m_glogs_harmonic_count = %d.\n",m_glogs_harmonic_count);
	#endif

    // size_t n = m_nfft/2 + 1;
	// double *tbuf = in_window;

	for (size_t i = 0; i < m_blockSize; i++) {
        LPF_time[i] = (double)(inputBuffers[0][i]) * m_timeWindow[i];
    }

//	#ifdef DEBUG
//		printf("	HASTA ACÁ ANDA!!!\n");
//		cout << flush;
//	#endif

	apply_LPF();
	// Señal filtrada queda en LPF_time

	Feature feature;
    feature.hasTimestamp = false;


	/* Solo a modo de prueba, voy a poner la salida del filtrado en «in» y
	voy a mostrar la FFT de eso, para ver el efecto del filtrado. */
//    for (size_t i = 0; i < m_nfft; i++) {
//        in[i][0] = tbuf[i];
//        in[i][1] = 0;
//    }
//	fftw_execute(planFFT);
//	double real, imag;
//	for (size_t i=0; i<n; ++i) {		// preincremento?? ver version de nacho
//		real = out[i][0];
//		imag = out[i][1];
//		feature.values.push_back(real*real + imag*imag);
//	}
//	fs[0].push_back(feature);

//	float real; 
//	float imag;
//	for (size_t i=0; i<m_blockSize/2+1; i++) {
//		real = (float)(LPF_frequency[i][0]);
//		imag = (float)(LPF_frequency[i][1]);
//		feature.values.push_back(real*real+imag*imag);
//		//feature.values.push_back((float)(mp_LPFWindow[i]));
//	}

// ----------------------------------------------------------------------------------------------
// 		Hanning window & FFT for all warp directions

	double max_glogs = -DBL_MAX;
	size_t ind_max_glogs = 0;

	for (size_t i_warp = 0; i_warp < m_warp_params.num_warps; i_warp++) {
		// Interpolate
		interp1q(LPF_time, (m_warpings.pos_int) + i_warp*m_warp_params.nsamps_twarp, m_warpings.pos_frac + i_warp*m_warp_params.nsamps_twarp, x_warping, m_warp_params.nsamps_twarp);

		// Apply window
		for (size_t i = 0; i < m_warp_params.nsamps_twarp; i++) {
			x_warping[i] *= mp_HanningWindow[i];
		}

		// Transform
		fftw_execute(plan_forward_xwarping);

		// Copy result
		//memcpy(m_absFanChirpTransform + i_warp*(m_warp_params.nsamps_twarp/2+1), m_auxFanChirpTransform, (m_warp_params.nsamps_twarp/2+1)*sizeof(fftw_complex)); asi como esta no funciona
		double *aux_abs_fcht = m_absFanChirpTransform + i_warp*(m_warp_params.nsamps_twarp/2+1);
		for (size_t i = 0; i < (m_warp_params.nsamps_twarp/2+1); i++) {
			aux_abs_fcht[i] = log10(1.0 + 10.0*sqrt(m_auxFanChirpTransform[i][0]*m_auxFanChirpTransform[i][0]+m_auxFanChirpTransform[i][1]*m_auxFanChirpTransform[i][1]));
			// smoothing high frequency values
			//aux_abs_fcht[i] *= m_glogs_hf_smoothing_window[i];
		}
		
//      -----------------------------------------------------------------------------------------
// 		GLogS
		interp1q(aux_abs_fcht, m_glogs_posint, m_glogs_posfrac, m_glogs_interp, m_glogs_harmonic_count);
		size_t glogs_ind = 0;
		for (size_t i = 0; i < m_glogs_num_f0s; i++) {
			double glogs_accum = 0;
			for (size_t j = 1; j <= m_glogs_n[i]; j++) {
				glogs_accum += m_glogs_interp[glogs_ind++];
			}
			m_glogs[i + i_warp*m_glogs_num_f0s] = glogs_accum/(double)m_glogs_n[i];
		}

//		Sub/super harmonic correction
		interp1q(m_glogs + i_warp*m_glogs_num_f0s, m_glogs_third_harmonic_posint, m_glogs_third_harmonic_posfrac, m_glogs_third_harmonic, (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
		interp1q(m_glogs + i_warp*m_glogs_num_f0s, m_glogs_fifth_harmonic_posint, m_glogs_fifth_harmonic_posfrac, m_glogs_fifth_harmonic, (m_f0_params.num_octs+1)*m_f0_params.num_f0s_per_oct);
		for (size_t i = m_glogs_num_f0s-1; i >= m_glogs_init_f0s; i--) {
			m_glogs[i + i_warp*m_glogs_num_f0s] -= MAX(MAX(m_glogs[i-m_f0_params.num_f0s_per_oct + i_warp*m_glogs_num_f0s],m_glogs_third_harmonic[i-m_glogs_init_f0s]),m_glogs_fifth_harmonic[i-m_glogs_init_f0s]);
			//m_glogs[i] -= MAX(m_glogs[i-m_f0_params.num_f0s_per_oct],m_glogs_third_harmonic[i-m_glogs_init_f0s]);
		}
		for (size_t i = m_glogs_init_f0s; i < m_glogs_num_f0s-m_f0_params.num_f0s_per_oct; i++) {
			m_glogs[i + i_warp*m_glogs_num_f0s] -= 0.3*m_glogs[i+m_f0_params.num_f0s_per_oct + i_warp*m_glogs_num_f0s];
			// Median, sigma $ weights correction
			m_glogs[i + i_warp*m_glogs_num_f0s] = (m_glogs[i + i_warp*m_glogs_num_f0s]-m_glogs_median_correction[i-m_glogs_init_f0s])*m_glogs_sigma_correction[i-m_glogs_init_f0s]*m_glogs_f0_preference_weights[i-m_glogs_init_f0s];
		}
	
		// Look for maximum value to determine best direction
		for (size_t i = m_glogs_init_f0s; i < m_glogs_num_f0s-m_f0_params.num_f0s_per_oct; i++) {
			if (m_glogs[i + i_warp*m_glogs_num_f0s] > max_glogs) {
				max_glogs = m_glogs[i + i_warp*m_glogs_num_f0s];
				ind_max_glogs = i_warp;
			}
		}
	}	
	
// ----------------------------------------------------------------------------------------------

	for (size_t i=m_glogs_init_f0s; i< m_glogs_num_f0s - m_f0_params.num_f0s_per_oct; i++) {
	//for (size_t i=0; i<(m_warp_params.nsamps_twarp/2+1); i++) {
		//feature.values.push_back((float)(m_warpings.pos_int[i])+ (float)(m_warpings.pos_frac[i]));
		//feature.values.push_back((float)(phi[i]*100000.0));
		//feature.values.push_back((float)(t_orig[i]));
		//feature.values.push_back((float)(pos1[i]));
		//feature.values.push_back((float)x_warping[i]);
		//feature.values.push_back(m_absFanChirpTransform[i + ind_max_glogs*(m_warp_params.nsamps_twarp/2+1)]);
		//feature.values.push_back((float)m_glogs[i+(long)ind_max_glogs*(long)m_glogs_num_f0s]);
		switch (m_f0gram_mode) {
			case 1:		
					max_glogs = -DBL_MAX;
					for	(size_t i_warp = 0; i_warp < m_warp_params.num_warps; i_warp++) {
						if (m_glogs[i + i_warp*m_glogs_num_f0s] > max_glogs) {
							max_glogs = m_glogs[i + i_warp*m_glogs_num_f0s];
							ind_max_glogs = i_warp;
						}
					}
					feature.values.push_back((float)max_glogs);
					break;
			case 0:
					feature.values.push_back((float)m_glogs[i+(size_t)ind_max_glogs*(size_t)m_glogs_num_f0s]);
					break;
		}
		//feature.values.push_back((float)m_glogs_hf_smoothing_window[i]);
	}

// ----------------------------------------------------------------------------------------------

	fs[0].push_back(feature);

	#ifdef DEBUG
		printf("	----------------------------- \n");
	#endif
	
	return fs;
//---------------------------------------------------------------------------

    //return FeatureSet();
}

FChTransformF0gram::FeatureSet
FChTransformF0gram::getRemainingFeatures() {
    return FeatureSet();
}

void
FChTransformF0gram::design_time_window() {

	size_t transitionWidth = (size_t)m_blockSize/128 + 1;;
    m_timeWindow = new double[m_blockSize];
	double *lp_transitionWindow = new double[transitionWidth];

	//memset(m_timeWindow, 1.0, m_blockSize);
	for (size_t i = 0; i < m_blockSize; i++) {
		m_timeWindow[i] = 1.0;
	}

	for (size_t i = 0; i < transitionWidth; i++) {
        lp_transitionWindow[i]=0.5*(1.0-cos(2*M_PI*(double)(i+1)/((double)transitionWidth+1.0)));
    }

	for (size_t i = 0; i < transitionWidth/2; i++) {
		m_timeWindow[i] = lp_transitionWindow[i];
		m_timeWindow[m_blockSize-1-i] = lp_transitionWindow[transitionWidth-1-i];
	}

	#ifdef DEBUG
		for (int i = 0; i < m_blockSize; i++) {
			if ((i<transitionWidth)) {
				printf("	m_timeWindow[%d] = %f.\n",i,m_timeWindow[i]);
			}
		}
	#endif

	delete [] lp_transitionWindow;
}

