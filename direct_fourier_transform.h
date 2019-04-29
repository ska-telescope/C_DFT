
// Copyright 2019 Adam Campbell, Seth Hall, Andrew Ensor
// Copyright 2019 High Performance Computing Research Laboratory, Auckland University of Technology (AUT)

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DIRECT_FOURIER_TRANSFORM_H_
#define DIRECT_FOURIER_TRANSFORM_H_

//=========================//
//   Algorithm Constants   //
//=========================//

// Pi (double precision)
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

// Speed of light
#ifndef C
	#define C 299792458.0
#endif

//=========================//
//        Structures       //
//=========================//

typedef struct Complex {
	double real;
	double imaginary;
} Complex;

typedef struct Source {
	double l;
	double m;
	double intensity;
} Source;

typedef struct Visibility {
	double u;
	double v;
	double w;
	Complex brightness;
	double intensity;
} Visibility;

typedef struct Config {
	int num_visibilities;
	int num_sources;
	char *source_file;
	char *vis_file;
	bool synthetic_sources;
	bool synthetic_visibilities;
	bool gaussian_distribution_sources;
	bool forceZeroWTerm;
	double min_u;
	double max_u;
	double min_v;
	double max_v;
	double min_w;
	double max_w;
	double grid_size;
	double cell_size;
	double uv_scale;
	double frequency_hz;
} Config;


//=========================//
//     Function Headers    //
//=========================//
void init_config (Config *config);

void load_sources(Config *config, Source **sources);

void load_visibilities(Config *config, Visibility **visibilities);

void extract_visibilities(Config *config, Source *sources, Visibility *visibilities, int num_visibilities);

void save_visibilities(Config *config, Visibility *visibilities);

double random_in_range(double min, double max);

double generate_sample_normal(void);

void unit_test_init_config(Config *config);

double unit_test_generate_approximate_visibilities(void);

#endif /* DIRECT_FOURIER_TRANSFORM_H_ */

#ifdef __cplusplus
}
#endif

