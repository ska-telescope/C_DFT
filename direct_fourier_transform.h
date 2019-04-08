
// Copyright 2019 Adam Campbell, Seth Hall, Andrew Ensor
// High Performance Computing Research Laboratory, 
// Auckland University of Technology (AUT)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

