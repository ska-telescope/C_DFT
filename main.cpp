
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

#include <cstdlib>
#include <cstdio>

#include "direct_fourier_transform.h"

int main(int argc, char **argv)
{
	printf("==============================================================\n");
	printf(">>> AUT HPC Research Laboratory - Direct Fourier Transform <<<\n");
	printf("==============================================================\n\n");

	// Initialise the DFT configuration
	Config config;
	init_config(&config);

	// Obtain Sources from file, or synthesize
	Source *sources = NULL;
	load_sources(&config, &sources);
	// Something went wrong during loading of sources
	if(sources == NULL)
		return EXIT_FAILURE;

	// Obtain Visibilities from file, or synthesize
	Visibility *visibilities = NULL;
	load_visibilities(&config, &visibilities);

	// Something went wrong during loading of visibilities
	if(visibilities == NULL)
	{
		if(sources) free(sources);
		return EXIT_FAILURE;
	}

	printf(">>> UPDATE: Performing extraction of visibilities from sources...\n\n");
	extract_visibilities(&config, sources, visibilities, config.num_visibilities);
	printf(">>> UPDATE: Visibility extraction complete...\n\n");
	
	// Save visibilities to file
	save_visibilities(&config, visibilities);

	// Clean up
	if(visibilities) free(visibilities);
	if(sources)      free(sources);

	printf(">>> UPDATE: Direct Fourier Transform operations complete, exiting...\n\n");
    
	return EXIT_SUCCESS;
}
