// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

void* GA(					// run the GA.  returns member.
        int cbMemberArg,			// size of population
        void (*pfnGenerateRandom)(void* pv),	// randomly generate a member
        void (*pfnMutateRandom)(void* pv, long cIter),	// mutate a member
        void (*pfnTweak)(void* pv),			// tweak a member
        double (*pfnComputeSuitability)(void* pv),	// fitness function
        double zSuitabilityMax,				// perfect fit
	int cBestArg,				// # of members to keep per gen.
	double tMaxSec				// timeout
        );
