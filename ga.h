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
