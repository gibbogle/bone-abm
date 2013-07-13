#include <qstring.h>
#include "params.h"

Params::Params()
{
	PARAM_SET params[] = {

{"MONOCYTE_DIAMETER", 10.0, 5.0, 25.0,
"Diameter of a monocyte",
"Currently osteocyte-precursor monocytes are assumed to be all the same size.\n\
[um]"},

{"MOTILITY_BETA", 0.2, 0.25, 0.5,
"Motility speed parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1.\n\
Median T cell speed is roughly proportional to MOTILITY_BETA."},

{"MOTILITY_RHO", 0.76, 0.5, 0.9,
"Motility persistence parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1.\n\
MOTILITY_RHO determines the extent to which motion is in the same direction from one time step to the next."},

{"CXCL12_CHEMOLEVEL", 1.0, 0.0, 1.0,
"Level of CXCL12 chemotaxis",
"Level of CXCL12-mediated chemotactic influence on a monocyte.  The overall CXCL12 chemotaxis influence factor is given by:\n\
CXCL12 influence factor = ..."},

{"CXCL12_KDIFFUSION", 100, 0.0, 0.0,
"CXCL12 diffusion coefficient",
"Coefficient of diffusion of CXCL12 through marrow."},

{"CXCL12_HALFLIFE", 30, 0.0, 0.0,
"CXCL12 half-life",
"Half-life of CXCL12.\n\
[min]"},

{"S1P_CHEMOLEVEL", 0.1, 0.0, 1.0,
"Level of S1P chemotaxis",
"Level of S1P-mediated chemotactic influence on a monocyte.  The overall S1P chemotaxis influence factor is given by:\n\
S1P influence factor = (Level of S1P chemotaxis)*min(1.0,(S1P gradient)/(Reference S1P gradient))*(cell S1P1 level)"},

{"S1P_KDIFFUSION", 100, 0.0, 0.0,
"S1P diffusion coefficient",
"Coefficient of diffusion of S1P through marrow."},

{"S1P_HALFLIFE", 30, 0.0, 0.0,
"S1P half-life",
"Half-life of S1P.\n\
[min]"},

{"S1P_GRADLIM", 0.0001, 0.0, 0.0,
"Reference S1P gradient",
"S1P gradient is scaled by this to arrive at the gradient influence factor.  The overall S1P chemotaxis influence factor is given by:\n\
S1P influence factor = (Level of S1P chemotaxis)*min(1.0,(S1P gradient)/(Reference S1P gradient))*(cell S1P1 level)"},

{"S1P1_THRESHOLD", 0.4, 0.0, 0.0,
"S1P1 threshold for egress",
"Monocyte egress is possible when the cell's S1P1 level exceeds this threshold.  The probability of egress (in a time step) of a monocyte in contact with a capillary is then given by:\n\
probability = (Max capillary egress probability)*(S1P1 - S1P1_THRESHOLD)/(1 - S1P1_THRESHOLD)\n\
[0-1]"},

{"S1P1_RISETIME", 12, 0.0, 0.0,
"S1P1 rise time",
"Time required for monocyte S1P1 level to rise from 0 to 1. (Currently the rate of increase is treated as constant).\n\
[hours]"},

{"OB_PER_MM3", 5, 0, 0,
"Osteoblasts/mm3",
"Average number of osteoblasts per cubic mm of marrow (initial)."},

{"OB_SIGNAL_FACTOR", 10, 0, 0,
"OB signal factor",
"Ratio of CXCL12 secretion to bone signal strength."},

{"X_SIZE", 500, 0, 0,
"Bone patch size",
"Dimension of the modelled bone region (in X and Z directions).\n\
[um]"},

{"Y_SIZE", 200, 0, 0,
"Slice thickness",
"Dimension of the modelled region in Y direction.\n\
[um]"},

{"CAPILLARY_DIAMETER", 30.0, 20.0, 200.0,
"Diameter of a capillary",
"(Currently we consider only a single capillary with a fixed diameter)\n\
[um]"},

{"MONO_PER_MM3", 6000, 0, 0,
"Monocytes/mm3",
"Average number of monocytes per cubic mm of marrow (initial)."},

{"IN_PER_HOUR", 5, 0, 0,
"Monocyte influx",
"Rate of influx pf osteoclast-precursor monocytes from the blood.\n\
[/hour]"},

{"STEM_PER_MM3", 20, 0, 0,
"Stem cells/mm3",
"Average number of stem cells per cubic mm of marrow."},

{"STEM_CYCLETIME", 6.0, 3.0, 12.0,
"Stem cell cycle time",
"Time taken for stem cell division.\n\
[hours]"},

{"CROSSING_TIME", 100.0, 5.0, 500.0,
"Capillary crossing time",
"Time taken for a monocyte the cross into a capillary\n\
[mins]"},

{"FUSING_TIME", 120.0, 10.0, 500.0,
"Monocyte -> osteoclast fusing time",
"Time taken for a group of monocytes to form an osteoclast, after they have aggregated\n\
[mins]"},

{"CLAST_LIFETIME", 12.0, 0, 0,
"Osteoclast lifetime",
"Length of time that an osteoclast remains active\n\
[days]"},

{"CLAST_DWELL_TIME", 360.0, 0.0, 0.0,
"Osteoclast dwell time",
"Minimum length of time that an osteoclast remains in one spot, before moving one grid\n\
[mins]"},

{"MAX_RESORPTION_RATE", 400, 0.0, 0.0,
"Maximum resorption rate",
"Maximum rate of bone removal by an osteoclast\n\
[um3/min]"},

{"MAX_RESORPTION_D", 50, 0, 0,
"Maximum resorption D",
"Excavation depth at which rate of bone removal goes to zero\n\
[um]"},

{"MAX_RESORPTION_N", 10, 0, 0,
"Maximum resorption N",
"Number of osteoclast monocytes corresponding to maximum resorption rate."},

{"CROSS_PROB", 0.5, 0.0, 1.0,
"Max capillary egress probability",
"The probability (per time step) of a monocyte beginning to cross the endothelium into the capillary, while next to a capillary site\n\
is given by: probability = (Max capillary egress probability)*(S1P1 - S1P1_THRESHOLD)/(1 - S1P1_THRESHOLD)\n\
[0-1]"},

{"NDAYS", 100.0, 0, 0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 3, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 20, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"RUNCASE", 0, 0, 1,
"Run case",
"Monocyte or osteoclast model run"}

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}
