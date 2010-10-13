#include <qstring.h>
#include "params.h"

Params::Params()
{
	PARAM_SET params[] = {
{"TC_AVIDITY_MEDIAN", 1.0, 0.1, 10.0,
"Receptor avidity median parameter",
"TCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
(TCR stimulation rate is proportional to the product of TC avidity and DC antigen density.)"},

{"TC_AVIDITY_SHAPE", 1.1, 1.01, 3.0,
"Receptor avidity shape parameter",
"TCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
The shape value must be greater than 1, and values close to 1 give distributions that are close to normal."},

{"TC_STIM_RATE_CONSTANT", 1, 0.0, 0.0,
"Stimulation rate constant",
"Rate constant Ks for TCR stimulation, where:\n\
rate of TCR stimulation = Ks*(TCR avidity)*(DC antigen density)\n\
[molecules/min]"},
	
{"TC_STIM_HALFLIFE", 24.0, 0.0, 100.0,
"Stimulation halflife",
"Integrated TCR stimulation decays with a specified halflife. \n\
[hours]"},

{"DIVIDE1_MEDIAN", 6.0, 0.0, 0.0,
"1st division time median parameter",
"The time taken for the first T cell division, after full activation, has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE1_SHAPE", 1.1, 1.01, 3.0,
"1st division time shape parameter",
"The time taken for the first T cell division, after full activation, has a lognormal distribution, described by the median and shape parameters."},

{"MOTILITY_BETA", 0.35, 0.25, 0.5,
"Motility speed parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1. Median T cell speed is roughly proportional to MOTILITY_BETA."},

{"MOTILITY_RHO", 0.76, 0.5, 0.9,
"Motility persistence parameter",
"T cell motility is described by speed and persistence parameters, each in the range 0 - 1. MOTILITY_RHO determines the extent to which motion is in the same direction from one time step to the next."},

{"DC_LIFETIME_MEDIAN", 3.0, 1.0, 20.0,
"Lifetime median parameter",
"DC lifetime has a lognormal distribution, described by the median and shape parameters.\n\
[days]"},

{"DC_LIFETIME_SHAPE", 1.2, 1.01, 3.0,
"Lifetime shape parameter",
"DC lifetime has a lognormal distribution, described by the median and shape parameters."},

{"IL2_THRESHOLD", 150, 10.0, 500.0,
"Stimulation threshold for IL-2",
"Integrated TCR stimulation needed to initiate IL-2/CD5 production."},

{"ACTIVATION_THRESHOLD", 150, 10.0, 500.0,
"Stimulation threshold for activation",
"Integrated TCR stimulation level needed for full activation."},

{"FIRST_DIVISION_THRESHOLD", 300, 10.0, 1000.0,
"Stimulation threshold for first division",
"Integrated TCR stimulation level needed for first division."},

{"DIVISION_THRESHOLD", 80, 10.0, 1000.0,
"Stimulation threshold for subsequent divisions",
"Integrated TCR stimulation level needed for subsequent divisions."},

{"EXIT_THRESHOLD", 480, 10.0, 1000.0,
"Stimulation threshold for exit",
"Integrated TCR stimulation level below which exit is permitted (using Exit Rule #2)."},

{"STIMULATION_LIMIT", 1000, 0.0, 0.0,
"Maximum stimulation level",
"Maximum integrated TCR stimulation level (saturation level)."},

{"NX", 50, 30, 300,
"Lattice size",
"Dimension of the lattice (number of sites in X, Y and Z directions).  Typically 4*BLOB_RADIUS is OK."},

{"EXIT_RULE", 1, 1, 2,
"Exit rule",
"T cell exit rule.  1 = use NGEN_EXIT, 2 = use EXIT_THRESHOLD, 3 = no restriction."},

{"EXIT_REGION", 1, 1, 2,
"Exit region",
"Determines blob region for cell exits: 1 = everywhere, 2 = lower half of blob, 3 = by chemotaxis, via discrete exits."},

{"CROSS_PROB", 0.02, 0.0, 1.0,
"Capillary egress probability",
"The probability (per time step) of a monocyte beginning to cross the endothelium into the capillary, while next to a capillary site\n\
[0-1]"},

{"CHEMO_RADIUS", 30.0, 10.0, 200.0,
"Radius of chemotactic influence",
"Range of chemotactic influence of an exit site or DC on T cell motion.  At this distance the influence is reduced to 5% of its maximum value.\n\
[um]"},

{"CHEMO_K_EXIT", 0.5, 0.0, 1.0,
"Exit chemotaxis influence parameter",
"Strength of chemotactic influence on T cell motion towards exits."},

{"CHEMO_K_DC", 0.0, 0.0, 10.0,
"DC chemotaxis influence parameter",
"Strength of chemotactic influence on T cell motion towards DCs."},

{"NDAYS", 4.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 1, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 20, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"SPECIES", 1, 0, 1,
"Species",
"Animal species (mouse or human)"}

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
