#ifndef RESULT_SET_H
#define RESULT_SET_H
struct result_set {
	QString casename;
	int nsteps;
 	double *tnow;      // simulation time (mins)
    double *nDC;       // number of DCs
	double *nborn;     // total number of monocytes created
	double *ntot;      // total monocyte population
    double *ncogseed;  // number of naive cognate T cells that have arrived
    double *ncog;      // current number of cognate T cells
    double *ndead;     // number of cognate T cells that have died
    double *teffgen;   // number of activated cognate T cells that have left the LN
 	double max_tnow;      // simulation time (mins)
    double max_nDC;       // number of DCs
	double max_nborn;     // total number of monocytes created
	double max_ntot;      // total monocyte population
    double max_ncogseed;  // number of naive cognate T cells that have arrived
    double max_ncog;      // current number of cognate T cells
    double max_ndead;     // number of cognate T cells that have died
    double max_teffgen;   // number of activated cognate T cells that have left the LN
};
typedef result_set RESULT_SET;
#endif
