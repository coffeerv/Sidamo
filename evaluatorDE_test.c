#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "de.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

double evaluator (double x[]);
//double correlation(const double *x1, const double *x2, const int n);

int
main (int argc, const char * argv [ ])
{
	if (argc < 2) {
		printf ("\nMust specify how many iterations the DE algorithm should run\n");
		exit(EXIT_FAILURE);
	}
	
	gsl_rng *ra;
	const gsl_rng_type *tipoRNG;
	gsl_rng_env_setup();
	tipoRNG = gsl_rng_ranlxs2;
	ra = gsl_rng_alloc(tipoRNG);
	/* Seed the generator */
	long seed = time(NULL);
	gsl_rng_set(ra,seed);



	DE (atoi(argv[1]),25,0.9,1,0.5,ra,evaluator);

	gsl_rng_free (ra);

	return EXIT_SUCCESS;
}

double evaluator (double x[]) 
{
	//int i;
 	//double res = 0.0;
 	//for (i = 0; i < 1; i++) res += (x[i] * x[i]); 
 	//return ( res );

	// Generate a time series with a known parameter
	double knownSinTS [100];
	for (int i = 0; i < 100; i++) {
		knownSinTS [i] = sin (3*i);
	}
	// Generate a time series with a random parameter
	double randomSinTS [100];
	for (int i = 0; i < 100; i++) {
		randomSinTS [i] = sin ( x[0] * i);
	}
	// Now see how correlated they are
	// Generate a time series
	// Correlation = -1 will yield a score of 2
	// Correlation =  0 will yield a score of 1
	// Correlation =  1 will yield a score of 0
	return 1 - gsl_stats_correlation (knownSinTS, 1, randomSinTS, 1, 100);
}
/*
double correlation(const double *x1, const double *x2, const int n){
	double corr = 0;
	for (int i = 0; i < n; i++) {
		corr += x1[i] * x2[i];
	}
	corr /= n;
	return corr;
}
*/

