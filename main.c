#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "utilities.h"
#include "networkDynamics.h"
#include "discreteFunction.h"

#include "de.h"
//#include "evaluatorMC_Glass.h"
#include "evaluatorHINTS_Glass.h"

#include <time.h>

gsl_rng * r;

// Shame on me. I have to resort to this quick and dirty solution...
// Global variables... eough! =(
// double * expCa, * expV, * costs;
double * costs;
gsl_matrix * solutions, * experiments;
//gsl_vector * costs;// * expCa, * expV;
int steps;


int main (int argc, const char * argv[]) {
	if (argc < 4) {
		printf ("\n\
###############################################################\n \
		Too few parameters pased!\n\
		USAGE: %s maxIterations populationSize dataPoints\n\
###############################################################\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// See how many points do we want to calculate
	steps = atoi (argv[3]);
	int dims = 26;
	int popSize = atoi (argv[2]);
	int maxGens = atoi (argv[1]);
	
	// We do need to set up the generator in order to pass it initialized to the DE routine
	const gsl_rng_type *tipoRNG;
	gsl_rng_env_setup();
	
	tipoRNG = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(tipoRNG);
	
	// Seed the generator 
	long seed = time(NULL);
	gsl_rng_set(r, seed);
	printf("Seed = %ld\n", seed);
	
	// Initialize the global matrix and vectors
	solutions = gsl_matrix_alloc (popSize, dims);	
	costs = (double *) malloc (steps * sizeof (*costs));//gsl_vector_alloc (popSize);

        experiments = gsl_matrix_alloc (steps, 22);

	//expCa = (double *) malloc (steps * sizeof (*expCa)); //gsl_vector_alloc (steps);
	//expV = (double *) malloc (steps * sizeof (*expCa)); //gsl_vector_alloc (steps);

	//Use Differential Evolution Algorithm to renew
	// DE (maxGenerations, populationSize, CR, dimensions, F, RNG, fitnessFunction);
	// 26 dimensions. 22 binary thresholds, 4 ternary thresholds.
	DE (maxGens, popSize, 0.9, dims, 0.5, r, evaluator);

	// At this point, we DO have results. Let's report them
	{
		FILE * fThrs = fopen ("dataPSM/thresholds.dat", "w");
		if (!fThrs) {
			fprintf(stderr, "Couldn't write to thresholds file. Aborting\n");
		       	exit(EXIT_FAILURE);
		}
		for (int i = 0; i < popSize; i++) {
			for (int j = 0; j < dims; j++) {
				double value = gsl_matrix_get (solutions, i, j);
				if (j<22)
					value = fmod(value,1);
				else
					value = fmod(value,2);
				fprintf (fThrs, "%lf\t", value);
				printf ( "[%.2lf]\t", value);
			}
			fprintf (fThrs, "%lf\n", costs[i]); //gsl_vector_get (costs, i) );
			printf ( "Cost = %g\n", costs[i] );
		}
		fclose (fThrs);
	}

	printf ("\nFinal population cost stats:\n\nMean = %lf\nStdDev = %lf\n\n", 
			//gsl_stats_mean (costs->data, costs->stride, costs->size),
			gsl_stats_mean (costs, 1, popSize),
			gsl_stats_sd (costs, 1, popSize)
			//gsl_stats_sd (costs->data, costs->stride, costs->size),
			//gsl_vector_min (costs), gsl_vector_max (costs)
	       );

	// Patching the quick and dirty solution, let's free the
	// memory allocated in the global variables
	gsl_matrix_free (solutions);
	//gsl_vector_free (costs);
	free (costs);
        gsl_matrix_free(experiments);
        //free (expCa); free (expV);
	//gsl_vector_free (expCa);
	//gsl_vector_free (expV);

	gsl_rng_free(r);
	
	return 0;
}
