#include "de.h"

extern gsl_matrix * solutions;
//extern gsl_vector * costs;
extern double * costs;

void DE ( int maxGenerations, int populationSize, double CR, int dimensions, double F, gsl_rng * myRNG, double (* evaluate) (double []) ) {
	int count = 0;
	double score = 50000.0;

	double trial [dimensions], cost [populationSize];
	double x1 [populationSize][dimensions];
	double x2 [populationSize][dimensions];


	// Initialize population
	for (int i = 0; i < populationSize; i++) {
		for (int j = 0; j < dimensions; j++) {
			x1 [i][j] = gsl_rng_uniform(myRNG);
			x2 [i][j] = 0.0;
		}
		cost [i] = 50000.0;
	}

	// Main loop
	while (count < maxGenerations) {
		for (int i = 0; i < populationSize; i++) {
			int a, b, c;
			// Mutate and recombine
			do a = gsl_rng_uniform (myRNG) * populationSize; while (a == i);
			do b = gsl_rng_uniform (myRNG) * populationSize; while (b == i || b == a);
			do c = gsl_rng_uniform (myRNG) * populationSize; while (c == i || c == a || c == b);
			int j = gsl_rng_uniform(myRNG) * dimensions;
			
			for ( int k = 0; k <= dimensions; k++) {
				if ( gsl_rng_uniform (myRNG) < CR || k == dimensions ) {
					trial [j] = x1 [c][j] + F * (x1 [a][j] - x1 [b][j]);
				}
				else trial [j] = x1 [i][j];
				j = (j+1) % dimensions;
			}
			// Evaluate and Select
			score = evaluate (trial);
			if ( score <= cost [i] ) {
				for (int m = 0; m < dimensions; m++) x2 [i][m] = trial [m];
				cost [i] = score;
			} 
			else for (int m = 0; m < dimensions; m++) x2 [i][m] = x1 [i][m];
		}
		for(int i = 0; i < populationSize; i++) {
			for (int j = 0; j < dimensions; j++) x1 [i][j] = x2 [i][j];
		}
		count++;
	}
	//End of main loop

	// Now report the results...
	for (int i = 0; i < populationSize; i++) {
		for (int j = 0; j < dimensions; j++) {
			gsl_matrix_set (solutions, i, j, x1 [i][j]);
			//printf ( "[%.2lf]\t", x1 [i][j] );
		}
		//gsl_vector_set (costs, i, cost [i]);
		costs[i] = cost [i];
		//printf ( "Cost = %g\n", cost [i] );
	}
} // End of DE

