#include <stdio.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void
DE (
	int maxGenerations,
	int populationSize,
	double CR,
	int dimensions,
	double F,
	gsl_rng * r,
	double (* evaluate) (double [])
	);

