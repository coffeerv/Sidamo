#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// This is the 21st century. Let's use C99
#include <stdbool.h>

#include "networkDynamics.h"
#include "discreteFunction.h"
#include "utilities.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

static const double condInicial [22] = {
		0.25,	//sp
		0.25,	//gc
		0.25,	//cGMP
		0.25,	//cKcGMP
		0.25,	//pK
		1.0,	//v
		0.25,	//NHE
		0.25,	//NCE
		0.25,	//HCN
		0.25,	//AC
		0.66,	//LVA
		0.25,	//uPH
		0.25,	//pNa
		0.25,	//cAMP
		0.66,	//HVA
		0.66,	//pCa
		0.25,	//CaP
		0.25,	//cClCa a.k.a CaCC
		0.25,	//cCacAMP a.k.a CatSper
		0.25,	//CKCa		¡ATENCIÓN! En esta representación, el nodo 26 es el 19
		0.25,	//pCl
		0.25,	//PDE		¡ATENCIÓN! En esta representación, el nodo 22 es el 21
	};
	
static const double alphas [22] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

double evaluator (double x[]);

