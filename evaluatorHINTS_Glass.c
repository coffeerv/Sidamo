#include "evaluatorHINTS_Glass.h"
#include <stdio.h>

extern gsl_matrix * experiments;
extern int steps;

double evaluator (double x[]) 
{
	// Once and for all code
	static bool isFirstTime = true;
	if (isFirstTime) {
		// Read the signal and statically store it once normalized into expCa and expV.
                gsl_matrix * expsRaw = gsl_matrix_alloc (steps, 22);
		{
			FILE * f = fopen ("dataPSM/experimentote.dat", "r");
			if(!f) {
				fprintf(stderr, "Calcium file not present. Aborting\n");
				exit(EXIT_FAILURE);
			}
			gsl_matrix_fscanf(f, expsRaw);
			fclose (f);
		}

		//Normalize
                normalize_gsl_matrix (expsRaw, experiments);
		gsl_matrix_free (expsRaw);
		
	
		//Write back the now normalized data
		{
			FILE * fn = fopen("dataPSM/N_experimentote.dat","w");
			if(!fn) {
				fprintf(stderr, "Couldn't write to normalized calcium file. Aborting\n");
				exit(EXIT_FAILURE);
		       	}
                        gsl_matrix_fprintf (fn, experiments, "%lf");
			fclose(fn);
		}	
		isFirstTime = false;
	}
	
	// We now have the experiments stored in expCa and expV.
	// Lets use the set of parameters and evaluate the solution against 
	// the experiments using the Slope Index.

	// If one of the thresholds is negative reject them like the plague!
	for (int i = 0; i < 26; i++) {
		if (x [i] < 0.0)
			return 500000.0;
	}

	// Tertiary nodes: 5, 10, 14, 15
	// x has 26 values. The first 22 are binary thresholds
	// the next 4 are the tertiary ones	
	double thr[44];
	// Populate the fake nodes
	for (int i = 0; i < 44; i+=2)  {
		thr[i] = 3.0;
	}
	// Populate the binary nodes
	for (int i = 1; i < 45; i+=2)  {
		thr[i] = fmod (x[i], 1);
	}
	// And now the ternary nodes
	thr[10] = fmod (x[22], 2);
	thr[20] = fmod (x[23], 2);
	thr[28] = fmod (x[24], 2);
	thr[30] = fmod (x[25], 2);
	
	// Check the ternary thresholds are actually >= binary ones
	if ( thr[10]<thr[11] || thr[20]<thr[21] || thr[28]<thr[29] || thr[30]<thr[31] )
		return 70000;

	// Calculate the Glass Dynamics
	gsl_matrix * dynamicsRaw = gsl_matrix_alloc (steps, 22);
	fullGlassDynamics2_withGSL_MATRIX (condInicial, 0.01, thr, alphas, chucherias, dynamicsRaw); 
	
	// Normalize the simulation
	gsl_matrix * dynamics = gsl_matrix_alloc (steps, 22);
	normalize_gsl_matrix (dynamicsRaw, dynamics);
	gsl_matrix_free (dynamicsRaw);

	// Normalized Experiment
	gsl_vector_const_view eV = gsl_matrix_const_column (experiments, 5);
        gsl_vector_const_view eLVA = gsl_matrix_const_column (experiments, 10);
        gsl_vector_const_view eHVA = gsl_matrix_const_column (experiments, 14);
        gsl_vector_const_view eC = gsl_matrix_const_column (experiments, 15);

        // Normalized Glass Simulation
	gsl_vector_const_view gV = gsl_matrix_const_column (dynamics, 5);
	gsl_vector_const_view gLVA = gsl_matrix_const_column (dynamics, 10);
	gsl_vector_const_view gHVA = gsl_matrix_const_column (dynamics, 14);
	gsl_vector_const_view gC = gsl_matrix_const_column (dynamics, 15);

	// Compare four channels of the experiment and the Glass simulation
	double thisSimulationScore = (
		0.25*MSE ( (&eV.vector)->data, (&gV.vector)->data, steps) +
		0.25*MSE ( (&eLVA.vector)->data, (&gLVA.vector)->data, steps) +
		0.25*MSE ( (&eHVA.vector)->data, (&gHVA.vector)->data, steps) +
		0.25*MSE ( (&eC.vector)->data, (&gC.vector)->data, (&gC.vector)->size)
		);

	gsl_matrix_free (dynamics);
	
	return thisSimulationScore;

	//        return 0.25*MSE( (&eGMP.vector)->data, cGMP->data, cGMP->size) +
	//      0.25*MSE( (&eV.vector)->data, voltage->data, voltage->size) +
        //       0.25*MSE( (&eAMP.vector)->data, cAMP->data, cAMP->size) +
        //       0.25*MSE( (&eC.vector)->data, calcium->data, calcium->size);

	// The score is computed using the Slope Index by Cho et al. 2006
	// Score = 2 - [SI(experimentalCalcium, simulationCalcium) + SI(experimentalVoltage, simulationVoltage)]
	// As we are optimizing with respect to two nodes, worst case is 2-(-2)=4. Best case is 2-2=0
//	return 2 - ( SlopeIndex (expCa->data, calcium->data, expCa->size) - SlopeIndex (expV->data, voltage->data, expV->size) );
       	//return 2 - ( SlopeIndex (expCa, calcium->data, calcium->size) + SlopeIndex (expV, voltage->data, voltage->size) );
}

