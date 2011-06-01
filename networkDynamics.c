/*
 *  networkDynamics.c
 *  pruebaGSL
 *
 *  Created by Rafael Verduzco Vázquez on 2/4/11.
 *  Copyright 2011 Universidad Nacional Autónoma de México. All rights reserved.
 *
 */

#include "networkDynamics.h"

//#include <stdio.h>
//#include <stdlib.h>
//#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_matrix.h>

//#include "utilities.h"


int Heaviside(const int i, const double * contState, const double * thresholds)
{
	return (contState[i] - thresholds[i]) > 0 ? 1 : 0;
}

int Heaviside2(const int i, const double * contState, const double thrH, const double thrL)
{
	if (contState[i] >= thrH) {
		return 2;
	}else if (contState[i] >= thrL) {
		return 1;
	}else {
		return 0;
	}
}

//double ** glassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), int (*Heaviside) (const int, const double *, const double *)) {
//void glassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), int (*Heaviside) (const int, const double *, const double *), double ** cS ) {
void glassDynamics
(const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, 
 const double * thresholds, const double * alphas, int (* discreteMapping) (const int, const int *),
 int nodeToTrack, double * cand)
{	
	double cS[2][numNodes];
	int dS[2][numNodes]; 
	
	// Copy the initial conditions 
	for (int i = 0; i < numNodes; i++) {
		cS[0][i] = initialConditions[i];
	}
	// Set the initial condition as the first point of the candidate 'cand'
	if (nodeToTrack > numNodes) {
		fprintf(stderr, "Invalid node to track. Aborting\n");
		exit(EXIT_FAILURE);
	}
	cand[0] = initialConditions[nodeToTrack];
	
	// Use Eulers method to solve the equations
	int t = 1; 
	
	while (t < dataPoints) {
		// Calcular la discretización de las variables continuas
		for (int i = 0; i < numNodes; i++) {
			dS[0][i] = Heaviside(i, cS[0], thresholds);
		}
		// Una vez discretizado, calcular el mapeo discreto.
		// Este valor es el que se usa para calcular la dinámica de Glass al tiempo t
		for (int i = 0; i < numNodes; i++) {
			dS[1][i] = discreteMapping(i,dS[0]);
		}
		// Resolver el sistema de ecuaciones diferenciales con el método de Euler simple
		for (int i = 0; i < numNodes; i++) {
			cS[1][i] = cS[0][i] + (alphas[i] * (dS[1][i] - cS[0][i]) * stepSize);
		}
		// Update the system
		cand[t] = cS[1][nodeToTrack];
		for (int i = 0; i < numNodes; i++) {
			cS[0][i] = cS[1][i];
		}
		t++;
	}
}

void glassDynamics2
(const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, 
 const double * thresholds, const double * alphas, int (* discreteMapping) (const int, const int *),
 int nodeToTrack, double * cand)
{	
	double cS[2][numNodes];
	int dS[2][numNodes]; 
	
	// Copy the initial conditions 
	for (int i = 0; i < numNodes; i++) {
		cS[0][i] = initialConditions[i];
	}
	// Set the initial condition as the first point of the candidate 'cand'
	if (nodeToTrack > numNodes) {
		fprintf(stderr, "Invalid node to track. Aborting\n");
		exit(EXIT_FAILURE);
	}
	cand[0] = initialConditions[nodeToTrack];
	
	// Use Eulers method to solve the equations
	int t = 1; 
	
	while (t < dataPoints) {
		// Calcular la discretización de las variables continuas
		for (int i = 0; i < numNodes; i++) {
			dS[0][i] = Heaviside2(i, cS[0], thresholds[i*2], thresholds[(i*2)+1]);
		}
		// Una vez discretizado, calcular el mapeo discreto.
		// Este valor es el que se usa para calcular la dinámica de Glass al tiempo t
		for (int i = 0; i < numNodes; i++) {
			dS[1][i] = discreteMapping(i,dS[0]);
		}
		// Resolver el sistema de ecuaciones diferenciales con el método de Euler simple
		for (int i = 0; i < numNodes; i++) {
			cS[1][i] = cS[0][i] + (alphas[i] * (dS[1][i] - cS[0][i]) * stepSize);
		}
		// Update the system
		cand[t] = cS[1][nodeToTrack];
		for (int i = 0; i < numNodes; i++) {
			cS[0][i] = cS[1][i];
		}
		t++;
	}
}


void fullGlassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), double *dyn) {	
	double cS[2][numNodes];
	int dS[2][numNodes]; 
	
	// Copy the initial conditions 
	for (int i = 0; i < numNodes; i++) {
		dyn[i] = initialConditions[i];
		cS[0][i] = initialConditions[i];
	}
		
	// Use Eulers method to solve the equations
	int t = 1, t2 = numNodes; 
	
	while (t < dataPoints) {
		// Calcular la discretización de las variables continuas
		for (int i = 0; i < numNodes; i++) {
			//dS[t-1][i] = Heaviside(i, cS[t-1], thresholds);
			dS[0][i] = Heaviside(i, cS[0], thresholds);
		}
		// Una vez discretizado, calcular el mapeo discreto.
		// Este valor es el que se usa para calcular la dinámica de Glass al tiempo t
		for (int i = 0; i < numNodes; i++) {
			dS[1][i] = discreteMapping(i, dS[0]);
		}
		// Resolver el sistema de ecuaciones diferenciales con el método de Euler simple
		for (int i = 0; i < numNodes; i++) {
			cS[1][i] = cS[0][i] + (alphas[i] * (dS[1][i] - cS[0][i]) * stepSize);
		}
		// Update the system
		for (int i = 0; i < numNodes; i++) {
			dyn[t2+i] = cS[1][i];
			cS[0][i] = cS[1][i];
		}
		t++;
		t2 += numNodes;
	}
}

void fullGlassDynamics2 (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), double *dyn) {	
	double cS[2][numNodes];
	int dS[2][numNodes]; 
	
	// Copy the initial conditions 
	for (int i = 0; i < numNodes; i++) {
		dyn[i] = initialConditions[i];
		cS[0][i] = initialConditions[i];
	}
	
	// Use Eulers method to solve the equations
	int t = 1, t2 = numNodes; 
	
	while (t < dataPoints) {
		// Calcular la discretización de las variables continuas
		for (int i = 0; i < numNodes; i++) {
			//dS[t-1][i] = Heaviside(i, cS[t-1], thresholds);
			dS[0][i] = Heaviside2(i, cS[0], thresholds[i*2], thresholds[(i*2)+1]);
		}
		// Una vez discretizado, calcular el mapeo discreto.
		// Este valor es el que se usa para calcular la dinámica de Glass al tiempo t
		for (int i = 0; i < numNodes; i++) {
			dS[1][i] = discreteMapping(i, dS[0]);
		}
		// Resolver el sistema de ecuaciones diferenciales con el método de Euler simple
		for (int i = 0; i < numNodes; i++) {
			cS[1][i] = cS[0][i] + (alphas[i] * (dS[1][i] - cS[0][i]) * stepSize);
		}
		// Update the system
		for (int i = 0; i < numNodes; i++) {
			dyn[t2+i] = cS[1][i];
			cS[0][i] = cS[1][i];
		}
		t++;
		t2 += numNodes;
	}
}

void fullGlassDynamics2_withGSL_MATRIX (const double * initialConditions, const double stepSize, const double * thresholds, const double * alphas, int (* discreteMapping) (const int, const int *), gsl_matrix * dyn) {
	int numNodes = (int) dyn->size2;
	int dataPoints = (int) dyn->size1;	
	double cS[2][numNodes];
	int dS[2][numNodes]; 
	
	// Copy the initial conditions 
	for (int node = 0; node < numNodes; node++) {
		//dyn[i] = initialConditions[i];
		gsl_matrix_set (dyn, 0, node, initialConditions[node]);
		cS[0][node] = initialConditions[node];
	}
	
	// Use Eulers method to solve the equations
	int t = 1;//, t2 = numNodes; 
	
	while (t < dataPoints) {
		// Calcular la discretización de las variables continuas
		for (int i = 0; i < numNodes; i++) {
			//dS[t-1][i] = Heaviside(i, cS[t-1], thresholds);
			dS[0][i] = Heaviside2(i, cS[0], thresholds[i*2], thresholds[(i*2)+1]);
		}
		// Una vez discretizado, calcular el mapeo discreto.
		// Este valor es el que se usa para calcular la dinámica de Glass al tiempo t
		for (int i = 0; i < numNodes; i++) {
			dS[1][i] = discreteMapping(i, dS[0]);
		}
		// Resolver el sistema de ecuaciones diferenciales con el método de Euler simple
		for (int i = 0; i < numNodes; i++) {
			cS[1][i] = cS[0][i] + (alphas[i] * (dS[1][i] - cS[0][i]) * stepSize);
		}
		// Update the system
		for (int node = 0; node < numNodes; node++) {
			gsl_matrix_set (dyn, t, node, cS[1][node]);
			cS[0][node] = cS[1][node];
		}
		//for (int i = 0; i < numNodes; i++) {
		//	dyn[t2+i] = cS[1][i];
		//	cS[0][i] = cS[1][i];
		//}
		t++;
		//t2 += numNodes;
	}
}




int ** discDynamics
(const int * initialConditions, const int numNodes, const int dataPoints, int (* discreteMapping) (int, int *) )
{
	int ** din;
	
	din = matrixInt(dataPoints, numNodes);
	
	for (int i = 0; i < numNodes; i++) {
		din[0][i] = initialConditions[i];
	}
	
	int t = 1;
	
	while (t < dataPoints) {
		for (int i = 0; i < numNodes; i++) {
			din[t][i] = discreteMapping(i, din[t-1]);
		}
		t++;
	}
	return din;
}


double sign (const double x) {
	if(x == 0.0)
		return 0.0;
	else 
		return x < 0.0 ? -1.0 : 1.0;
}

double PearsonCorrelation(const double * data1, const double * data2, const int n) {
	return gsl_stats_correlation(data1, 1, data2, 1, n);
}

double SlopeIndex (const double * data1, const double * data2, const int n) {
	double si = 0.0;
	for(int i = 0; i < n-2; i++)
		si += sign((data2[i+1] - data2[i]) / (data1[i+1] - data1[i]));
		return si / (n - 2);
}

double fitnessPearsonCorrelation (const double * data1, const double * data2, const int n) {
	return 1 - gsl_stats_correlation(data1, 1, data2, 1, n);
}

double fitnessSlopeIndex (const double * data1, const double * data2, const int n) {
	double si = 0.0;
	for(int i = 0; i < n-2; i++)
		si += sign((data2[i+1] - data2[i]) / (data1[i+1] - data1[i]));
	return 1 - (si / (n - 2));
}

//double fitnessPowerSpectra (const double * x1, const double * x2, const int n) //{
//	printf("Yet to be implemented\nTerminating program...\n");
//	exit(0);
//}

double correlation(const double *x1, const double *x2, const int n){
	double corr = 0;
	for (int i = 0; i < n; i++) {
		corr += x1[i] * x2[i];
	}
	corr /= n;
	return corr;
}

double MSE (const double * realData, const double * estimator, const int n) {
	double _mse = 0.0;
	for (int i = 0; i < n; i++) {
		_mse += (estimator[i] - realData[i])*(estimator[i] - realData[i]);
	}
	return _mse/n;
}


