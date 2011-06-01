/*
 *  networkDynamics.h
 *  pruebaGSL
 *
 *  Created by Rafael Verduzco Vázquez on 2/4/11.
 *  Copyright 2011 Universidad Nacional Autónoma de México. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>

#include "utilities.h"

int Heaviside(const int i, const double * contState, const double * thresholds);
int Heaviside2(const int i, const double * contState, const double thrH, const double thrL);

/* Dynamics solvers */
//double ** glassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const double * alphas, int (* discreteMapping) (const int, const int *), int (* Heaviside) (const int, const double *, const double *));
//void glassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const double * alphas, int (* discreteMapping) (const int, const int *), int (* Heaviside) (const int, const double *, const double *), double ** cS);
void glassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), int nodeToTrack, double * cand );
void glassDynamics2 (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), int nodeToTrack, double * cand );
void fullGlassDynamics (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), double *dyn);	
void fullGlassDynamics2 (const double * initialConditions, const int numNodes, const int dataPoints, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), double *dyn);
void fullGlassDynamics2_withGSL_MATRIX (const double * initialConditions, const double stepSize, const double * thresholds, const	double * alphas, int (* discreteMapping) (const int, const int *), gsl_matrix * dyn);	
int ** discDynamics (const int * initialConditions, const int numNodes, const int dataPoints, int (* discreteMapping) (int, int *) );

/* Fitness functions */


double sign (const double x);
double PearsonCorrelation(const double * data1, const double * data2, const int n);
double SlopeIndex (const double * data1, const double * data2, const int n);
double fitnessPearsonCorrelation (const double * data1, const double * data2, const int n);
double fitnessSlopeIndex (const double * data1, const double * data2, const int n);
//double fitnessPowerSpectra (const double * x1, const double * x2, const int n);
double correlation(const double *x1, const double *x2, const int n);
double MSE (const double *realData, const double *estimator, const int n);


//double sign (const double n);
//double fitnessPearsonCorrelation (const double * data1, const double * data2, const int n);
//double fitnessSlopeIndex (const double * data1, const double * data2, const int n);
//double fitnessPowerSpectra (const double * data1, const double * data2, const int n);
