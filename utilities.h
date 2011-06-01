/*
 *  utilities.h
 *  pruebaGSL
 *
 *  Created by Rafael Verduzco Vázquez on 1/26/11.
 *  Copyright 2011 Universidad Nacional Autónoma de México. All rights reserved.
 *
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

int compareints (const void * a, const void * b);

double * vectorDouble (const int n);
void matrixDouble (const int rows, const int cols, double ** matrix);
void freeMatrixDouble (double ** matrix, int rows);
int * vectorInt(const int n);
int ** matrixInt(const int rows, const int cols);

double * readSignal(const char * dataFile, const int dataSize); 

void normalize_gsl_vector (const gsl_vector * input, gsl_vector * output);
void normalize_gsl_matrix (gsl_matrix * input, gsl_matrix * output);
