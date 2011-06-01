/*
 *  utilities.c
 *  pruebaGSL
 *
 *  Created by Rafael Verduzco Vázquez on 1/26/11.
 *  Copyright 2011 Universidad Nacional Autónoma de México. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
#include <gsl/gsl_rng.h>

#include "utilities.h"

extern gsl_rng * r;

int compareints (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}

double * vectorDouble (const int n) {
	//double * ptr = calloc(n, sizeof(*ptr) );
	double * ptr = (double *) malloc((size_t) (n * sizeof(*ptr)));
	if (ptr == NULL) {
		fprintf(stderr, "Failed to allocate memory for vector\n");
		exit(1);
	}
	return ptr;
}

void matrixDouble (const int rows, const int cols, double ** matrix) {
	// double ** matrix = calloc(rows, sizeof(double *) );
	// double ** 
	matrix = (double **) malloc((size_t) (rows * sizeof (double *)));
	if (matrix == NULL) {
		fprintf(stderr, "Failed to allocate memory for matrix (failure at first dim)\n");
		exit(2);
	}
	for (int i = 0; i < rows; i++) {
		//matrix[i] = calloc(cols, sizeof(double) );
		matrix[i] = (double *) malloc((size_t) (cols * sizeof (double) ));
		if (matrix[i] == NULL) {
			fprintf(stderr, "Failed to allocate memory for matrix (1st dim ok, failure at 2nd dim)\n");
			exit(3);
		}
	}
}

int * vectorInt (const int n) {
	int * ptr = (int *) malloc((size_t) (n * sizeof (*ptr)) );
	if (ptr == NULL) {
		fprintf(stderr, "Failed to allocate memory for vector (integer)\n");
		exit(6);
	}
	return ptr;
}

int ** matrixInt (const int rows, const int cols) {
	int ** matrix = (int **) malloc((size_t) (rows * sizeof (int *)) );
	if (matrix == NULL) {
		fprintf(stderr, "Failed to allocate memory for matrix (failure at first dim)\n");
		exit(2);
	}
	for (int i = 0; i < rows; i++) {
		matrix[i] = (int *) malloc((size_t) (cols * sizeof (int)) );
		if (matrix[i] == NULL) {
			fprintf(stderr, "Failed to allocate memory for matrix (1st dim ok, failure at 2nd dim)\n");
			exit(3);
		}
	}
	return matrix;
}

double * readSignal (const char * dataFile, const int dataSize) {
	
	double * _signal, fake;
	_signal = vectorDouble(dataSize);
	
	FILE * fp;
	fp = fopen(dataFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "The file containing the signal could not be opened. Aborting.\n");
		exit(4);
	}
	
	for (int i = 0; i < dataSize; i++) {
		fscanf(fp, "%lf\t%lf\n", &fake, &_signal[i]);
	}
	
	printf("Successfully read data from file:\n%s\n", dataFile);
	fclose(fp);
	
	return _signal;
}

void freeMatrixDouble (double ** matrix, int rows)
{
    for(int i = 0; i < rows; i++)
    {
        free (matrix[i]);
    }
    free (matrix);
}

void normalize_gsl_vector (const gsl_vector * input, gsl_vector * output) {
    double m, sd;
    m = gsl_stats_mean (input->data, input->stride, input->size);
    sd = gsl_stats_sd_m (input->data, input->stride, input->size, m);
    for (int i = 0; i < (int)input->size; i++) {
    	    gsl_vector_set (output, i, (gsl_vector_get(input, i) - m)/sd);
    }   
}

void normalize_gsl_matrix (gsl_matrix * input, gsl_matrix * output) {
    for (int col = 0; col < (int)input->size2; col++) {
        gsl_vector_const_view currCol = gsl_matrix_const_column (input, col);
	gsl_vector * theCol = gsl_vector_alloc (input->size1);
	normalize_gsl_vector (&currCol.vector, theCol);
	gsl_matrix_set_col (output, col, theCol);
        gsl_vector_free (theCol);
    }

}

