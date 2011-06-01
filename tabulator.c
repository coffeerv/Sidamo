/*
 *  tabulator.c
 *  glass Dynamics CL
 *
 *  Created by Rafael Verduzco Vázquez on 3/14/11.
 *  Copyright 2011 Universidad Nacional Autónoma de México. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "tabulator.h"

void readTable(const char * input)
{
	FILE * ptr = fopen(input, "r");
	if (ptr == NULL) exit(1);
	
	// initCond, numregs, binary/ternary, regsIndex, table
	char ic[10], nr[10], ns[10];
	fscanf(ptr, "%s\n", ic);
	fscanf(ptr, "%s\n", nr);
	fscanf(ptr, "%s\n", ns);
	
	int initialCondition = atoi(ic);
	int numRegs = atoi(nr);
	int numStates = atoi(ns);

	int regsIndex[numRegs];
	for (int i = 0; i < numRegs; i++) {
		char reg[10];
		fscanf(ptr, "%s\t", reg);
		regsIndex[i] = atoi(reg);
	}	
	
	
	discreteSettings d;
	d.initialConditions = atoi(ic);
	d.numRegulators = atoi(nr);
	d.numStates = atoi(ns);
	for (int i = 0; i < d.numRegulators; i++) {
		d.regulatorsIndex[i] = regsIndex[i];
	}
	
		
	printf("IC = %d, NR = %d, NS = %d\nRegulators = ", initialCondition, numRegs, numStates);
	for (int i = 0; i < numRegs; i++) {
		printf("%d\t", regsIndex[i]);
	}
	printf("\n");
	
	
	
	
	printf("IC = %d, NR = %d, NS = %d\nRegulators = ", d.initialConditions, d.numRegulators, d.numStates);
	for (int i = 0; i < numRegs; i++) {
		printf("%d\t", d.regulatorsIndex[i]);
	}
	printf("\n");
	
	
	fclose(ptr);
}