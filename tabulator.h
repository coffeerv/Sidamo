/*
 *  tabulator.h
 *  glass Dynamics CL
 *
 *  Created by Rafael Verduzco Vázquez on 3/14/11.
 *  Copyright 2011 Universidad Nacional Autónoma de México. All rights reserved.
 *
 */

typedef struct ds{
	int initialConditions;
	int numRegulators;
	int numStates;
	int * regulatorsIndex;
}discreteSettings;

void readTable(const char * pt);//, int ** regulatoryTable);
