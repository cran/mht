/*
 *  scale.cpp
 *  mhtpp
 *
 *  Created by Florian on 06/09/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "scale.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>

void scale(float *source, int source_line, int source_column) {
	
	int i,j;
//	float moyenne[source_column];
//	float var[source_column];

	
//	float moyenne,var;
	for(j=0;j<source_column;j++)
	{float moyenne;
		float var;
		var=(float) 0;
		moyenne=(float) 0;//[j]=0;
		//var[j]=0;
		for(i=0;i<source_line;i++)
			moyenne= moyenne + (float) source[i*source_column+j] / (float) source_line;
		//printf("moyenne j %lg\n",moyenne);
		for(i=0;i<source_line;i++)
		{	source[i*source_column+j] -=moyenne;
			var= var + (float) source[i*source_column+j]* (float)source[i*source_column+j];
		}
	//	printf("var j %lg\n",var);
		
		for(i=0;i<source_line;i++)
		{	source[i*source_column+j]= (float) source[i*source_column+j] / (float) sqrt(fabs(var));
		}
	}
	
//	free(moyenne);
//	free(var);
	
}
