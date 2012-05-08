/*
 *  transpose.cpp
 *  mhtp
 *
 *  Created by Florian on 02/09/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

/*
#include "transpose.h"

float transpose(float *a,float *b,int lig_a,int col_a)
{int i,j;
	for(i=0;i< lig_a;i++)
		for(j=0;j<col_a;j++) //lig_b=lig_a
			b[i*col_a +j]= a[j*lig_a +i];
	return 0;
}
*/

#include "transpose.h"

float transpose(float *a,float *b,int lig_a,int col_a)
{int i,j;
	for(i=0;i< lig_a;i++)
		for(j=0;j<col_a;j++) //lig_b=lig_a
			b[i*col_a +j]= a[j*lig_a +i];
	return 0;
}

