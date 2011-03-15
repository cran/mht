/*
 *  prod.h
 *  try
 *
 *  Created by Florian on 12/01/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
float transpose(float *a,float *b,int lig_a,int col_a)
{int i,j;
for(i=0;i< lig_a;i++)
	for(j=0;j<col_a;j++) //lig_b=lig_a
		b[i*col_a +j]= a[j*lig_a +i];
	return 0;
}
