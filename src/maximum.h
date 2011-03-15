/*
 *  prodmat.h
 *  try
 *
 *  Created by Florian on 12/01/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

float maximum(float *a, int taille, int *aa, float *bb)
{int i,indic;
	float maxim;
	
	indic=0;
	maxim=a[0];
	for(i=1;i<taille;i++)
	{if(a[i]>a[indic])
		{indic=i;
			maxim=a[i];
		}
	}
	*aa=indic;
	*bb=maxim;
	//printf("%ld\n",*aa);
	//printf("%lg\n",*bb);
	
	return 0;
}


