/*
 *  gaussrand.h
 *  try
 *
 *  Created by Florian on 17/01/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

double gaussrand()
{
	int NSUM=50;
	double x = 0;
	int i;
	for(i = 0; i < NSUM; i++)
		x += (double)rand() / RAND_MAX;
	
	x -= NSUM / 2.0; /*on retranche l'esperance de la somme de NSUM variable uniforme entre 0 et 1 (d'esperance 0+1/2)*/
	x /= sqrt(NSUM / 12.0); /*on divise par la racine carré de la variance(la variance d'une uniforme est (a+b)^2/12)*/
	
	return x;
}
