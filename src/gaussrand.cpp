/*
 *  gaussrand.cpp
 *  mhtp
 *
 *  Created by Florian on 02/09/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "gaussrand.h"
#include <vector>
#include <iostream>

using namespace std;

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
float gaussrand()
{
	int NSUM=50;
	float x = 0;
	int i;
	for(i = 0; i < NSUM; i++)
		x += (float)rand() / RAND_MAX;
	
	x -= NSUM / 2.0; /*on retranche l'esperance de la somme de NSUM variable uniforme entre 0 et 1 (d'esperance 0+1/2)*/
	x /= sqrt(NSUM / 12.0); /*on divise par la racine carré de la variance(la variance d'une uniforme est (a+b)^2/12)*/
	
	return x;
}
