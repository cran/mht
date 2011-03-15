float produit_mat(float *a, float *b,float *c, int lig_a, int col_a, int col_b)
{int i,j;

	
	for (i = 0; i < lig_a; i++)
		for (j = 0; j < col_b; j++)
		{
			int z;
			c[i* col_b+j] = 0;
			for (z = 0; z < col_a; z++)
				c[i* col_b+j] +=a[i* col_a+z] * b[z* col_b+j];
		}

	return *c;
}




