#include <vector>
#include <iostream>

#include <string.h> 
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
//#include <omp.h>

#include "produit_mat.h"
#include "transpose.h"
#include "maximum.h"
#include "gaussrand.h"
#include "scale.h"
#include "mhtquant.h"

void delete_kth_column(float *source, float *target, int source_line, int source_column, int k) {
	int i= k + 1;
	int j= k;
	int source_size = source_line * source_column;
	// first block if any
	if(k > 0) {
		memcpy(target, source, k * sizeof(float));
	}
	for(; i < (source_column*(source_line-1)) ; i += source_column, j += (source_column - 1)) {
		memcpy(&target[j], &source[i], (source_column - 1) * sizeof(float));
	}
	
	// last block if any
	if(i < source_size) {
		memcpy(&target[j], &source[i], (source_column - k-1) * sizeof(float));
	}
}


int mhtquant(int *Ktest,double *a,double *b, double *c,int *lig_a,int *p,double *F,int *IT,int *maxq,double *sigma)//,double *a, double *b,double *c,double *d, int *lig_a, int *col_a, int *col_b,double *epsilon,int *SUP)//,int *ijk)
{
	int iter;
	int IT2= *IT;
	double sigma2= *sigma; 
	
	
	/*on ne doit pas toucher a  XI2dep, U2dep, XII,*/
//#pragma omp parallel for
//	for(iter=0;iter<IT2;iter++)
//	{	
	
		int ijk,ktest,nl,nc_init,i,j,ij;
		int constante= *maxq+1;
		 
		/*a=XII
		 b=XI2dep, de dimension n*(p-ktest)	existe si ktest>0
		 c=U2dep, de dimension n*ktest		existe si ktest>0
		 */

	nl= *lig_a; /* =n*/
	ktest= *Ktest;
	nc_init= *p-ktest;	/*dimension de XI2dep */
	ijk= ktest;
	
	float *XII=(float*)malloc(nl* *p * sizeof(float));// float[nl* *p];
	float *XI2dep=(float*)malloc(nl* nc_init * sizeof(float));//new float[nl*nc];
	float *U2dep=(float*)malloc((nl*ktest*(ktest>0)+nl*(ktest==0)) * sizeof(float));//new float[nl*ktest*(ktest>0)+nl*(ktest==0)];
	
	int sup;
	sup=(int)ktest+pow(2.0, *maxq);
	
	float *U2final=(float*)malloc(nl* sup * sizeof(float));//new float[nl*sup];
	float *XI2=(float*)malloc(nl* nc_init * sizeof(float));//new float[nl*nc];
	float *U3=(float*)malloc(nl*sup * sizeof(float));//new float[nl*ncb];
//	float *U3=(float*)malloc((nl*ktest*(ktest>0)+nl*(ktest==0)) * sizeof(float));//new float[nl*ncb];
	
	float *epsilon=(float*)malloc(nl * sizeof(float));
	float *eps2=(float*)malloc(nl* sizeof(float));
	float *cor=(float*)malloc(nc_init * sizeof(float));//new float[nc];

	float *x=(float*)malloc(nl* nc_init * sizeof(float));//new float[nl*nc];
	float *cd=(float*)malloc((nc_init*ktest*(ktest>0)+nc_init*(ktest==0)) * sizeof(float));//new float[nc*ncb];
	float *cc=(float*)malloc((ktest*nc_init*(ktest>0)+nc_init*(ktest==0)) * sizeof(float));//new float[ncb*nc];
	float *dd=(float*)malloc(nl*nc_init * sizeof(float));//new float[nl*nc];
	float *colonne_dd=(float*)malloc(nl * sizeof(float));//new float[nl];
	float *vect=(float*)malloc(nl* sizeof(float));
	
	float *reg=(float*)malloc(sup * sizeof(float));
	float *FF=(float*)malloc(constante * sizeof(float));
	float *coef=(float*)malloc(nl* sup * sizeof(float));//new float[nl*maxit];  

	for(iter=0;iter<IT2;iter++)
	{	
		ijk=ktest;
		int nc=nc_init;
		if(ktest>0) /*si ktest>0, on se sert de b et c (XI2dep et U2dep) fourni en entrée*/
		{	
			for(i=0;i<nl*ktest;i++)
				U2dep[i]=c[i];
			
			for(i=0;i<nl*nc;i++)
			{	XI2dep[i]=b[i];
			}
			
		}
		else{
			for(i=0;i<nl* *p;i++)/*si ktest=0 on a besoin de a =XII*/
				XII[i]=a[i];}

		
		for(i=0;i<nl;i++)
			if(sigma2==0)
			{epsilon[i]=gaussrand();
			}else{
				epsilon[i]=sqrt(sigma2)*gaussrand();
			}

		float eps_moy=0;	
		for(i=0;i<nl;i++)
			eps_moy +=epsilon[i];
		
		for(i=0;i<nl;i++)
		{		
			eps2[i]=epsilon[i]-eps_moy/nl; /* eps2 est le epsilon moyenné*/
		}
	

		if(ijk==0)
		{	/*la matrice XII est deja centré reduite*/
			/* on calcule les correlations de epsilon avec XII, donc avec d scalé*/
			for(j=0;j<nc;j++)
			{	cor[j]=0;
				for(i=0;i<nl;i++)
					cor[j] +=eps2[i]*XII[i*nc+j];
				cor[j]=fabs(cor[j]);
			}
			int indice=0;
			float maxi=0;
			maximum(cor,nc,&indice,&maxi);		   
			for(i=0;i<nl;i++)
				U2dep[i]=XII[i*nc+indice]; /*contient la colonne qui nous interesse, en autre celle a enlever de XII pur avoir XI2dep */
			/*on enlève la colonne "indice" du tableau XII et ca devient le nouveau XI2dep, qui est de dimension une colonne en moins donc nc--*/
			delete_kth_column(XII, XI2dep, nl, nc, indice); 

			nc--;
			ijk=1;
		}
		
		
		/*a ce stade on a XI2dep et U2dep et c'est darty!!*/
		
		int ncb=ijk;		/*ncb est la dimension de U2final*/
		int ncc=ncb;
		
		/*nombre de colonne de U2dep, =ijk*/
		for(i=0;i<nl;i++)
		{for(j=0;j<ncb;j++)
			{U2final[i*(sup)+j]=U2dep[i*ncb+j];		/*U2final de dimension nl*(ktest+sup), mais avec juste de dimension nl*ncb pour le moment*/
				U3[i*ncb+j]=U2dep[i*ncb+j];			/*U3 est de dimension nl*ncb */
			}
		}
	
		
		memcpy(XI2, XI2dep, nl*nc * sizeof(float));
		
		/*on garde XI2dep et U2dep et on fait les calculs avec U2 et XI2*/
		for(ij=ijk;ij<sup;ij++)
		{			
			
			/*on transpose a pour faire le produit transp(a) * b*/
			transpose(XI2,x,nc,nl);

			/*x est la transpose de XI2, on peut donc faire le produit x * U2=transp(XI2) * U2  */
			produit_mat(x,U3,cd,nc,nl,ncb);	/*le resultat est dans cd : en ligne sa decomp dans la base b*/
			
			/*on transpose cd pour faire le produit b * transp(c)=u * transp(c) */
			transpose(cd,cc,ncb,nc);			/*le resultat est dans cc*/
			
			produit_mat(U3,cc,dd,nl,ncb,nc); /* resultat dans dd (d=X) : la decomp de chaque variable dans Vk (=b)*/
			

			
			for(i=0;i<nl*nc;i++)
				dd[i]=XI2[i]-dd[i];
		
/*			double sortie[ncc];
			for(j=0;j<ncc;j++)
			{
				sortie[j]=0;
				for(i=0;i<nl;i++)
				{
					sortie[j]=sortie[j]+ (double) dd[i*nc+4]*(double) U2final[i*sup+j];
				}
				printf("J= %ld\n",j);
				printf(" %lg \n",sortie[j]);
			}
			
*/
			
			/*dd=X_, il reste ce qu'il y a dans XI2 qui est dans l'orthogonal de Vk (=b)*/
			/*on cherche le max de fabs(d)=fabs(X_) par colonne, afin de mettre a 0 la colonne si tout est petit, 
			 ce qui correspondrait a des 0 en réalité*/
			for(j=0;j<nc;j++)
			{	int indice=0;
				float maxi=0;
				for(i=0;i<nl;i++)
					colonne_dd[i]=fabs(dd[i*nc+j]);
				maximum(colonne_dd,nl, &indice, &maxi); /*maximum de d (qui est de taille nl*nc), resultat dans indice et maxi*/
				if(maxi<(10^-10))
				{for(i=0;i<nl;i++)
				{colonne_dd[i]=0;
				dd[i*nc+j]=0;}
				}
			}
			

			/*on enleve la moyenne de chaque colonne de dd avant de scaler pour avoir sum carré =1*/
			
			scale(dd,nl,nc);
		

			/* on cherche a calculer les correlations de epsilon avec X_, donc avec d scalé*/
			produit_mat(eps2, dd, cor, 1, nl, nc);
			
			for(j=0;j<nc;j++)
			{	
				cor[j]=fabs(cor[j]);
			}
		
			
			/*on cherche le max des valeurs abs des correlations, ie le max de cor*/
			int indice=0;
			float maxi=0;
			maximum(cor,nc,&indice,&maxi);

			for(i=0;i<nl;i++)
			{	U3[i]=dd[i*nc+indice]; /*contient la colonne qui nous interesse, celle a enlever de XI2(=dd) et a rajouter dans U(=c) */
				U2final[i*(sup)+(ncc)]=U3[i];	/*on rajoute la colonne indice dans le tableau U2final, qui a maintenant une colonne en plus donc ncc++*/
			}
			ncb=1;


			/* U3 contient la colonne qui nous interesse, celle a enlever de XI2(=dd) et a rajouter dans U
			 on cherche l'orthogonal de XI2 par rapport a U3, puisqu'on est deja dans l'orthogonal de Vk*/	
			/* U2final est maintenant de dimension nl*ncc, nl ligne et ncc colonne, correspond au nouveau U*/

			ncc++; //=ij

			/*on enlève la colonne "indice" du tableau dd et ca devient le nouveau XI2, qu'on met dans XI2 de dimension une colonne en moins donc nc--*/
			delete_kth_column(dd, XI2, nl, nc, indice); 

			nc--;
			/*XI2 est maintenant de dimension nl*nc, nl ligne et nc colonne, correspond au nouveau XI2*/			
		}//fin ij<sup

		/*U2final est de dimension nl*ncc (=n* (ktest+SUP))
		 on veut faire la regression de epsilon sur U2final, donc de eps2(le epsilon centré) sur c
		ncc=sup
		 */
		
/*		double mean[ncc];
		double sortie[ncc];
		double sortie1[ncc];
		for(j=0;j<ncc;j++)
		{
			sortie[j]=0;
			sortie1[j]=0;
			mean[j]=0;
			for(i=0;i<nl;i++)
			{
				mean[j]=mean[j]+ (double) U2final[i*sup+j];
				sortie[j]=sortie[j]+ (double) U2final[i*sup+j]*(double) U2final[i*sup+j];
				sortie1[j]=sortie1[j]+ (double) U2final[i*sup+1]*(double) U2final[i*sup+j];
			}
			printf("J= %ld\n",j);
			printf("mean   %lg \n",mean[j]);
			printf("sortie   %lg \n",sortie[j]);
			printf("sortie 1   %lg \n",sortie1[j]);
		}
 
*/		
	//	printf("REG\n");
		for(j=0;j<ncc;j++)
		{	reg[j]=0;
			for(i=0;i<nl;i++)
				reg[j] +=eps2[i]*U2final[i*ncc+j];
	//	printf("%lg\n",reg[j]);
		}
		int m;
		for(m=0;m<constante;m++)
		{
			/*on calcule sum(beta^2) et le produit mat et la sum du produit mat*/
			float sum_bet=0;
			for(i=0;i<pow(2.0,m);i++)
				sum_bet+=reg[ktest+i]*reg[ktest+i];
			
	//		printf("sum_beta %lg\n",sum_bet);
			if(sigma2==0)
			{
				int maxit;
				maxit=(int)ktest+pow(2.0,m);

				for(j=0;j<maxit;j++)
				{
					for(i=0;i<nl;i++)
						coef[i*maxit+j]=U2final[i*ncc+j];
				}
				
				produit_mat(coef,reg ,vect, nl, maxit, 1);
				
				/*dd est la somme de la matrice de epsilon - U%*%beta au carré*/
				float som=0;
				for(i=0;i<nl;i++)
					som+=(eps2[i]-vect[i])*(eps2[i]-vect[i]);
				FF[m]=(nl-ktest-pow(2.0,m))*sum_bet/(pow(2.0,m)*som);			
			}else{FF[m]=sum_bet/sigma2;			
			}
		}
		/*il faut mettre FF dans F=rbind(F,FF)*/
	
		for(m=0;m<constante;m++)
		{F[iter*constante+m]=FF[m];
		}		
	}
	free(XII);//delete[] XII;
	free(XI2dep);//delete[] XI2dep;
	free(U2dep);//delete[] U2dep;
	free(U2final);//delete[] U2final;
	free(XI2);//delete[] XI2;
	free(U3);//delete[] U3;
	free(epsilon);
	free(eps2);
	free(cor);
	
	free(x);//delete[] x;
	free(cd);//delete[] cd;
	free(cc);//delete[] cc;
	free(dd);//delete[] dd;

	free(colonne_dd);//delete[] colonne_dd;
	free(vect);
	free(reg);
	free(FF);
	free(coef);
	return 0;
	
}





