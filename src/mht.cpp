#include <vector>
#include <iostream>

using namespace std;

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>


#include "produit_mat.h"
#include "transpose.h"
#include "maximum.h"
#include "gaussrand.h"

int main(int *Ktest,double *a,double *b, double *c,int *lig_a,int *p,double *F,int *IT,int *maxq)//,double *a, double *b,double *c,double *d, int *lig_a, int *col_a, int *col_b,double *epsilon,int *SUP)//,int *ijk)
{
	int iter;
	int IT2= *IT;
	
	int nbr;
	nbr=2;
	/*on ne doit pas toucher a  XI2dep, U2dep, XII,*/
	for(iter=0;iter<IT2;iter++)
	{	
		int ijk,ktest,nl,nc,i,j,ij;
		int constante;
		
		/*a=XII
		 b=XI2dep, de dimension n*(p-ktest)	existe si ktest>0
		 c=U2dep, de dimension n*ktest		existe si ktest>0
		 */
		nl= *lig_a; /* =n*/
		ktest= *Ktest;
		nc= *p-ktest;	/*dimension de XI2dep */

		constante= *maxq+1;
		float *XII=new float[nl* *p];
		float *XI2dep=new float[nl*nc];
		float *U2dep=new float[nl*ktest*(ktest>0)+nl*(ktest==0)];
		
		
		if(ktest>0) /*si ktest>0, on se sert de b et c (XI2dep et U2dep) fourni en entrée*/
		{	
			for(i=0;i<nl*ktest;i++)
				U2dep[i]=c[i];
			for(i=0;i<nl*nc;i++)
				XI2dep[i]=b[i];
		}
		else{
			for(i=0;i<nl* *p;i++)/*si ktest=0 on a besoin de a =XII*/
				XII[i]=a[i];}
		
		nl= *lig_a; /* =n*/
		ktest= *Ktest;
		nc= *p-ktest;	/*dimension de XI2dep */
		ijk= ktest;
		std::vector<float> epsilon(nl);
		for(i=0;i<nl;i++)
			epsilon[i]=gaussrand();
		float norm_eps=0;
		for(i=0;i<nl;i++)
			norm_eps +=epsilon[i]*epsilon[i];
		float eps_moy=0;	
		for(i=0;i<nl;i++)
			eps_moy +=epsilon[i];

		std::vector<float> eps2(nl);
		for(i=0;i<nl;i++)
			eps2[i]=epsilon[i]-eps_moy/nl; /* eps2 est le epsilon moyenné*/

		float *cor=new float[nc];

		if(ijk==0)
		{	/*la matrice XII est deja centré reduite*/
			/* on calcule les correlations de epsilon avec XII, donc avec d scalé*/
			for(j=0;j<nc;j++)
			{	cor[j]=0;
				for(i=0;i<nl;i++)
					cor[j] +=eps2[i]*XII[i*nc+j];
				cor[j] /=sqrt(norm_eps);
				cor[j]=fabs(cor[j]);
			}
			int indice=0;
			float maxi=0;
			maximum(cor,nc,&indice,&maxi);			
			for(i=0;i<nl;i++)
				U2dep[i]=XII[i*nc+indice]; /*contient la colonne qui nous interesse, en autre celle a enlever de XII pur avoir XI2dep */
			/*on enlève la colonne "indice" du tableau XII et ca devient le nouveau XI2dep, qui est de dimension une colonne en moins donc nc--*/
			float h;
			j=0;
			for(i=0 ; i< (nc*nl) ; i++)
			{	h=(i-indice);
				h/=(nc);
				if (fmod(h,1)!=0) 
				{XI2dep[j] = XII[i];
					j++; 
				}
			}
			nc--;
			ijk=1;
		}
		delete[] cor;
		
		
		/*a ce stade on a XI2dep et U2dep et c'est darty!!*/
		
		int sup;
		
		sup=ktest+pow(2, *maxq);
		int ncb=ijk;		/*ncb est la dimension de U2final*/
		int ncc=ncb;
		
		float *U2final=new float[nl*sup];
		float *XI2=new float[nl*nc];
		float *U3=new float[nl*ncb];
		
		/*nombre de colonne de U2dep, =ijk*/
		for(i=0;i<nl;i++)
		{for(j=0;j<ncb;j++)
		{U2final[i*(sup)+j]=U2dep[i*ncb+j];		/*U2final de dimension nl*(ktest+sup), mais avec juste de dimension nl*ncb pour le moment*/
			U3[i*ncb+j]=U2dep[i*ncb+j];			/*U3 est de dimension nl*ncb */
		}
			
			for(j=ncb;j<sup;j++)
				U2final[i*(sup)+j]=0;
			
			for(j=0;j<nc;j++)
				XI2[i*nc+j]=XI2dep[i*nc+j]; /*XI2 de dimension nl*nc*/
		}
		
		/*on garde XI2dep et U2dep et on fait les calculs avec U2 et XI2*/
		
		
		for(ij=ijk;ij<sup;ij++)
		{
			/*on transpose a pour faire le produit transp(a) * b*/
			float *x=new float[nl*nc];
			transpose(XI2,x,nc,nl);
			
			/*x est la transpose de XI2, on peut donc faire le produit x * U2=transp(XI2) * U2  */
			float *cd=new float[nc*ncb];
			produit_mat(x,U3,cd,nc,nl,ncb);	/*le resultat est dans cd : en ligne sa decomp dans la base b*/
			
			/*on transpose cd pour faire le produit b * transp(c)=u * transp(c) */
			float *cc=new float[ncb*nc];
			transpose(cd,cc,nc,ncb);			/*le resultat est dans cc*/
			
			float *dd=new float[nl*nc];
			produit_mat(U3,cc,dd,nl,ncb,nc); /* resultat dans dd (d=X) : la decomp de chaque variable dans Vk (=b)*/
			
			float *ddd=new float[nl*nc];
			for(i=0;i<nl*nc;i++)
				ddd[i]=XI2[i]-dd[i];
			for(i=0;i<nl*nc;i++)
				dd[i]=ddd[i];	
			/*dd=X_, il reste ce qu'il y a dans XI2 qui est dans l'orthogonal de Vk (=b)*/
			delete[] cc;
			delete[] ddd;
			delete[] cd;
			delete[] x;
			/*on cherche le max de fabs(d)=fabs(X_) par colonne, afin de mettre a 0 la colonne si tout est petit, 
			 ce qui correspondrait a des 0 en réalité*/
			for(j=0;j<nc;j++)
			{	int indice=0;
				float maxi=0;
				float *colonne_dd=new float[nl];
				for(i=0;i<nl;i++)
					colonne_dd[i]=fabs(dd[i*nc+j]);
				maximum(colonne_dd,nl, &indice, &maxi); /*maximum de d (qui est de taille nl*nc), resultat dans indice et maxi*/
				if(maxi<(10^-10))
				{for(i=0;i<nl;i++)
				{colonne_dd[i]=0;
				dd[i*nc+j]=0;}
				}
				delete[] colonne_dd;
			}
			/*on enleve la moyenne de chaque colonne de dd avant de scaler pour avoir sum carré =1*/
			std::vector<float> moyenne_dd(nc);
			for(j=0;j<nc;j++)
			{moyenne_dd[j]=0;
				for(i=0;i<nl;i++)
					moyenne_dd[j] +=dd[i*nc+j];
				moyenne_dd[j] /=nl;
				for(i=0;i<nl;i++)
					dd[i*nc+j] -=moyenne_dd[j];
			}
			
			/*on scale la matrice dd par colonne pour avoir somme de la colonne au carré =1*/
			std::vector<float> scale (nc);
			for(j=0;j<nc;j++)
			{scale[j]=0;
				for(i=0;i<nl;i++)
					scale[j] +=dd[i*nc+j]*dd[i*nc+j];
				for(i=0;i<nl;i++)
					dd[i*nc+j]=dd[i*nc+j]/sqrt(fabs(scale[j]));
			}
			
			/* on cherche a calculer les correlations de epsilon avec X_, donc avec d scalé*/
		
			float *cor=new float[nc];

			for(j=0;j<nc;j++)
			{	cor[j]=0;
				for(i=0;i<nl;i++)
					cor[j] +=eps2[i]*dd[i*nc+j];
				cor[j] /=sqrt(norm_eps);
				cor[j]=fabs(cor[j]);
			}
			
			/*on cherche le max des valeurs abs des correlations, ie le max de cor[]*/
			
			int indice=0;
			float maxi=0;
			maximum(cor,nc,&indice,&maxi);
			
			delete[]cor;
			std::vector<float> vect (nl);

			for(i=0;i<nl;i++)
				vect[i]=dd[i*nc+indice]; /*contient la colonne qui nous interesse, celle a enlever de XI2(=dd) et a rajouter dans U(=c) */
			
			/*le nouveau U3 est vect, on cherche l'orthogonal de XI2 par rapport a vect, puisqu'on est deja dans l'orthogonal de Vk*/
			ncb=1;
			for(i=0;i<nl*ncb;i++)
			{	U3[i]=vect[i];
			}
			
			/*on rajoute la colonne indice dans le tableau U2final, qui a maintenant une colonne en plus donc ncc++*/
			
			for(i=0;i<nl;i++)
				U2final[i*(sup)+ncc]=vect[i];
			
			
			ncc++;
			
			/* U2final est maintenant de dimension nl*ncc, nl ligne et ncc colonne, correspond au nouveau U*/
			
			
			
			
			/*on enlève la colonne "indice" du tableau XI2 et ca devient le nouveau XI2, qu'on met dans XI2 de dimension une colonne en moins donc nc--*/
			float h;
			j=0;
			for(i=0 ; i< (nc*nl) ; i++)
			{	h=(i-indice);
				h/=(nc);
				
				if (fmod(h,1)!=0) 
				{XI2[j] = dd[i];
					j++; 
				}
			}
			for(i=(nc-1)*nl;i<nc*nl;i++) 
				XI2[i] = 0;
			
			nc--;
			
			/*d est maintenant de dimension nl*nc, nl ligne et nc colonne, correspond au nouveau XI2*/
			delete[] dd;
		}//fin ij<sup
		
		/*U2final est de dimension nl*ncc (=n* (ktest+SUP))
		 
		 on veut faire la regression de epsilon sur U2final, donc de eps2(le epsilon centré) sur c
		 */
		
		std::vector<float> reg (ncc);
		
		for(j=0;j<ncc;j++)
		{	reg[j]=0;
			for(i=0;i<nl;i++)
				reg[j] +=eps2[i]*U2final[i*ncc+j];
		}
		
		std::vector<float> FF (constante);
		int m;
		for(m=0;m<constante;m++)
		{
			/*on calcule sum(beta^2) et le produit mat et la sum du produit mat*/
			float sum_bet=0;
			for(i=0;i<pow(2,m);i++)
				sum_bet+=reg[ktest+i]*reg[ktest+i];
			
			int maxit;
			maxit=ktest+pow(2,m);
			
			float *dd=new float[nl*maxit];  
			float *bet=new float[maxit];
			for(j=0;j<maxit;j++)
			{
				for(i=0;i<nl;i++)
					dd[i*maxit+j]=U2final[i*ncc+j];
				bet[j]=reg[j];
			}
			
			float *cc=new float[nl];
			produit_mat(dd,bet ,cc, nl, maxit, 1);
			
			
			
			/*dd est la somme de la matrice de epsilon - U%*%beta au carré*/
			float som=0;
			for(i=0;i<nl;i++)
				som+=(eps2[i]-cc[i])*(eps2[i]-cc[i]);
			FF[m]=(nl-ktest-pow(2,m))*sum_bet/(pow(2,m)*som);			
			delete[] cc;
			delete[] dd;
			delete[] bet;
		}
		/*il faut mettre FF dans F=rbind(F,FF)*/
		
		for(m=0;m<constante;m++)
		{F[iter*constante+m]=FF[m];
		}
		
		delete[] U2final;
		delete[] XI2;
		delete[] U3;
		
		delete[] XI2dep;
		delete[] U2dep;
		delete[] XII;
	}
	return 0;
	
}





