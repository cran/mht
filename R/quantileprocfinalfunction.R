quantileprocbol=function(data,k,alpha,IT,maxq,sigma)
{
n=nrow(data)
p=ncol(data)
	
	if(missing(alpha)){alpha=c(0.1,0.05)}
	if(missing(IT)){IT=1000}
	if(missing(maxq)){maxq=log(min(n,p)-1,2)}
	if(missing(sigma)){sigma=0}

b=0
nbrprob=0
seuil=2
while((sum(b==0)!=0)&(nbrprob<seuil)) # on recommence si à la fin on n'a pas assez de simulations pour obtenir une estimation du quantile
{
ijk=k


if(ijk>0) #on doit construire une base de Vk
{
	liste1=1:ijk
	dataa=as.matrix(data[,liste1]) 	#on doit construire une base ortho de dataa
	XI2dep=as.matrix(data[,-liste1])
		
	if(ijk>1)
	{
		dec=decompbaseortho(dataa)
		U=dec$U
			}else{U=svd(dataa)$u}
	U2dep=U
}	

sup=floor(maxq) #nombre max de tests que l'on fait
#print(sup)
F=matrix(0,IT*sup)
if(k==0){XI2dep=0
	U2dep=0}
a=quantiletest(k,data,XI2dep,U2dep,n,p,F,IT,maxq-1,sigma=sigma)	
F=t(matrix(a$retour,sup,IT))
#print(F[,sup])
long=200
lg_alpha=length(alpha)
bV=array(0,c(long+1,sup,lg_alpha))
aV=matrix(0,lg_alpha,sup)

alphatest=matrix(0,long+1,lg_alpha)

alphatest[1,]=alpha/log(n,2)
for(j in 1:lg_alpha)
{
a=(alpha[j]-alpha[j]/(log(n,2)))/long
for(i in 1:long)
{
alphatest[i+1,j]=alphatest[i,j]+a}
}
#alphatest

FF=matrix(0,long+1,lg_alpha)
for(m in 0:(maxq-1))
{for(j in 1:lg_alpha)
	{for(i in 1:(long+1))
		{
		FF[i,j]=quantile(F[,m+1],1-alphatest[i,j]) #quantile de U_{k,m}
		}
	bV[,m+1,j]=FF[,j] #chaque colonne=quantile, le numero de colonne est le m, la troisieme dim est le alpha
	}
}

FFF=matrix(0,long+1,lg_alpha)
for(i in 1:IT)
{
FF=matrix(0,long+1,lg_alpha)
	for(m in 0:(maxq-1))
	{for(j in 1:lg_alpha)
		{
		FF[,j]=FF[,j]+(F[i,m+1]>=bV[,m+1,j])	#=0 si la ieme iteration n'est pas rejet», >1 sinon
		}
	}

#print(F)
for(j in 1:lg_alpha)
{FFF[,j]=FFF[,j]+(FF[,j]>0)}
}

b=numeric(0)
for(j in 1:lg_alpha)
{ind=which(FFF[,j]<=alpha[j]*IT)
a=ind[length(ind)]
#print(a)
#b=c(b,a)
if(length(a)==0){b=c(b,0)}else{b=c(b,a)}
}
			
	if(sum(b==0)!=0){nbrprob=nbrprob+1}#on incremente de 1 si on recommence

}#fin while

if(nbrprob==seuil){print(paste("number of simulations too low to estimate the ",alpha[which(b==0)],"-quantile, please increase IT"))}

#break

if(nbrprob<seuil)
{for(j in 1:lg_alpha)
{
ind=which(FFF[,j]<=alpha[j]*IT)
a=ind[length(ind)]
	aV[j,]=bV[a,,j]
}
}
return(list(quantile=aV,nbrprob=nbrprob))
}#fin function
