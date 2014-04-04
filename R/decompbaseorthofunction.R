decompbaseortho=function(data)
{
	#safety
	if(ncol(data)<2){stop("only one dimension, it is useless to run this algorithm.. ")}
	   
U=svd(data[,1])$u
p=ncol(data)
n=nrow(data)
nonind=numeric(0)
for(k in 1:(p-1))
{	
	m=0
	u=U
	a=dim(u)[2]#base ortho de Vk

#on decompose Xk+1,..,Xk+2^m (=dataak+1,..,dataak+2^m) dans cette base :
x=matrix(0,a,1)
	for(i in 1:a)
	{
		x[i,1]=sum(data[,k+1]*u[,i]) #en colonne on a la variable qui est decomposÈ, et en ligne sa decomp dans la base
	}
X=u%*%x#la decomp de chaque variable dans Vk

#on peut donc calculer l'orthogonal:

X_=data[,k+1]-X
	if(max(abs(X_))<10^-10)		#on met un seuil pour dire que la colonne devient 0 (à cause de la précision des calculs)
	{X_=0
	u1=0}else{u1=svd(X_)$u}	
	U=cbind(U,u1)
	if(sum(t(u1)%*%u1)==0){nonind=c(nonind,k+1)}
	if((dim(U)[2]-length(nonind))==n){break}
	}
a=dim(U)[2]
U=cbind(U,matrix(0,n,p-a))
return(list(U=U,nonind=nonind))
}