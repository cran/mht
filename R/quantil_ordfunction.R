quantil_ord=function(n,p,k,alpha,IT)
{
	if(missing(alpha)){alpha=c(0.1,0.05)}
	if(missing(IT)){IT=20000}
		
#simulation pour trouver \alpham
FF=numeric(0)

for(i in 1:IT)
{
eps=as.matrix(rnorm(n))

F=numeric(0)
for(m in 0:(log(min(n,p)-k-1,2)-1))
{
A=(n-k-2^m)*sum(eps[(k+1):(k+2^m)]^2)/(2^m*sum((eps[(k+2^m+1):n])^2))
F[m+1]=1-pf(A,2^m,n-k-2^m)
}
FF[i]=min(F)
}

alph=numeric(0)
for(i in 1:length(alpha))
{alph=c(alph,quantile(FF,alpha[i]))}
	

return(list(quantile=alph))
}
