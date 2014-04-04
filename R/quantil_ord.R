quantil_ord=function(n,p,k,alpha,IT,sigma)
{
	if(missing(alpha)){alpha=c(0.1,0.05)}
	if(missing(IT)){IT=20000}
	if(missing(sigma)){sigma=0}
	if(sigma<0){stop("sigma<0? Really?")}
	if(missing(n)||(missing(p))||(missing(k))){stop("something is missing, either n or p or k")}

#simulation pour trouver \alpham
FF=numeric(0)


if(sigma==0)
{
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
}else{
for(i in 1:IT)
{
	eps=as.matrix(rnorm(n,sd=sqrt(sigma)))
	F=numeric(0)
	for(m in 0:(log(min(n,p)-k-1,2)-1))
	{
	A=sum(eps[(k+1):(k+2^m)]^2)/sigma
	F[m+1]=1-pchisq(A,2^m)
	}
	FF[i]=min(F)
}
}

alph=numeric(0)
for(i in 1:length(alpha))
{alph=c(alph,quantile(FF,alpha[i]))}
	

return(list(quantile=alph))
}
