#refit.proctest=function(object,...){UseMethod("refit")}

refit.proctest_ord=function(object,Ynew,ordrenew,IT,var_nonselect,sigma,showresult,...)
{
	#print("refit.proctest_ord")
	n=nrow(object$data$X)
	p=ncol(object$data$X)

	if(missing(Ynew)){stop('Ynew is missing')}
	
	alpha=as.numeric(colnames(object$coefficients))  	     
	data=object$data$X
	Y=object$data$Y	
	
	p=ncol(data)
	n=nrow(data)
	if(length(Ynew)!=n){stop(" 'data' and 'Ynew' must have the same length ")}
	if(missing(showresult)){showresult=TRUE}
	if(missing(sigma)){sigma=0}
	if(missing(IT)){IT=2000}


if(missing(ordrenew))
{ordrenew=object$ordre	
	dataa=data
}else{
	if(length(ordrenew)<p)
	{for(i in 1:p)
		{if(sum(i==ordrenew)==0){ordrenew=c(ordrenew,i)}} #on complete l'ordre par les variables restantes
		}
	dataa=data[,ordrenew]
	}
	
ORDREBETA=ordrenew
		
dec=decompbaseortho(dataa)
nonind=dec$nonind
U=dec$U
if(p>(n+length(nonind)))
{nonind2=c(nonind,(n+1+length(nonind)):p)	
Uchap=U[,-nonind2]}else{if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}}
dim_X=ncol(Uchap)			#nombre de variables utiles
			


	aV2=array(0,c(length(alpha),dim_X-1,nrow(object$ordrebeta)+1)) 	
	aV2[,1:dim(object$quantile)[2],1:nrow(object$ordrebeta)]=object$quantile
	#print(I)
INDICE=matrix(0,nrow=2,ncol=nrow(object$ordrebeta))
#pour chaque ordre i on a INDICE[1,i] qui nous donne var_nonselect de l'ordre i et INDICE[2,i] qui nous donne le ktest max
	for(i in 1:nrow(object$ordrebeta))
	{
		j=1
		while(sum(aV2[,j,i]^2)==0) {j=j+1}
			INDICE[1,i]=j-1
		while(sum(aV2[,j,i]^2)!=0) {j=j+1
			if(j>ncol(object$quantile)){break}}
			INDICE[2,i]=j-1			
	}
						  
if(missing(var_nonselect)) 
{var_nonselect=INDICE[1,i]+1
}else{	
	if(var_nonselect<1){stop("var_nonselect has to be greater than 1")}}			

	#n=nrow(data)	
	#p=ncol(data)
	
	Y=Ynew
	ORDREBETA2=matrix(object$ordrebeta,ncol=p)
		
beta=lm(Y~Uchap-1)$coefficients #decomposition de Y2 dans la base orthonormal (X(1),..,X(p))
beta[-which(beta!=0)]=0

		#		----------------------------------------------------------------------------
		#			on test pour separer les variables en deux: les bonnes et les autres
		#		----------------------------------------------------------------------------
		
NBR=numeric(0) #on y met kchap
NBR_effect=numeric(0) #on y met kchap
relevant_variables=numeric(0) #on y met les bonnes variables
aV=matrix(0,length(alpha),min(n,dim_X)-1) #on y met tous les quantiles
indice=0
calcul=numeric(0)
for(alph in 1:length(alpha)) #boucle sur alpha
{
	ktest=var_nonselect-sum(nonind<=var_nonselect) # on commence par la première variable à selectionner
	
	nbr_test=var_nonselect

#indice=indice2[alph]
	T=1
	while((T>0)&&(dim_X>ktest+1))
		{
		if(showresult){print(paste("ktest=",ktest,"alpha=",alpha[alph]))}
		
		#on a indice de longueur compteur-1, aV_compteur avec les quantiles, var_select de dim compteur-1, et ORDREBETA2 
if(ktest>indice)
{abc=nrow(ORDREBETA2)#dim de indice/ORDREBETA2/var_select
	i=0
	I=numeric(0)
	TT=0
	#A=numeric(0)
	while((TT==0)&&(i<abc))#dim(ORDREBETA2)=abc*maxordre
	{i=i+1
	a=numeric(0)
	K=numeric(0)
	for(j in 1:nbr_test)
	{a=c(a,sum(ORDREBETA[j]==ORDREBETA2[i,1:nbr_test]))}
		if((sum(a)==nbr_test)&&(INDICE[2,i]>=ktest)&&(INDICE[1,i])<ktest){TT=1
#		K=c(K,ktest)
		I=c(I,i)}
	}

	if(length(I)!=0)
	{	aV[,ktest]=aV2[,ktest,I[length(I)]]#object$quantile[,,ktest]#get(paste("aV",I,sep="_"))[,,ktest]
		#blabla=c(blabla,0)#on met 0 si le quantile est deja calculé
		calcul=c(calcul,0)}else{
			print("blabla")
			#if(showresult){print(paste("ktest=",ktest))}
		quant=quantil_ord(n,dim_X,k=ktest,alpha,IT,sigma=sigma)#calcul du alpham
		
			#quant=quantileprocbol(XI_ord,k=ktest,alpha,IT,maxq=maxq,sigma=sigma)
		aV[,ktest]=quant$quantile
		calcul=c(calcul,1)
		#blabla=c(blabla,1)#on met 1 si on a calculé un quantile manquant
			#if(showresult){print(aV[,ktest])}
		}
indice=indice+1
}

alpha2=aV[alph,ktest]

		bb=numeric(0)
		for(m in 0:(log(min(n,dim_X)-ktest-1,2)-1))
		{
			if(sigma==0)
			{a=(n-ktest-2^m)*sum(beta[(ktest+1):(ktest+2^m)]^2)/(2^m*sum((Y-as.matrix(Uchap[,1:(ktest+2^m)])%*%beta[1:(ktest+2^m)])^2))
				F=qf(1-alpha2,2^m,n-ktest-2^m)
			}else{a=sum(beta[(ktest+1):(ktest+2^m)]^2)/sigma
			F=qchisq(1-alpha2,2^m)}
				b=a>F #1 si on doit rejeter le test, 0 sinon
				bb=c(bb,b)
				if(showresult)
				{	
				print(a)
				print(F)
				}
			}
		if(showresult){print(bb)}

		if(sum(bb)>0){T=1
			ktest=ktest+1
		
			nbr_test=nbr_test+1
			if(length(nonind)>0)
			{for(i in 1:length(nonind))
				if(sum(nonind==(nbr_test+1))==1){nbr_test=nbr_test+1}}
					
			if(showresult)
			{print(paste("number of variables:",nbr_test))
			print(paste("dimension of the space:",ktest))}	
		}else{
			T=0
			#print(ktest)
		} #on rejete Hk s'il y a au moins un test qui rejete
						
}
if(ktest==dim_X){k0=dim_X}else{k0=ktest}
							
NBR=c(NBR,nbr_test) #rÈsultat contenant le nombre de variables sÈlectionnÈes
NBR_effect=c(NBR_effect,k0)		
relevant_variables=rbind(relevant_variables,c(ORDREBETA[1:nbr_test],rep(0,NBR[1]-nbr_test)))
	#}else{
	#	OR=ORDREBETA[1:(k0+sum(nonind<=var_nonselect))]-1
	#	if(NBR[1]>1){relevant_variables=rbind(relevant_variables,c(OR[-1],rep(0,NBR[1]-(length(OR)))))}else{relevant_variables=rbind(relevant_variables,0)}}

}#fin boucle sur alpha
	
	if(showresult){print("relevant variables:")
	print(relevant_variables)}

if(sum(calcul)!=0)
{
#aV2[,,1:dim(object$quantile,1)[2]]=object$quantile
#aV2=array(0,c(length(alpha),maxqdep,dim_X)) #on y met tous les quantiles
#aV2[,,1:max(NBR+(var_nonselect==0)),nrow(object$ordrebeta)+1]=aV[,,(1:max(NBR+(var_nonselect==0)))]
aV2[,,nrow(object$ordrebeta)+1]=aV

	#aV2=aV2[,1:max(NBR_effect,dim(object$quantile)[2]),]		
if(length(alpha)==1){aV2=array(aV2[,1:max(NBR_effect,dim(object$quantile)[2]),],c(1,max(NBR_effect,dim(object$quantile)[2]),nrow(object$ordrebeta)+1))}else{
aV2=aV2[,1:max(NBR_effect,dim(object$quantile)[2]),]}

			
ORDREBETA=rbind(ORDREBETA2,ORDREBETA)

rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("ktest=",1:max(NBR_effect,dim(object$quantile)[2]))	
dimnames(aV2)[[3]]=paste("ordrebeta=",1:dim(aV2)[3])
}else{aV2=object$quantile
	ORDREBETA=object$ordrebeta}
#aV2=aV2[,,1:max(NBR+(var_nonselect==0))]

#fitted.values


Y.fitted=NULL
coefficients=matrix(0,nrow=p,ncol=length(alpha))

#if(intercept)
#{
for(i in 1:length(alpha))
{reg=lm(Y~data[,relevant_variables[i,]]-1)
	coefficients[relevant_variables[i,],i]=reg$coefficients
	reg$coefficients[-which(reg$coefficients!=0)]=0
	Y.fitted=cbind(Y.fitted,data[,relevant_variables[i,]]%*%reg$coefficients)
}
#}else{
#	for(i in 1:length(alpha))
#{reg=lm(Y~data[,c(1,1+relevant_variables[i,])]-1)
#	coefficients[c(1,1+relevant_variables[i,]),i]=reg$coefficients
#	print(coefficients)
#	Y.fitted=cbind(Y.fitted,data[,c(1,1+relevant_variables[i,])]%*%reg$coefficients)
#}
	
	
#}
colnames(Y.fitted)=alpha
rownames(coefficients)=c("intercept",paste("X",2:p,sep=""))
colnames(coefficients)=alpha

out=list(data=list(X=data,Y=Y),coefficients=coefficients,relevant_var=relevant_variables,fitted.values=Y.fitted,ordre=object$ordre,ordrebeta=ORDREBETA,kchap=NBR,quantile=aV2,call=match.call(),call.old=object$call)

out
structure(out,class="proctest_ord")	

		
}#fin refit proc_ord
