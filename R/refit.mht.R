refit.mht=function(object,Ynew,var_nonselect,sigma,maxordre,ordre=c("bolasso","pval","pval_hd","FR"),m,show,IT,maxq,...)
{
	n=nrow(object$data$X)
	p=ncol(object$data$X)
	
	if(missing(m)){m=100}
	if(missing(maxordre)){maxordre=min(n/2-1,p/2-1)}
	if(missing(ordre)){ordre="bolasso"}
	if(missing(IT)){IT=1000}
	if(missing(maxq)){maxq=log(min(n,p)-1,2)}
	if(missing(sigma)){sigma=0}
	if(missing(show)){show=c(1,0,1)}
		showordre=show[1]
		showtest=show[2]
		showresult=show[3]
	if(missing(Ynew)){stop('Ynew is missing')}
	choix_ordre=ordre
	
alpha=as.numeric(colnames(object$coefficients))  	     
data=object$data$X
Y=object$data$Y

dec=decompbaseortho(data)
nonind=dec$nonind
U=dec$U
if(p>(n+length(nonind)))
{nonind2=c(nonind,(n+1+length(nonind)):p)	
	Uchap=U[,-nonind2]}else{if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}}
dim_X=ncol(Uchap)			#nombre de variables utiles

	aV2=array(0,c(length(alpha),dim(object$quantile)[2],dim_X,nrow(object$ordrebeta)+1)) 	
	aV2[,,1:dim(object$quantile)[3],1:nrow(object$ordrebeta)]=object$quantile
	
if(missing(var_nonselect)) var_nonselect=1
i=var_nonselect
while(sum(aV2[,,i,1]^2)==0) 
{var_nonselect=var_nonselect+1
i=i+1
if(i>dim(object$quantile)[3]) stop('bug')}


#safety
if(length(Ynew)!=n){stop(paste(" 'data' and 'Ynew' must have the same length:",n))}
	
	Y=Ynew
	ORDREBETA2=matrix(object$ordrebeta,ncol=p)
	
	indice2=NULL
	for(i in 1:nrow(object$ordrebeta))
	{j=var_nonselect
	while(sum(aV2[,,j,i]^2)!=0) 
	{j=j+1
		#	print(j)
	}
	j=j-1
	
	indice2=cbind(indice2,j)
	}
		

#		-------------------------------------
#			on ordonne les variables
#		-------------------------------------
if(choix_ordre=="pval")
{
	if(p<n) 	#on calcule pval avec une seule regression comportant toutes les variables
	{reg=lm(Y~data-1)
		a=summary(reg)$coefficients[,4]
		a[1]=0	#on ne selectionne pas l'intercept
		b=sort(a,index.return=TRUE)
		ORDRE=b$ix[1:maxordre]
		ORDREBETA=b$ix
		XI_ord=data[,b$ix] #on a ainsi les XI ordonnÈes
		if(showordre){print(b$ix)}
		}else{		#on calcule pval avec FDR2: une regression pour chaque variable

		print("p>n, the order 'pval' is not possible, 'pval_hd' is used instead")
		a=numeric(0)
		for(i in 1:p)
		{
			reg=lm(Y~data[,i]-1)
			c=summary(reg)$coefficients[,4]
			a=c(a,c)
		}
		a[1]=0
		b=sort(a,index.return=TRUE)
		ORDRE=b$ix[1:maxordre]
		ORDREBETA=b$ix
		XI_ord=data[,b$ix] #on a ainsi les XI ordonnÈes
		if(showordre){print(b$ix)}
		} 
}
if(choix_ordre=="pval_hd")
{a=numeric(0)
	for(i in 1:p)
	{
		reg=lm(Y~data[,i]-1)
		c=summary(reg)$coefficients[,4]
		a=c(a,c)
	}
	a[1]=0
	b=sort(a,index.return=TRUE)
	ORDRE=b$ix[1:maxordre]
	ORDREBETA=b$ix
	XI_ord=data[,b$ix] #on a ainsi les XI ordonnÈes
	if(showordre){print(b$ix)}
}
if(choix_ordre=="bolasso")
{
	ordre=dyadiqueordre(data,Y,m,maxordre,var_nonselect,showtest=showtest,showordre=showordre)# donne l'ordre (dans ordre) et le nombre de fois ou l'algo a redemarré (dans prob)	
	b=ordre$ordre
	ORDRE=b
	if(showordre){print(b)}

	for(i in 1:p)
{if(sum(i==b)==0){b=c(b,i)}} #on complete l'ordre par les variables restantes
ORDREBETA=b
XI_ord=data[,b] #on a ainsi les XI ordonnées dans XI_ord
}

if(choix_ordre=="FR")
{
    ordre=1
    ind.left=2:p
    while(length(ordre)<min(n-1,p,maxordre))
    {
        res=NULL
        for(i in 1:length(ind.left))
        {
            X=data[,c(ordre,ind.left[i])]
            I=diag(1,n)
            H=X%*%solve(t(X)%*%X)%*%t(X)
            res[i]=t(Y)%*%(I-H)%*%Y	#same as reg=lm(Y~XI[,VRAINDICE]-1); sum(reg$residuals^2)
        }
        a=which.min(res)
        ordre=c(ordre,ind.left[a])
        ind.left=ind.left[-a]
        
    }
    b=ORDRE=ordre
    for(i in 1:p)
    {if(sum(i==b)==0){b=c(b,i)}} #on complete l'ordre par les variables restantes
    ORDREBETA=matrix(b,nrow=1)
    XI_ord=data[,b] #on a ainsi les XI ordonnées dans XI_ord
}


#		-------------------------------------
#		  on décompose dans une base ortho
#		-------------------------------------
dec=decompbaseortho(XI_ord)#nous donne une base ortho de X(1),..X(p) ou U[,i] est l'apport de Xi dans la decomp en base	
#on rajoute dans nonind2 les dernieres variables, celle qui n'ont pas d'utilités puisque dans Rn
nonind=dec$nonind
U=dec$U
if(p>(n+length(nonind)))
{nonind2=c(nonind,(n+1+length(nonind)):p)	
	Uchap=U[,-nonind2]
}else{
	if(length(nonind)==0){Uchap=U}else{Uchap=U[,-nonind]}
		}
dim_X=ncol(Uchap)			#nombre de variables utiles


beta=lm(Y~Uchap-1)$coefficients #decomposition de Y2 dans la base orthonormal (X(1),..,X(p))
beta[-which(beta!=0)]=0

#		----------------------------------------------------------------------------
#			on test pour separer les variables en deux: les bonnes et les autres
#		----------------------------------------------------------------------------

maxqdep=min(log(min(n,dim_X)-1,2),maxq)
maxq=maxqdep
		

aV=array(0,c(length(alpha),maxq,dim_X)) #on y met tous les quantiles
NBR=numeric(0) #on y met kchap
NBR_effect=numeric(0) #on y met kchap
relevant_variables=numeric(0) #on y met les bonnes variables
indice=var_nonselect-1-sum(nonind<=var_nonselect) #on fait les calculs de quantiles pour tous les alpha en même temps, donc pas besoin de le refaire pour chaque avec indice
calcul=numeric(0)
for(alph in 1:length(alpha)) #boucle sur alpha
{
	ktest=var_nonselect-sum(nonind<=var_nonselect) # on commence par la première variable à selectionner
	
	nbr_test=var_nonselect

	T=1
	while((T>0)&&(dim_X>ktest+1))
		{
		if(showresult){print(paste("ktest=",ktest,"alpha=",alpha[alph]))}
		maxq=min(log(min(n,dim_X)-ktest-1,2),maxqdep)
		
		#on a indice de longueur compteur-1, aV_compteur avec les quantiles, var_select de dim compteur-1, et ORDREBETA2 
if(ktest>indice)
{abc=nrow(ORDREBETA2)#dim de indice/ORDREBETA2/var_select
	i=0
	I=numeric(0)
	TT=0
	while((TT==0)&&(i<abc))#dim(ORDREBETA2)=abc*maxordre
	{i=i+1
	a=numeric(0)
	K=numeric(0)
	for(j in 1:nbr_test)
	{a=c(a,sum(ORDREBETA[j]==ORDREBETA2[i,1:nbr_test]))}
	if((sum(a)==nbr_test)&&(indice2[i]>=ktest)){TT=1
#		K=c(K,ktest)
		I=c(I,i)}
	}

	if(length(I)!=0)
	{	aV[,,ktest]=aV2[,,ktest,I[length(I)]]#object$quantile[,,ktest]#get(paste("aV",I,sep="_"))[,,ktest]
		#blabla=c(blabla,0)#on met 0 si le quantile est deja calculé
		calcul=c(calcul,0)}else{
		if(showresult){print(paste("ktest=",ktest))}
		quant=quantilemht(XI_ord,k=ktest,alpha,IT,maxq=maxq,sigma=sigma)
		aV[,1:maxq,ktest+(var_nonselect==0)]=quant$quantile
		#blabla=c(blabla,1)#on met 1 si on a calculé un quantile manquant
		if(showresult){print(aV[,,ktest])
			calcul=c(calcul,1)}
		}
indice=indice+1
}

		#test S	
		bb=numeric(0)
	
		for(m in 0:(maxq-1))
			{
				if(sigma==0){
					a=(n-ktest-2^m)*sum(beta[(ktest+1):(ktest+2^m)]^2)/(2^m*sum((Y-as.matrix(Uchap[,1:(ktest+2^m)])%*%beta[1:(ktest+2^m)])^2))
				}else{a=sum(beta[(ktest+1):(ktest+2^m)]^2)/sigma}
			b=a>aV[alph,m+1,ktest+(var_nonselect==0)]#F #1 si on doit rejeter le test, 0 sinon
			if(showresult)
			{print(a)
			print(aV[alph,m+1,ktest+(var_nonselect==0)])
				}
			bb=c(bb,b) #on met tous les tests de Hk
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

}#fin boucle sur alpha
	
	if(showresult){print("relevant variables:")
	print(relevant_variables)}



if(sum(calcul)!=0)
{
aV2[,,,nrow(object$ordrebeta)+1]=aV

	aV2=aV2[,,1:max(NBR_effect,dim(object$quantile)[3]),]		
				
ORDREBETA=rbind(ORDREBETA2,ORDREBETA)

rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("Hk,",0:(maxqdep-1))
dimnames(aV2)[[3]]=paste("ktest=",(var_nonselect!=0):(max(NBR_effect,dim(object$quantile)[3])-(var_nonselect==0)))
dimnames(aV2)[[4]]=paste("ordrebeta=",1:dim(aV2)[4])
}else{aV2=object$quantile
	ORDREBETA=object$ordrebeta}

#fitted.values
Y.fitted=NULL
coefficients=matrix(0,nrow=p,ncol=length(alpha))

for(i in 1:length(alpha))
{reg=lm(Y~data[,relevant_variables[i,]]-1)
	coefficients[relevant_variables[i,],i]=reg$coefficients
	reg$coefficients[-which(reg$coefficients!=0)]=0
	Y.fitted=cbind(Y.fitted,data[,relevant_variables[i,],drop=FALSE]%*%reg$coefficients)
}

colnames(Y.fitted)=alpha
rownames(coefficients)=c("intercept",paste("beta",1:(p-1),sep=""))
colnames(coefficients)=alpha

out=list(data=list(X=data,Y=Y),coefficients=coefficients,relevant_var=relevant_variables,fitted.values=Y.fitted,ordre=ORDRE,ordrebeta=ORDREBETA,kchap=NBR,quantile=aV2,call=match.call(),call.old=object$call)

out
structure(out,class="mht")

		
}#fin refit procbol
