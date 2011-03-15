procbol=function(data,Y,var_nonselect,alpha,maxordre,choix_ordre=c("bolasso","pval","pval_hd"),m,maxit,showordre,IT,maxq,showtest,showresult)
{
#-----------------------------------	
#	data=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=data*beta
#	var_nonselect= le nombre de variables qu'on ne veut pas séléctionner, ce sont les premières colonnes dans data
#	alpha=erreur de première espèce du test
#	maxordre= nombre de variables max qu'on ordonne
#	choix_ordre=avec quelle méthode on ordonne (pval, pval_hd ou bolasso)
#	m=nombre d'iteration lasso pour le bolasso
#	maxit=nombre maximum d'itération de l'algorithme d'ordre que l'on fait
#	showordre=affiche l'ordre au fur et à mesure
#	IT=nombre de simulations pour le calcul du quantile
#	maxq=nombre max d'hypotheses alternative testees
#	showtest=affiche le nombre de mu tester dans l'algorithme d'ordre
#	showresult=affiche les resultats des tests
#-----------------------------------------
	
n=nrow(data)	
p=ncol(data)
	#safety
	if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}
	
	if(missing(m)){m=100}
	if(missing(maxordre)){maxordre=min(n-1,p-1)}
	if(missing(alpha)){alpha=c(0.1,0.05)}
	if(missing(choix_ordre)){choix_ordre="bolasso"}
	if(missing(IT)){IT=1000}
	if(missing(maxit)){maxit=5}#nombre de fois ou on redemarre l'algo
	if(missing(maxq)){maxq=log(min(n,p)-1,2)}
	if(missing(showtest)){showtest=FALSE}
	if(missing(showordre)){showordre=TRUE}
	if(missing(showresult)){showresult=TRUE}



TOT=numeric(0)
if(missing(var_nonselect)){var_nonselect=1
	}else{
		if(sum(data[,1])!=n){var_nonselect=var_nonselect+1
			}
	}			
if(var_nonselect<0){stop("var_nonselect has to be nonnegative")}


#		-------------------------------------
#			on scale la matrice de départ
#		-------------------------------------

#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute et on rajoute 1 dans var_nonselect s'il n'etait pas manquant
if(sum(data[,1])==n){data=cbind(data[,1],scale(data[,-1])/sqrt(n-1))
		data[which(is.na(data))]=0
	intercept=TRUE
			}else{
				data=scale(data)/sqrt(n-1)
				data[which(is.na(data))]=0
				data=cbind(1,data)
				intercept=FALSE
				p=p+1
			}


maxordre=min(n-1,p-1,maxordre+!intercept)
maxq=min(maxq,log(min(n,p)-1,2))



#		-------------------------------------
#			on ordonne les variables
#		-------------------------------------

if(choix_ordre=="pval")
{
	if(p<n) 	#on calcule pval avec une seule regression comportant toutes les variables
	{reg=lm(Y~data-1)
		a=summary(reg)$coefficients[,4]
		b=sort(a,index.return=TRUE)
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
		b=sort(a,index.return=TRUE)
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
	b=sort(a,index.return=TRUE)
	ORDREBETA=b$ix
	XI_ord=data[,b$ix] #on a ainsi les XI ordonnÈes
	if(showordre){print(b$ix)}
}


if(choix_ordre=="bolasso")
{
	ordre=dyadiqueordre(data,Y,m,maxordre,var_nonselect,maxit,showtest=showtest,showordre=showordre)# donne l'ordre (dans ordre) et le nombre de fois ou l'algo a redemarré (dans prob)	
	if(ordre$prob>=maxit)
	{stop(paste("can't achieve to order", maxordre, "variables in ",maxit, "iterations"))}
		
	b=ordre$ordre
	if(showordre){print(b)}

	for(i in 1:p)
{if(sum(i==b)==0){b=c(b,i)}} #on complete l'ordre par les variables restantes
ORDREBETA=b
XI_ord=data[,b] #on a ainsi les XI ordonnées dans XI_ord
}


#		-------------------------------------
#		  on décompose dans une base ortho
#		-------------------------------------
dec=decompbaseortho(XI_ord)#nous donne une base ortho de X(1),..X(p) ou U[,i] est l'apport de Xi dans la decomp en base
#U est une base orthonormal de X(1),...,X(p).
#ex: nonind=4 signifie que la variable X(4) est dans la base V(3)=vect(X(1),X(2),X(3))

		
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
relevant_variables=numeric(0) #on y met les bonnes variables
indice=var_nonselect-1-sum(nonind<=var_nonselect) #on fait les calculs de quantiles pour tous les alpha en même temps, donc pas besoin de le refaire pour chaque avec indice
for(alph in 1:length(alpha)) #boucle sur alpha
{
	ktest=var_nonselect-sum(nonind<=var_nonselect) # on commence par la première variable à selectionner
	

	T=1
	while((T>0)&&(dim_X>ktest+1))
		{
		print(paste("ktest=",ktest,"alpha=",alpha[alph]))
		maxq=min(log(min(n,dim_X)-ktest-1,2),maxq)
		
		if(ktest>indice)
		{quant=quantileprocbol(XI_ord,k=ktest,alpha,IT,maxq=maxq)
			aV[,1:maxq,ktest]=quant$quantile
		indice=indice+1}
		if(quant$nbrprob==3){break}

		#test S	
		bb=numeric(0)
	
		for(m in 0:(maxq-1))
			{
			a=(n-ktest-2^m)*sum(beta[(ktest+1):(ktest+2^m)]^2)/(2^m*sum((Y-as.matrix(Uchap[,1:(ktest+2^m)])%*%beta[1:(ktest+2^m)])^2))
			b=a>aV[alph,m+1,ktest]#F #1 si on doit rejeter le test, 0 sinon
			if(showresult)
			{print(a)
			print(aV[alph,m+1,ktest])
				}
			bb=c(bb,b) #on met tous les tests de Hk
			}
			if(showresult){print(bb)}

		if(sum(bb)>0){T=1}else{T=0} #on rejete Hk s'il y a au moins un test qui rejete

		ktest=ktest+1
		}
	if(ktest+1==dim_X){k0=dim_X-1}else{k0=ktest-1}
	NBR=c(NBR,(k0+sum(nonind<=var_nonselect))) #resultat contenant le nombre de variables selectionnees
if(k0>(var_nonselect-sum(nonind<=var_nonselect)))
{	if(intercept)
	{relevant_variables=rbind(relevant_variables,c(ORDREBETA[1:(k0+sum(nonind<=var_nonselect))],rep(0,NBR[1]-(k0+sum(nonind<=var_nonselect)))))
	}else{relevant_variables=rbind(relevant_variables,c(ORDREBETA[2:(k0+sum(nonind<=var_nonselect))]-1,rep(0,NBR[1]-(k0+sum(nonind<=var_nonselect)))))}
}else{relevant_variables=rbind(relevant_variables,ORDREBETA[1:var_nonselect])}
	
}#fin boucle sur alpha
	
if(max(NBR)>1)
{aV2=aV[,,(1:max(NBR))]
rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("Hk,",0:(maxqdep-1))
dimnames(aV2)[[3]]=paste("ktest=",1:max(NBR))	
}else{aV2=as.matrix(aV[,,1])
	rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("Hk,",0:(maxqdep-1))}
	
return(list(data=data,relevant_var=relevant_variables,ordre=ordre$ordre,ordrebeta=ORDREBETA,kchap=NBR,quantile=aV2))
}#fin procbol
