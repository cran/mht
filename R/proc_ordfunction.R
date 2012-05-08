#proc_ord=function(data,...){UseMethod("proc_ord")}
#proc_ord.default=function(data,Y,ordre,var_nonselect,alpha,IT,sigma,showresult,...)
proc_ord=function(data,Y,ordre,var_nonselect,alpha,IT,sigma,showresult)
{
#-----------------------------------	
#	da=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=XI*beta
#	ordre=ordre qu'on veut donner à la matrice data, peut etre partiel, a defaut on considere la matrice deja ordonnée
#	var_nonselect= le nombre de variables qu'on ne veut pas séléctionner, ce sont les premières variables de "ordre"
#	alpha=erreur de première espèce du test
#	IT=nombre de simulations pour le calcul du quantile
#	sigma= si on travaille à variance connue, pas la même stat de test
#	showresult=affiche les resultats des tests
#-----------------------------------------

	
# verifications de depart
p=ncol(data)
n=nrow(data)
if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}
if(missing(ordre)){ordre=1:p}
if(missing(alpha)){alpha=c(0.1,0.05)}
if(missing(showresult)){showresult=TRUE}
if(missing(IT)){IT=20000}
if(missing(sigma)){sigma=0}

ordreinit=ordre
if(length(ordre)<p)
{for(i in 1:p)
	{if(sum(i==ordre)==0){ordre=c(ordre,i)}} #on complete l'ordre par les variables restantes
	}	
		
#		-------------------------------------
#			on scale la matrice de départ
#		-------------------------------------

#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute 
	if(sum(data[,1])==n){data=cbind(data[,1],scale(data[,-1])/sqrt(n-1))
	data[which(is.na(data))]=0
		intercept=TRUE
		dataa=data[,ordre]
		}else{
			message("intercept has been added")
			data=scale(data)/sqrt(n-1)
			data[which(is.na(data))]=0
			data=cbind(1,data)
			p=p+1
			intercept=FALSE
			ordre=c(1,1+ordre)
			dataa=data[,ordre]#pour compenser l'ajout de l'intercept
		}
ordre=matrix(ordre,nrow=1)
if(missing(var_nonselect)){var_nonselect=1
}else{if(!intercept){var_nonselect=var_nonselect+1}
	if(var_nonselect<1){stop("var_nonselect has to be greater than 1")}}			


#		-------------------------------------
#		  on décompose dans une base ortho
#		-------------------------------------
dec=decompbaseortho(dataa) #U, nonind

U=dec$U
nonind=dec$nonind

if(p>(n+length(nonind)))
{nonind2=c(nonind,(n+1+length(nonind)):p)
	UUU=U[,-nonind2]
}else{
	if(length(nonind)==0){UUU=U}else{UUU=U[,-nonind]}
	}
dim_X=ncol(UUU)		#nombre de variables utiles

Uchap=UUU
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

for(alph in 1:length(alpha)) #boucle sur alpha
{

	ktest=var_nonselect-sum(nonind<=var_nonselect) # on commence par la première variable à selectionner
	
	nbr_test=var_nonselect
	
	T=1
	while((T>0)&&(dim_X>ktest+1))
	{
		if(showresult){print(paste("ktest=",ktest,"alpha=",alpha[alph]))}
		if(ktest>indice)
		{ 	quant=quantil_ord(n,dim_X,k=ktest,alpha,IT,sigma=sigma)#calcul du alpham
			aV[,ktest]=quant$quantile
		indice=indice+1}
		
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
	#if(intercept)
	#{
	relevant_variables=rbind(relevant_variables,c(ordre[1:nbr_test],rep(0,NBR[1]-nbr_test)))
	#}else{relevant_variables=rbind(relevant_variables,c(ordre[1:(k0-1)],rep(0,NBR[1]-k0)))}

}#fin boucle sur alpha

	if(showresult){print("relevant variables:")
	print(relevant_variables)}

if(length(alpha)==1){aV2=matrix(aV[,1:max(NBR_effect)],nrow=1)}else{
aV2=aV[,1:max(NBR_effect)]}

rownames(aV)=paste("alpha=",alpha)
rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("ktest=",1:max(NBR_effect))	
print(aV2)
#fitted.values
Y.fitted=NULL
coefficients=matrix(0,ncol=length(alpha),nrow=p)
for(i in 1:length(alpha))
{reg=lm(Y~dataa[,relevant_variables[i,]]-1)
	coefficients[relevant_variables[i,],i]=reg$coefficients
	reg$coefficients[-which(reg$coefficients!=0)]=0
	Y.fitted=cbind(Y.fitted,dataa[,relevant_variables[i,]]%*%reg$coefficients)

}

colnames(Y.fitted)=alpha
rownames(coefficients)=c("intercept",paste("X",2:p,sep=""))
colnames(coefficients)=alpha
	
out=list(data=list(X=dataa,Y=Y),coefficients=coefficients,relevant_var=relevant_variables,fitted.values=Y.fitted,ordre=ordreinit,ordrebeta=ordre,kchap=NBR,quantile=aV2,call=match.call())

out
structure(out,class="proctest_ord")
	
}