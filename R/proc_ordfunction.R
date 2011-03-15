proc_ord=function(data,Y,ordre,alpha,IT,showresult)
{
#-----------------------------------	
#	da=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=XI*beta
#	ordre=ordre qu'on veut donner à la matrice data, peut etre partiel, a defaut on considere la matrice deja ordonnée
#	alpha=erreur de première espèce du test
#	IT=nombre de simulations pour le calcul du quantile
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
			data=scale(data)/sqrt(n-1)
			data[which(is.na(data))]=0
			data=cbind(1,data)
			p=p+1
			intercept=FALSE
			dataa=data[,c(1,1+ordre)]#pour compenser l'ajout de l'intercept
		}

	

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
relevant_variables=numeric(0) #on y met les bonnes variables
aV=matrix(0,length(alpha),min(n,dim_X)-1) #on y met tous les quantiles
indice=0

for(alph in 1:length(alpha)) #boucle sur alpha
{

	ktest=1 #on ne selectionne pas l'intercept
	
	T=1
	while((T>0)&&(dim_X>ktest+1))
	{
		print(paste("ktest=",ktest,"alpha=",alpha[alph]))
		if(ktest>indice)
		{ 	quant=quantil_ord(n,dim_X,k=ktest,alpha,IT)#calcul du alpham
			aV[,ktest]=quant$quantile
		indice=indice+1}
		
		alpha2=aV[alph,ktest]

		bb=numeric(0)
		for(m in 0:(log(min(n,dim_X)-ktest-1,2)-1))
		{
			a=(n-ktest-2^m)*sum(beta[(ktest+1):(ktest+2^m)]^2)/(2^m*sum((Y-as.matrix(Uchap[,1:(ktest+2^m)])%*%beta[1:(ktest+2^m)])^2))
			F=qf(1-alpha2,2^m,n-ktest-2^m)
			b=a>F #1 si on doit rejeter le test, 0 sinon
			bb=c(bb,b)
			if(showresult)
			{	
				print(a)
				print(F)
			}
		}
		if(showresult){print(bb)}

		if(sum(bb)>0){T=1}else{T=0} #on rejete Hk s'il y a au moins un test qui rejete

		if(T==0){print(ktest)}
		ktest=ktest+1
		}
if(ktest+1==dim_X){k0=dim_X-1}else{k0=ktest-1}

k0=k0+sum(nonind<k0)#on reprend les indices qu'on a laissé de coté parce que n'apportant aucune dimension
NBR=c(NBR,k0) #rÈsultat contenant le nombre de variables sÈlectionnÈes

if(intercept)
{relevant_variables=rbind(relevant_variables,c(ordre[1:k0],rep(0,NBR[1]-k0)))
}else{relevant_variables=rbind(relevant_variables,c(ordre[1:(k0-1)],rep(0,NBR[1]-k0)))}

}#fin boucle sur alpha

aV2=aV[,1:max(NBR)]
rownames(aV)=paste("alpha=",alpha)
rownames(aV2)=paste("alpha=",alpha)
colnames(aV2)=paste("ktest=",1:max(NBR))	
	
return(list(data=dataa,relevant_var=relevant_variables,kchap=NBR,quantile=aV2))
	
}