bolasso<-function(data,Y,mu,m,probaseuil,var_nonselect)
{################################################################################
###########################   Bolasso   #######################################
################################################################################

if(missing(probaseuil)){probaseuil=0.5}
if(missing(m)){m=100}

n=nrow(data)
p=ncol(data)

if(length(Y)!=n){stop(" 'data' and 'Y' must have the same length ")}

if(missing(var_nonselect)){var_nonselect=1
	}else{
		if(sum(data[,1])!=n){var_nonselect=var_nonselect+1}
	}
if(var_nonselect<0){stop("var_nonselect has to be nonnegative")}
				
#		-------------------------------------
#			on scale la matrice de départ
#		-------------------------------------

#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute et on rajoute 1 dans var_nonselect s'il n'etait pas manquant
if(sum(data[,1])==n){data=cbind(data[,1],scale(data[,-1])/sqrt(n-1))
		data[which(is.na(data))]=0
			}else{
				data=scale(data)/sqrt(n-1)
				data[which(is.na(data))]=0
				data=cbind(1,data)
				p=p+1
			}


dimu=length(mu)
#on fait une premiere estimation de notre modèle, avec un lasso
lasso1=glmnet(data,Y,family="gaussian",alpha=1,lambda=mu,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)))

mat=matrix(0,n)
for(i in 1:dimu)
{
	if(var_nonselect>0)
	{beta_ind=c(1:var_nonselect,which(lasso1$beta[,i]!=0))
		}else{beta_ind=which(lasso1$beta[,i]!=0)}

	reg=lm(Y~data[,beta_ind]-1)
	beta0=reg$coefficients
	beta0[-which(beta0!=0)]=0
	mat0=as.matrix(data[,beta_ind])%*%beta0
	mat=cbind(mat,mat0)
	}
mat=mat[,-1]

eps0=as.matrix(Y)%*%matrix(1,1,dimu)-mat
eps1=scale(eps0,center=TRUE,scale=FALSE)


compteur=matrix(0,p,dimu)
for (j in 1:m)
{a=runif(n,0,n)
	b=ceiling(a)
	Y33=as.matrix(mat+eps1[b,])	#on bootstrap les residus
	
	betaboot2=matrix(0,p,dimu)
	for(i in 1:dimu)
	{mu1=mu[i]
	lasso1=glmnet(data,Y33[,i],family="gaussian",alpha=1,lambda=mu1,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)))
	
		betaboot2[,i]=(lasso1$beta[,1]!=0)
		if(var_nonselect>0){betaboot2[1:var_nonselect,i]=1}
    }
    compteur=compteur+betaboot2 #on recence les coefficients non nuls

}
compteur2=as.matrix(compteur)


probavariable=compteur2/m
beta_ind=(probavariable>probaseuil)


return(list(data=data,ind=beta_ind,frequency=probavariable))
}