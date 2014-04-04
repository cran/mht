bolasso<-function(data,Y,mu,m,probaseuil,penalty.factor,random)
{################################################################################
###########################   Bolasso   #######################################
################################################################################

if(missing(probaseuil)){probaseuil=1}
if(missing(m)){m=100}

ntot=nrow(data)
p=ncol(data)

if(length(Y)!=ntot){stop(" 'data' and 'Y' must have the same length ")}

			
#		-------------------------------------
#			on scale la matrice de départ
#		-------------------------------------

#si la matrice de depart ne contient pas l'intercept (colonne de 1) on la rajoute et on rajoute 1 dans penalty.factor s'il n'etait pas manquant
if(sum(data[,1])==ntot){data=cbind(data[,1],scale(data[,-1])*sqrt(ntot)/sqrt(ntot-1))
		data[which(is.na(data))]=0
	intercept=TRUE
			}else{
				data=scale(data)*sqrt(ntot)/sqrt(ntot-1)
				data[which(is.na(data))]=0
				data=cbind(1,data)
				intercept=FALSE
				p=p+1
			}
			
if(missing(penalty.factor)){penalty.factor=c(0,rep(1,p-1))
	}else{
		if(!intercept){penalty.factor=c(0,penalty.factor)}
	}
	

dimu=length(mu)
#on fait une premiere estimation de notre modèle, avec un lasso
lasso1=glmnet(data,Y,alpha=1,lambda=mu,penalty.factor=penalty.factor)

mat=NULL


for(i in 1:dimu)
{
	if(penalty.factor[1]==0)
	{beta_ind=c(1,which(lasso1$beta[,i]!=0))
		}else{beta_ind=which(lasso1$beta[,i]!=0)}

	reg=lm(Y~data[,beta_ind]-1)
	beta0=reg$coefficients
	beta0[-which(beta0!=0)]=0
	mat0=as.matrix(data[,beta_ind])%*%beta0
	mat=cbind(mat,mat0)
	}

eps0=as.matrix(Y)%*%matrix(1,1,dimu)-mat
eps1=scale(eps0,center=TRUE,scale=FALSE)

if(missing(random))
{a=matrix(runif(ntot*m,0,ntot),nrow=ntot)
	random=ceiling(a)}
b=random

compteur=matrix(0,p,dimu)
for (j in 1:m)
{
	Y33=as.matrix(mat+eps1[b[,j],])	#on bootstrap les residus
	
	betaboot2=matrix(0,p,dimu)
	for(i in 1:dimu)
	{mu1=mu[i]
	lasso1=glmnet(data,Y33[,i],alpha=1,lambda=mu1,penalty.factor=penalty.factor)
	
		betaboot2[,i]=(lasso1$beta[,1]!=0)
		betaboot2[which(penalty.factor==0),i]=1
		#if(var_nonselect>0){betaboot2[1:var_nonselect,i]=1}
    }
    compteur=compteur+betaboot2 #on recence les coefficients non nuls

}


compteur2=as.matrix(compteur)

probavariable=compteur2/m
colnames(probavariable)=mu
beta_ind=(probavariable>=probaseuil)

out=list(data=data,ind=beta_ind,frequency=probavariable,call=match.call())		
		out
		structure(out,class="bolasso")
		
}