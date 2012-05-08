dyadiqueordre=function(data,Y,m,maxordre,var_nonselect,maxit,showtest,showordre)
{
	
#-----------------------------------	
#	data=matrice data, the first column should be 1, for the intercept
#	Y=observation	Y=XI*beta
#	m=nombre d'iteration lasso pour le bolasso
#	maxordre= nombre de variables max qu'on ordonne
#	var_nonselect= le nombre de variables qu'on ne veut pas séléctionner, ce sont les premières colonnes dans XI
#	maxit=nombre maximum d'itération de l'algorithme d'ordre que l'on fait
#	showtest=affiche le nombre de mu tester 
#	showordre=affiche l'ordre au fur et à mesure
#-----------------------------------------

p=ncol(data)
n=nrow(data)

if(missing(maxordre)){maxordre=min(n/2-1,p/2-1)}
if(missing(maxit)){maxit=5}
if(missing(showtest)){showtest=FALSE}
if(missing(showordre)){showordre=TRUE}
if(missing(m)){m=100}

			
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
	
if(missing(var_nonselect)){var_nonselect=1
}else{
	if(!intercept){var_nonselect=var_nonselect+1}
	}
		
				

numberproblem=0
problem=1
while((problem==1)&&(numberproblem<maxit)) #on recommence l'algorithme maxit fois au plus
{
problem=0

MU=numeric(0) #on va y mettre tous les mu que l'on teste
COMPTEUR2=numeric(0) #on va y mettre tous les compteurs que l'on teste

#on initialise le decoupage par le deuxième mu d'un lasso classique (le premier étant très grand)
lasso1=glmnet(data,Y,family="gaussian",alpha=1,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)))
	
mu1=lasso1$lambda[2]
muinf=0 #correspond au mu de gauche
musup=mu1 #=mulasso #correspond au mu de droite


if(var_nonselect>0){ordre=1:var_nonselect}else{ordre=numeric(0)}
nonordonne=numeric(0)
mu=musup
mumil=musup #correspond au mu du milieu, celui qu'on teste pour savoir si on regarde a gauche, dans ce cas musup devient mumil, ou si on regarde a droite, dans ce cas muinf devient mumil.

while(length(ordre)<maxordre)#on ordonne les variables jusqu'à en avoir maxordre
{

muinf=0
musup=mu1
mumil=musup
longordo=0
toto=0
while(longordo!=1)
{

if(sum(MU==mu)==0)#si on a pas encore testé le mu
{
bol=bolasso(data,Y,mu,m,probaseuil=1,penalty.factor=c(rep(0,var_nonselect),rep(1,p-var_nonselect)))
compteur2=bol$frequency

MU=c(MU,mu)#on met tous les mu qu'on a deja testé
COMPTEUR2=cbind(COMPTEUR2,compteur2)#on met les compteur2 de chaque mu testé

#un compteur qui affiche le nombre de mu testé
if(showtest)
{print(dim(COMPTEUR2))}

}else{a=which(MU==mu) 			#si on a deja testé le mu
	compteur2=COMPTEUR2[,a]}	#on evite de refaire un calcul deja fait


nonordonne=which(compteur2==1)#nous donne les indices des variables selectionnées

#on regarde maintenant lesquelles étaient deja ordonnées
ordonne=numeric(0)
longnon=length(nonordonne)
if(longnon>0)
{	for(i in 1:longnon)
	{if(sum(nonordonne[i]==ordre)==0)
		{ordonne=c(ordonne,nonordonne[i])}#on recence les variables séléctionnées qui n'ont pas encore été ordonnées
	}
	longordo=length(ordonne) # nombre de difference entre ordre et nonordonne
}else{longordo=0}

#si il n'y a pas de nouvelles variables à ordonner, on diminue la pénalité mu
if(longordo==0){
	musup=mumil
	mumil=(musup+muinf)/2
	mu=mumil}
	
#si il y a plus qu'une unique variable à ordonner, on augmente la pénalité mu pour n'en trouver qu'une seule	
if(longordo>1){muinf=mumil
	mumil=(muinf+musup)/2
	mu=mumil}

#on met une sécurité
if((musup-muinf)<10e-15){
	#print("break")
	problem=1
	break #nous sort du while s'il y a un problème
	}

}#fin while, on a donc une seule variable dans ordonne, c'est la variable suivante 

ordre=c(ordre,ordonne)

if(showordre){print(ordre)}

if((musup-muinf)<10e-15){
	if(showtest){print("break")}
	break #nous sort de l'algorithme s'il y a un problème
	}
}


numberproblem=numberproblem+1*(problem==1) #recence le nombre de problèmes

}#fin while numberproblem<maxit
#print(paste("number of iteration(s):",numberproblem+1))
#print(numberproblem+1)

return(list(data_used=bol$data,ordre=ordre,prob=numberproblem,mu=MU,compteur=COMPTEUR2))
}#fin function