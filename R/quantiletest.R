quantiletest<-function(k,a,b,d,n,p,F,I,maxq)
{.C("main",as.integer(k),as.numeric(a),as.numeric(b),as.numeric(d),as.integer(n),as.integer(p),retour=as.numeric(F),as.integer(I),as.integer(maxq))}
