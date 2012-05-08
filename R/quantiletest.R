quantiletest<-function(k,a,b,c,n,p,F,I,maxq,sigma)
{.C("main",as.integer(k),as.numeric(t(a)),as.numeric(t(b)),as.numeric(t(c)),as.integer(n),as.integer(p),retour=as.numeric(F),as.integer(I),as.integer(maxq),as.numeric(sigma))}
