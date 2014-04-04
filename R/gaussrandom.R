gaussrandomtest<-function(x,y)
{.C("gaussrand",as.numeric(x),as.integer(y))}

gaussrandom=function()
{
    y=1:10
    out=gaussrandomtest(y,10)
return(out)
}

