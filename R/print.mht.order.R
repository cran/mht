print.mht.order <-
function(x,...) {
    # 0. part
    print(x$call)

    # 2a. part
    alpha=as.numeric(colnames(x$coefficients))  	     
    cat("","\n")
    for(i in 1:length(alpha))
{    cat("Results for alpha=",alpha[i],"\n")
    cat("Number of relevant variables:",x$kchap[i],"\n")
  	cat("which are:",x$relevant_var[i,],"\n\n")

}
}

