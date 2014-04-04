predict.mht.order=function(object,newx,...)
{
	 
	 if(missing(newx)){out=object$fitted.values
	  	
	  }else{
	  	    
          n=nrow(object$data$X)
          p=ncol(object$data$X)
          
		  if(sum(newx[,1])==n){newx=cbind(newx[,1],scale(newx[,-1])/sqrt(n-1))
			  newx[which(is.na(newx))]=0
			  intercept=TRUE
		  }else{
			  message("intercept has been added")
			  
			  newx=scale(newx)/sqrt(n-1)
			  newx[which(is.na(newx))]=0
			  newx=cbind(1,newx)
			  intercept=FALSE
			  p=p+1
		  }
		  
	alpha=as.numeric(colnames(object$coefficients))  	     
	out=NULL
	for(i in 1:length(alpha))
	{
	out=cbind(out,newx%*%object$coefficients[,i])
	}
	colnames(out)=paste("alpha=",alpha)	

	}
out
}
