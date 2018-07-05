### neg log lik function for P(data | parameters)

## 

poisson
nll<-function(data,predict){
	
if(any(is.na(predict)) | any(is.nan(predict)) | min(predict) < 0){
		
return(NA)
	
}else{

return(- sum (stats::dpois(data, lambda=predict, log=TRUE)))
	
}
}




## poisson

#nll<-function(data,predict){

#- sum(dpois(data, lambda=predict, log=TRUE))

#}
