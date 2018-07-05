
 ########## neg log lik function for P(data | parameters)

## normal errors
nll<-function(data, predict, SD){
if(any(is.na(predict)) | any(is.nan(predict)) | min(predict) < 0){
return(NA)
}else{
return(-sum(stats::dnorm(x = data, mean = predict, sd = SD,  log=TRUE)))
}
}
