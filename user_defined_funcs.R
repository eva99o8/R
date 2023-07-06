## Visualization 

gen_diseaseplot=function(Y,start=1,end=ncol(Y),region_names){
  Time= seq(start,end,by=1)
  n_region= nrow(Y)
  return(tibble(Y = c(Y[,Time]),time= rep(Time,each = n_region),region = factor(rep(region_names,length(Time)),levels = region_names))%>%
    ggplot()+
    geom_line(aes(x=time,y= Y,group=region,color=region))+
    facet_wrap(~region))
}
