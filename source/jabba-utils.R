require(JABBA)
require(FLCore)
require(FLBRP)

require(plyr)
require(dplyr)


auxFn<-function(lag=0,obsE=0.3,sigma=TRUE,type="",...){
  
  args=list(...)
  
  auxilary      =NULL
  
  if (length(args)>0)
    if(names(args)%in%c("z", "f", "ffmsy", "bbmsy", "bk")){
      type     =names(args)[names(args)%in%c("z", "f", "ffmsy", "bbmsy", "bk")][1]
      auxilary=args[[type]]}
  
  list(
    auxiliary      =auxilary, 
    auxiliary.type =type,
    auxiliary.sigma=sigma, # estimated?
    auxiliary.obsE =obsE,   
    auxiliary.lag  =lag)    # lag effect between impact and Z pop structure
}


jabFn<-function(catch,prior,index=NULL,model="Pella_m",fix=NULL,...){
  
  ## priors
  r        =unlist(c(prior[,c("r")]))
  r.prior  =c(r,  0.3)

  psi      =unlist(c(prior[,c("ssb.minyr")]/prior[,c("b0")]))
  psi.prior=c(psi,0.3)
  
  shape    =unlist(c(prior[,c("shape")]))
  
  b.prior  =unlist(c(prior[,"ssb.maxyr"]/prior[,"bmsy"]))
  b.prior  =c(b.prior,1e-6,max(catch$year),"bbmsy") 
  
  f.prior  =unlist(c(prior[,"f.maxyr"]/prior[,"fmsy"]))
  f.prior  =c(f.prior,1e-6,max(catch$year),"ffmsy") 
  
  ## auxillary index
  aux=auxFn(...)
  
  args=c(list(scenario  ="",
              model.type=model,
              BmsyK     =shape,
               
              catch     =catch,
              cpue      =index,
               
              r.prior   =r.prior,
              psi.prior =psi.prior,
               
              verbose   =FALSE),aux)
  
  if (!is.null(fix))
     if (fix=="b") args=c(args,list(b.prior=b.prior)) else
       if (fix=="f") args=c(args,list(b.prior=f.prior))
          
  ## Fit with Catch + Index: Simple Fox with r = Fmsy
  input=try(do.call("build_jabba", args))
  
  if ("try-error"%in%is(input)) return(NULL)
  
  fit=try(fit_jabba(input,quickmcmc = T,verbose=F))
  
  if ("try-error"%in%is(fit)) return(NULL)
  
  list(input=input,fit=fit)}


hindJabba<-function (jbinput, fit, ni = NULL, nt = NULL, nb = NULL, nc = NULL, 
          quickmcmc = TRUE, init.values = TRUE, peels = 1:5, verbose = FALSE, save.all=TRUE){
  runs = as.list(peels)
  if (is.null(ni) & !quickmcmc) 
    ni = fit$settings$mcmc$ni
  if (is.null(nt) & !quickmcmc) 
    ni = fit$settings$mcmc$nt
  if (is.null(nb) & !quickmcmc) 
    ni = fit$settings$mcmc$nb
  if (is.null(nc) & !quickmcmc) 
    ni = fit$settings$mcmc$nc
  if (init.values) {
    Kin = fit$pars[1, 1]
    rin = fit$pars[2, 1]
    qin = fit$pars[(3:(2 + fit$settings$nq)), 1]
  }
  else {
    Kin = NULL
    rin = NULL
    qin = NULL
  }
  fithc = lapply(runs, function(x) {
    fit_jabba(jbinput, save.trj = TRUE, ni = ni, nt = nt, 
              nb = nb, nc = nc, init.values = init.values, init.K = Kin, 
              init.r = rin, init.q = qin, quickmcmc = T, peels = x, 
              verbose = verbose, save.all=save.all)
  })
  retro = c(list(fit), fithc)
  names(retro) = c(fit$scenario, paste0("-", max(fit$yr) - 
                                          peels + 1))
  retro = Map(function(x, y) {
    x$settings$scenario = y
    x
  }, x = retro, y = as.list(names(retro)))
  return(retro)
}

hindFn1<-function(catch,prior,index=NULL,model="Pella_m",fix=NULL,...){
  
  ## priors
  r        =unlist(c(prior[,c("r")]))
  psi      =unlist(c(prior[,c("ssb.minyr")]/prior[,c("b0")]))
  shape     =unlist(c(prior[,c("shape")]))
  r.prior  =c(r,  0.3)
  psi.prior=c(psi,0.3)
  
  b.prior  =unlist(c(prior[,"ssb.maxyr"]/prior[,"bmsy"]))
  b.prior  =c(b.prior,1e-6,max(catch$year),"bbmsy") 
  f.prior  =unlist(c(prior[,"f.maxyr"]/prior[,"fmsy"]))
  f.prior  =c(f.prior,1e-6,max(catch$year),"ffmsy") 
  
  ## auxillary index
  aux=auxFn(...)
  
  args=c(list(scenario  ="",
              model.type=model,
              BmsyK     =shape,
              
              catch     =catch,
              cpue      =index,
              
              r.prior   =r.prior,
              psi.prior =psi.prior,
              
              verbose   =FALSE),aux)
  
  if (!is.null(fix))
    if (fix=="b") args=c(args,list(b.prior=b.prior)) else
      if (fix=="f") args=c(args,list(b.prior=f.prior))
  
  ## Fit with Catch + Index: Simple Fox with r = Fmsy
  jbInput=try(do.call("build_jabba", args))
  
  if ("try-error"%in%is(jbInput)) return(NULL)
  
  jbFit=try(fit_jabba(jbInput,quickmcmc = T,verbose=F))
  
  if ("try-error"%in%is(jbFit)) return(NULL)

  hnd=try(hindcast_jabba(jbInput, jbFit, verbose=FALSE,save.all=TRUE))
  
  if ("try-error"%in%is(hnd)) return(NULL)
  
  list(input=jbInput,fit=jbFit,hind=hnd)}


source("C:/active/FLCandy/R/OMstats.R")


benchmarks<-function(x) {
  if ("logical"%in%is(attributes(x)$benchmark))
    return(FLPar(Fmsy=NA,Flim=NA,Fpa=NA,Blim=NA,Bpa=NA,Btrigger=NA))
  
  as(attributes(x)$benchmark,"FLPar")}

fishlife2lhPar<-function(x) {
  res=attributes(x)$fishlife
  
  if ("lm"%in%names(res))
    names(res)[seq(length(res))[(names(res)=="lm")]]="l50"
  
  res=FLPar(res,units="NA")
  
  lhPar(res[c("linf","k","l50","s")])}

priorFn<-function(x,nmin=0:2,nmax=0:2){
  .
  fmsy     =unlist(c(attributes(x)$benchmark["Fmsy"]))
  bmsy     =unlist(c(attributes(x)$eqsim["BMSY"]))
  b0       =unlist(c(attributes(x)$eqsim["B0"]))
  
  ssb.minyr=mean(ssb( x)[,1+nmin])
    f.minyr=mean(fbar(x)[,1+nmin])
    h.minyr=mean(hr(  x)[,1+nmin])
  ssb.maxyr=mean(ssb( x)[,dim(ssb( x))[2]-nmax])
    f.maxyr=mean(fbar(x)[,dim(fbar(x))[2]-nmax])
    h.maxyr=mean(hr(  x)[,dim(fbar(x))[2]-nmax])
  
  shape=bmsy/b0
  if (is.na(shape)) shape=0.4
  m    =NA
  r    =NA
  
  mi   =seq(0.01,2,0.001) 
  m    =(mi^(-1/(mi-1))-shape)^2
  m    =mi[m==min(m)]
  r    =(1-exp(-fmsy))*(m-1)/(1-m^-1)
  
  rtn=c(r=r,mpar=m,fmsy=fmsy,bmsy=bmsy,b0=b0,shape=shape,
        ssb.minyr=ssb.minyr,ssb.maxyr=ssb.maxyr,
          f.minyr=  f.minyr,  f.maxyr=  f.maxyr,
          h.minyr=  h.minyr,  h.maxyr=  h.maxyr)
  names(rtn)=c("r","mpar","fmsy","bmsy","b0","shape",
               "ssb.minyr","ssb.maxyr",
               "f.minyr",  "f.maxyr",
               "h.minyr",  "h.maxyr")
  
  rtn}

ebiomass<-function(object){
  sel  =harvest(object)
  wt   =catch.wt(object)%*%sel%/%fapex(sel)
  eb.wt=qmax(wt,0.000001)
  
  apply(eb.wt%*%stock.n(object),2:6,sum)}

dataFn<-function(object){
  
  model.frame(FLQuants(catch   =catch(object),
                       eb      =ebiomass(object),
                       ssb     =ssb(object),
                       p       =production(object),
                       f       =fbar(object),
                       h       =catch(object)/ebiomass(object),
                       m       =FLQuant(aaply(m(object)[ac(range(object)["minfbar"]:range(object)["maxfbar"])],2,mean),
                                        dimnames=dimnames(fbar(object)))),drop=TRUE)}

eqlFn<-function(x){
  
  sr=fmle(as.FLSR(x,model="bevholtSV"), fixed=list(s=attributes(x)$fishlife["s"],spr0=spr0(x)),
            control=list(silent=TRUE))
  
  brp(FLBRP(x,nyears=dim(x)[2],sr=list(model =bevholt()$model,
                                  params=ab(params(sr),model="bevholt"))))}

peFn<-function(stk,eq,stock=FLCore:::ssb){
  (stock(stk)%-%
     window(stock(stk)[,-1],end=dims(stk)$maxyear+1)-
     catch(stk)%+%sp(stk,eq,stock))%/%stock(stk)}


ldfFn<-function(object){
  
  ak =invALK(iter(par,1),cv=0.1,age=an(dimnames(object)$age),bin=1)
  lfd=lenSamples(catch.n(object),ak,n=5000)
  
  units(lfd)="cm"
  
  lfd}

tsFn<-function(x) ldply(names(x), function(.id) { 
  if (is.null(x[[.id]][["fit"]])) return(NULL)
  
  cbind(.id=.id,
        year=as.numeric(dimnames(x[[.id]][["fit"]][["timeseries"]][,"mu",])[[1]]),
        as.data.frame(x[[.id]][["fit"]][["timeseries"]][,"mu",]),
        x[[.id]][["fit"]][["refpts"]][1,c("k","bmsy","fmsy","msy")])})

kbFn<-function(x) ldply(names(x), function(.id) { 
  if (is.null(x[[.id]][["fit"]])) return(NULL)
  
  cbind(.id=.id,
        x[[.id]][["fit"]]$kobe,
        x[[.id]][["fit"]]$refpts_posterior,
        x[[.id]][["fit"]]$pars_posterior)[,-(2:4)]})


ldfFn<-function(object){
  
  ak =invALK(iter(par,1),cv=0.1,age=an(dimnames(object)$age),bin=1)
  lfd=lenSamples(catch.n(object),ak,n=5000)
  
  units(lfd)="cm"
  
  lfd}


priorPostFn<-function(x) {
  rtn=try({
    prior=data.frame(par  =dimnames(x$fit$pars)[[1]],
                     prior=x$fit$pars[,"Median"])
    post =melt(x$fit$pars_posterior)
    rtn=merge(prior,post,by.y="variable",by.x="par")
    names(rtn)[3]="posterior"
    rtn})
  
  if ("try-error"%in%is(rtn)) return(NULL)
  rtn}

maseBoot<-function(data, indices) {
  # select bootstrap sample
  d=data[indices, ] 
  return(mase(d$observed, d$predicted))}

if(FALSE){
# Load required packages
require(Metrics)
require(boot)

# Data is in a data frame  with columns:
# 'observed', 'predicted', 'sd'
df=data.frame(observed=rlnorm(12),predicted=rlnorm(12,0,0.2),sd=rep(0.1,12))

# Calculate MASE
val=mase(df$observed, df$predicted)

# Function to calculate MASE for bootstrapping

# Perform bootstrap with 1000 replicates
set.seed(123) # for reproducibility
bootRes<-boot(df, statistic=maseBoot, R=1000)

# Get 95% confidence intervals
ci     =boot.ci(bootRes, type="basic")
lower95=ci$basic[4] 
upper95=ci$basic[5]

cat("MASE:",   round(val, 3), "\n")
cat("95% CI:", round(lower95, 3), "-", round(upper95, 3), "\n")

rtn=ddply(subset(pp,par%in%c("K","r")), .(.id,par), with, 
            data.frame(posterior=median(posterior),
                       prior    =median(prior)))
rtn=merge(cast(rtn,.id~par,value="prior"),
            cast(rtn,.id~par,value="posterior"),
            by=".id")
  
ggplot(rtn)+geom_point(aes(r.x,K.y))
}




