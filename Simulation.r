#! /bin/Rscript

library("MASS")

mkdf <- function(cm,n){
  df <- mvrnorm(n,mu=rep(0,ncol(cm)),Sigma=cm)
  return(df)
}

preCorrection <- function(){
    qv <<- c(1:ceiling(RsnpNum*0.001),seq(ceiling(RsnpNum*0.001)+1,ceiling(RsnpNum*0.01),round(RsnpNum*0.00001)),seq(ceiling(RsnpNum*0.01)+1,ceiling(RsnpNum*0.999)-1,round(RsnpNum*0.0002)),ceiling(RsnpNum*0.999):RsnpNum)
}

GC <- function(pv,df){
  lamed <- qchisq(quantile(pv,0.5),df,lower.tail=F)/qchisq(0.5,df,lower.tail=F)
  pv <- pchisq(qchisq(pv,df,lower.tail=F)/lamed,df,lower.tail=F)
  return(pv)
}

preGPC <- function(){
    ae <<- -log10(ppoints(Esnp))
    indhp <<- ceiling(Esnp*0.0001):floor(Esnp*0.5)
    indhhhp <<- 20:(ceiling(Esnp*0.01)-1)
    indhhhhp <<- 1:(ceiling(Esnp*0.01)-1)
    indlp <<- (floor(Esnp*0.5)+1):Esnp

    ind151p <<- ceiling(Esnp*0.05):floor(Esnp*0.5)
    indhap <<- indh15p <<- 1:floor(Esnp*0.5)
    indhmp <<- ind11p <<- ceiling(Esnp*0.01):ceiling(Esnp*0.1)
    indhhp <<- indh1p <<- 1:ceiling(Esnp*0.1)

    ind252p <<- ceiling(Esnp*0.005):floor(Esnp*0.05)
    indh25p <<- 1:floor(Esnp*0.05)
    ind22p <<- ceiling(Esnp*0.001):ceiling(Esnp*0.01)
    indh2p <<- 1:ceiling(Esnp*0.01)

    ind353p <<- ceiling(Esnp*0.0005):floor(Esnp*0.005)
    indh35p <<- 1:floor(Esnp*0.005)
    ind33p <<- ceiling(Esnp*0.0001):ceiling(Esnp*0.001)
    indh3p <<- 1:ceiling(Esnp*0.001)
}

GPC1 <- function(p,ml){
    ind <- order(p,decreasing=F)
    ao <- -log10(p[ind])

    o = ao[ind151p]
    e = ae[ind151p]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    o = ao[indh15p]
    ao[indh15p] <- (o-min(o))/(lamea/ml)+min(o)
    o = ao[indh1p]
    e = ae[indh1p]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    ao[indh1p] <- (o-min(o))/(lamea/ml)+min(o)

    e = ae[indlp]
    ao[indlp] <- ae[indlp]

    p[ind] <- 10^(-ao)
    return(p)
}

GPC2 <- function(p,ml){
    ind <- order(p,decreasing=F)
    ao <- -log10(p[ind])

    o = ao[ind252p]
    e = ae[ind252p]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    o = ao[indh25p]
    ao[indh25p] <- (o-min(o))/(lamea/ml)+min(o)
    o = ao[indh2p]
    e = ae[indh2p]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    ao[indh2p] <- (o-min(o))/(lamea/ml)+min(o)

    p[ind] <- 10^(-ao)
    return(p)
}

GPC3 <- function(p,ml){
    ind <- order(p,decreasing=F)
    ao <- -log10(p[ind])

    o = ao[ind353p]
    e = ae[ind353p]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    o = ao[indh35p]
    ao[indh35p] <- (o-min(o))/(lamea/ml)+min(o)
    o = ao[indh3p]
    e = ae[indh3p]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    ao[indh3p] <- (o-min(o))/(lamea/ml)+min(o)

    p[ind] <- 10^(-ao)
    return(p)
}

GPC <- function(p,ml){
    ind <- order(p,decreasing=F)
    ao <- -log10(p[ind])

    o = ao[indhp]
    e = ae[indhp]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    o = ao[indhap]
    ao[indhap] <- (o-min(o))/(lamea/ml)+min(o)

    o = ao[indhhp]
    e = ae[indhhp]
    lamea <- sum((e-min(e))*(o-min(o)))/sum((e-min(e))*(e-min(e)))
    mo = ao[indhmp]
    oldmaxo <- max(mo)
    ao[indhmp] <- mo <- (mo-min(mo))/(lamea/ml)+min(mo)

    ho = ao[indhhhp]
    he = ae[indhhhp]
    hlamea <- sum((he-min(he))*(ho-min(ho)))/sum((he-min(he))*(he-min(he)))
    ho = ao[indhhhhp]
    ao[indhhhhp] <- (ho-min(ho))/(hlamea/ml)+min(ho)-(oldmaxo-max(mo))

    e = ae[indlp]
    ao[indlp] <- ae[indlp]

    p[ind] <- 10^(-ao)
    return(p)
}

corterm3all <- function(p,b,fc){
  stepw <- order(p)
  p <- p[stepw]
  b <- b[stepw]
  corm <- fc[stepw,stepw]
  temv <- which(b<0)
  corm[temv,] <- -corm[temv,]
  corm[,temv] <- -corm[,temv]
  corm <- 3.2630398097*corm + 0.7095678755*corm^2 + 0.0268257772*corm^3 + 0.0005732151*corm^4

  stat <- -2*log(p)
  curcov <- rep(0,phenoNum)
  for(i in 2:phenoNum){
    stat[i] <- stat[i]+stat[i-1]
    curcov[i] <- curcov[i-1]+sum(corm[1:(i-1),i])
  }

  V <- 2*curcov + 2*E
  scaf <- V/(2*E)
  dff <- 2*E^2/V
  tdff <<- tdff+dff
  cp <- pchisq(stat/scaf,dff,lower.tail=F)
  return(cp)
}

# main function
    rprofPath = 'D:/work/HPC419/0_CGWAStest/rprof/step2_rprof.out'
    Rprof(rprofPath)
    t1 <- Sys.time()
    Esnp <- snpNum <- 1e+05
    RsnpNum <- 110000
    ocorm <- as.matrix(read.csv("D:/work/HPC419/0_CGWAStest/Result/4StatsCor.csv",header=T))
    phenoNum <- ncol(ocorm)
    E <- 2*(1:phenoNum)
    simulateNum <- 20

    preGPC()
    preCorrection()
    setwd("D:/work/HPC419/0_CGWAStest/Tempdata/Simudata/Allp/")
    fminm <- matrix(0,ncol=phenoNum,nrow=simulateNum*length(qv))
    fgpcm <- matrix(0,ncol=5,nrow=simulateNum*length(qv))
    fdffm <- matrix(0,ncol=phenoNum,nrow=simulateNum)
    for(i in 1:simulateNum){
        resm <- matrix(0,ncol=phenoNum,nrow=snpNum)
        fillind <- (length(qv)*(i-1)+1):(length(qv)*i)
        dm <- signif(mkdf(ocorm,snpNum),6)
        pm <- pnorm(-abs(dm))*2
        tdff <- rep(0,phenoNum)
        for(j in 1:snpNum){
            resm[j,] <- corterm3all(pm[j,],dm[j,],ocorm)
        }
        fdffm[i,] <- tdff/snpNum
        fminm[fillind,] <- 10^(apply(log10(resm),2,quantile,probs=((qv-1)/(RsnpNum-1))))
        opv <- GC(resm[,phenoNum],fdffm[i,phenoNum])
        fgpcm[fillind,1] <- 10^(quantile(log10(GPC(opv,1)),(qv-1)/(RsnpNum-1)))
        gpcv <- GPC1(opv,1)
        fgpcm[fillind,2] <- 10^(quantile(log10(gpcv),(qv-1)/(RsnpNum-1)))
        gpcv <- GPC2(gpcv,1)
        fgpcm[fillind,3] <- 10^(quantile(log10(gpcv),(qv-1)/(RsnpNum-1)))
        gpcv <- GPC3(gpcv,1)
        fgpcm[fillind,4] <- 10^(quantile(log10(gpcv),(qv-1)/(RsnpNum-1)))
        unifp <- pchisq(rchisq(snpNum,fdffm[i,phenoNum]),fdffm[i,phenoNum],lower.tail=F)
        fgpcm[fillind,5] <- 10^(quantile(log10(unifp),(qv-1)/(RsnpNum-1)))
        print(paste0('finish ', i))
        print(difftime(Sys.time(), t1, units = 'secs'))
    }
    for(i in 1:phenoNum){
        write.table(t(matrix(fminm[,i],nrow=length(qv))),paste0("qutm_",i),row.names=F,col.names=F,quote=F)
    }
    for(i in 1:5){
        write.table(t(matrix(fgpcm[,i],nrow=length(qv))),paste0("qutm_",phenoNum+i),row.names=F,col.names=F,quote=F)
    }
    write.table(fdffm,"dfm",row.names=F,col.names=F,quote=F)
    Rprof(NULL)

