step5 <- function(){
    if(HPC){
        setwd(paste0(CGPath,"Tempdata/Rowdata/"))
        system("rm x*")
        for(n in 1:Sepblock){
            system(paste0("rm Cp",n,".r"))
        }
    }
    preCorrection()
    
    setwd(paste0(CGPath,"Result/"))
    system(paste0("cat 1.o > All.sortsp"))
    #system(paste0("rm 1.o"))
    system(paste0("cat 1.lo > All.labsp"))
    #system(paste0("rm 1.lo"))
    system(paste0("cat 1.ro > All.sp"))
    #system(paste0("rm 1.ro"))
    print(paste0("Merge GWAS Data 1"))
    for(i in 2:Sepblock){
        system(paste0("cat ",i,".o >> All.sortsp"))
        #system(paste0("rm ",i,".o"))
        system(paste0("cat ",i,".lo >> All.labsp"))
        #system(paste0("rm ",i,".lo"))
        system(paste0("cat ",i,".ro >> All.sp"))
        #system(paste0("rm ",i,".ro"))
        print(paste0("Merge GWAS Data ",i))
    }

    pm <- matrix(0,ncol=length(Traitname),nrow=Rsnpnum)
    for(i in 1:Sepblock){
        tem <- as.matrix(fread(paste0(i,".res"),header=F))
        pm[(round(Rsnpnum/Sepblock)*(i-1)+1):(round(Rsnpnum/Sepblock)*(i-1)+nrow(tem)),] <- tem
        #system(paste0("rm ",i,".res"))
        print(paste0(i," CGWAS Data Extracted"))
    }

    basecf <- as.numeric(read.table(paste0(CGPath,"Tempdata/Simudata/Minp/",length(Traitname)+1),header=F,stringsAsFactors=F)[,1])
    for(i in 1:length(Traitname)){
        curcf <- as.numeric(read.table(paste0(CGPath,"Tempdata/Simudata/Minp/",i),header=F,stringsAsFactors=F)[,1])
        baseper <- basecf/curcf
        simuper <- matrix(0,1,Rsnpnum)
        simuper[1,qv] <- baseper
        temq <- qv-c(0,qv[-length(qv)])
        temb <- c(baseper,0)-c(0,baseper)
        seqtack <- which(temq!=1)

        for(ii in seqtack){
            simuper[,(qv[ii-1]+1):(qv[ii]-1)] <- temb[ii]/temq[ii]*c(1:(temq[ii]-1))+baseper[ii-1]
        }
        tind <- order(pm[,i])
        pm[tind,i] <- pm[tind,i]*simuper
        print(paste0(i," CGWAS Data Corrected"))
        if(i==1){
            sptind <- tind
            spsimuper <- simuper
        }
    }

    pm[pm>1] <- 1
    pm <- signif(pm,4)
    mincp <- apply(pm,1,min)
    write.table(as.data.frame(mincp),"Min.corrcp",row.names=F,col.names=F,quote=F)
    print(paste0("Derived CGWAS Minp"))
    sugghitind2 <- which(mincp<St)
    HJind2 <- which(mincp<Ht)

    spo <- as.matrix(fread("All.sortsp",header=F))
    #system(paste0("rm All.sortsp"))
    
    baseseq <- c((1:length(Traitname))-1)/(length(Traitname)-1)
    for(i in 1:length(sptind)){
        spo[sptind[i],] <- spo[sptind[i],]*spsimuper[i]/(baseseq*(spsimuper[i]-1)+1)
    }
    spo[spo>1] <- 1
    spo <- signif(spo,4)
    fwrite(as.data.frame(spo),"All.corrsortsp",row.names=F,col.names=F,quote=F,sep=" ")
    print(paste0("corrGWAS Result Saved"))
    minsp <- apply(spo,1,min)
    write.table(as.data.frame(minsp),"Min.corrsp",row.names=F,col.names=F,quote=F)
    print(paste0("Derived GWAS Minp"))
    sugghitind <- which(as.numeric(minsp)<St)
    HJind <- which(minsp<Ht)

    sugghit3 <- union(sugghitind,sugghitind2)
    HJ3 <- union(HJind,HJind2)
    sugghit <- matrix(0,length(sugghit3)*length(Traitname),3)
    HJM <- matrix(0,length(HJ3),2)
    sugghit[,1] <- rep(sugghit3,each=length(Traitname))
    sugghit[,2] <- as.vector(t(pm[sugghit3,]))
    sugghit[,3] <- rep(c(1:length(Traitname)),length(sugghit3))
    HJM[,1] <- as.numeric(mincp)[HJ3]
    HJM[,2] <- as.numeric(minsp)[HJ3]
    HJpm <- pm[HJ3,]
    fwrite(as.data.frame(pm),paste0("All.corrcp"),row.names=F,col.names=F,quote=F,sep=" ")
    print(paste0("corrCGWAS Result Saved"))
    rm(pm)

    Sind <- as.data.frame(fread("SnpIndex",header=T))
    temres <- cbind(Sind[as.numeric(sugghit[,1]),],sugghit)
    temHJ <- cbind(Sind[HJ3,],HJM)
    temres <- temres[order(temres[,5]),]
    temres <- temres[order(temres[,2]),]
    temres <- temres[order(temres[,1]),]
    temHJ <- temHJ[order(temHJ[,2]),]
    temHJ <- temHJ[order(temHJ[,1]),]
    colnames(temres)[4:6] <- c("SNPid","CorrCP","CombNum")
    colnames(temHJ)[4:5] <- c("CP","BP")
    write.csv(temres,"CgwasRawSuggHits.csv",row.names=F,quote=F)
    write.table(as.data.table(temHJ),"HJ_man.txt",row.names=F,col.names=F,quote=F)
    print(paste0("CGWAS Suggestive Hits Saved"))

    sspo <- spo[sugghit3,]
    HJspo <- spo[HJ3,]
    rm(spo)
    for(i in 1:nrow(sspo)){
        for(ii in (length(Traitname)-1):1){
            if(sspo[i,ii]>sspo[i,ii+1]){
                sspo[i,ii] <- sspo[i,ii+1]
            }
        }
    }
    spsignum <- c()
    cpsignum <- c()
    for(i in 1:nrow(HJspo)){
        for(ii in (length(Traitname)-1):1){
            if(HJspo[i,ii]>HJspo[i,ii+1]){
                HJspo[i,ii] <- HJspo[i,ii+1]
            }
        }
        spsignum <- c(spsignum,sum(HJspo[i,]<=0.05))
        cpsignum <- c(cpsignum,which(HJpm[i,]==min(HJpm[i,]))[1])
    }
    spl <- as.matrix(fread("All.labsp",header=F))
    sspl <- spl[sugghit3,]
    HJspl <- spl[HJ3,]
    rm(spl)
    sugghit <- matrix(0,length(sugghit3)*length(Traitname),3)
    sugghit[,1] <- rep(sugghit3,each=length(Traitname))
    sugghit[,2] <- as.vector(t(sspo[,]))
    sugghit[,3] <- as.vector(t(sspl[,]))
    
    temres <- cbind(Sind[as.numeric(sugghit[,1]),],sugghit)
    temHJall <- cbind(Sind[HJ3,],cpsignum,spsignum,HJpm,HJspl,HJspo)
    temres <- temres[order(temres[,2]),]
    temres <- temres[order(temres[,1]),]
    temHJall <- temHJall[order(temHJall[,2]),]
    temHJall <- temHJall[order(temHJall[,1]),]
    temHJall <- temHJall[,-c(1:2)]
    colnames(temres)[4:6] <- c("SNPid","SP","Trait")
    write.csv(as.matrix(temres),"GwasRawSuggHits.csv",row.names=F,quote=F)
    write.table(as.data.table(temHJall),"HJ_bar.txt",row.names=F,col.names=F,quote=F)
    print(paste0("GWAS Suggestive Hits Saved"))
    print(paste0("Step5 Completed"))
}