library(plyr); library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
JS<-function(x){length(x[!is.na(x)])}
header<-function(x, name  ) {
write.table(x, name, quote = F, col.names = T, row.names = F, sep = " " )	}
noheader<-function(x, name  ) {
write.table(x, name, quote = F, col.names = F, row.names = F, sep = " " )	}

## Sampling recombination breakpoints
breakPoints <- function(pms,pme,rpbMap){
  rce  <- rbinom(length(rpbMap),size=1,prob=rpbMap)
  nbrk <- sum(rce)
  if(nbrk>0){
    l   <- which(rce==1)
    brk <- sort( runif(n=nbrk,min=pms[l],max=pme[l]) )
  }else{
    brk <- NULL
  }
  return(brk)
}

## Recombine
recomb <- function(x,y,posChip,pms,pme,rpbMap){
  brks <- breakPoints(pms,pme,rpbMap)
  nbks <- length(brks)
  if(nbks>0){
    boundaries <- unique( c(0,brks,max(c(brks,posChip))) )
  }else{
    boundaries <- c(0,max(posChip))
  }
  nsegs  <- length(boundaries)-1
  offset <- rbinom(n=1,size=1,prob=0.5)
  ((1:nsegs) + offset)%%2
  
  z <- rep(NA,length(x))
  for(j in 1:nsegs){
    segment_j    <- which(posChip>=boundaries[j] & posChip<=boundaries[j+1])
    w            <- (j+offset)%%2
    if(w==1){
      z[segment_j] <- x[segment_j]
    }else{
      z[segment_j] <- y[segment_j]
    }
  }
  return(z)
}

FullSibs  <- function(posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,chr){
  p <- length(posChip)
  ## P and M
  p1 <- rep("p1",p)
  p2 <- rep("p2",p)
  m1 <- rep("m1",p)
  m2 <- rep("m2",p)
  
  ## Sib 1
  p12a <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  m12a <- recomb(m1,m2,posChip,pms,pme,rpbMap)
  
  ## Sib B
  p12b <- recomb(p1,p2,posChip,pms,pme,rpbMap)
  m12b <- recomb(m1,m2,posChip,pms,pme,rpbMap)
  
  ## Check IBD
  ## IBD 1
  SIBD1 <- 0
  SIBD2 <- 0
  segs  <- NULL
  
  comp  <- (p12a==p12b) & (m12a==m12b)
  nsg2  <- length(table(diff(which(comp))))

  if(nsg2>0){
    ecomp <- c(0,comp,0)
    epos  <- c(posChip[1],posChip,posChip[p])/1e6
    strs  <- epos[which(diff(ecomp)==+1)]
    ends  <- epos[which(diff(ecomp)==-1)]
    segs  <- rbind(segs,cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs,IBD=2))
    SIBD2 <- SIBD2 + sum(ends-strs)
  }

  #comp  <- (p12a==p12b) | (m12a==m12b)
  comp  <-  ( (p12a==p12b) | (m12a==m12b) ) & (!((p12a==p12b) & (m12a==m12b)))
  nsg1  <- length(table(diff(which(comp))))

  if(nsg1>0){ 
    ecomp  <- c(0,comp,0)
    epos   <- c(posChip[1],posChip,posChip[p])/1e6
    strs   <- epos[which(diff(ecomp)==+1)]
    ends   <- epos[which(diff(ecomp)==-1)]
    segs   <- rbind(segs,cbind(CHR=chr,START=strs,END=ends,LENGTHinMB=ends-strs,IBD=1))
    SIBD1  <- SIBD1 + sum(ends-strs) 
  }
  
  ## Simulate genotypes
  simGeno <- function(o12a,o12b){
    oa  <- rep(NA,p)
    ob  <- rep(NA,p)
    for(i in c("m1","m2","p1","p2")){
      li <- which(o12a==i)
      if(length(li)>0){
        oa[li] <- get(paste0(i,"h"))[li]
      }
      li <- which(o12b==i)
      if(length(li)>0){
        ob[li] <- get(paste0(i,"h"))[li]
      }
    }
    go <- oa + ob
  }
  geno_sibA <- simGeno(p12a,m12a)
  geno_sibB <- simGeno(p12b,m12b)
  
  ## IBD measures
  CL    <- (posChip[p]/1e6)
  FIBD1 <- SIBD1 / CL
  FIBD2 <- SIBD2 / CL

  nsg <- nrow(segs)
  if(length(nsg)==0){
    nsg <- 0
  }else{
   if(nsg>0){
     segs[order(segs[,"START"]),]
   }
  }
  ## Return
  return( list(genotype1=geno_sibA,
               genotype2=geno_sibB,
               segments=segs,
               fstats=c(CHR=chr,nSEG=nsg,SIBD1inMB=SIBD1,SIBD2inMB=SIBD2,Fibd1=FIBD1,Fibd2=FIBD2,CL=CL)
             )
  )
}

simulateGenotype <- function(indexParent1a,indexParent1b,
                             indexParent2a,indexParent2b,
                             Chromosome,verbose=FALSE){
  chr <- Chromosome
  ## Read genomic map
  if(verbose){
    cat(paste0("\tRead genomic map from chromosome ",ifelse(chr<10,paste0("0",chr),chr),".\n"))
  }
  
  mapfile <- paste0(dirMap,"/genetic_map_chr",chr,"_combined_b37.txt")
  map     <- read.table(mapfile,header=TRUE)
  rpbMap  <- diff(map[,3]) * 0.01
  pms     <- map[1:(nrow(map)-1),1]
  pme     <- map[2:nrow(map),1]
  
  ## Loading haplotypes
  if(verbose){
    cat("\tLoading haplotypes.\n")
  }
  load(paste0(dirHaps,"/chrom",chr,".haps.RData"))
  H12      <- cbind(haps$H1,haps$H2)
  posChip  <- as.numeric( gsub(paste0("chr",chr,":"),"",haps$ids,fixed=TRUE) )
  
  ## Get haplotypes
  p1h     <- H12[,indexParent1a]
  p2h     <- H12[,indexParent1b]
  m1h     <- H12[,indexParent2a]
  m2h     <- H12[,indexParent2b]
  
  results <- FullSibs(posChip,pms,pme,rpbMap,p1h,p2h,m1h,m2h,chr)
  return(results)
}

setwd("/gpfs1/scratch/30days/uqjsidor/projects/SibPairs/SIM/01_data_bcf_map")
dirMap  <- "/QRISdata/Q2970/within-family/ToolsForSimulations/bcftools_maps"
dirHaps <- "/QRISdata/Q2970/within-family/ToolsForSimulations/pool-haps"

## Things required to run this script
## Maps (folder)
## Haplotypes (folder) 
## map file (file.map)
## bim file (file.bim) corresponding to that map file

options       <- commandArgs(trailingOnly = TRUE)
prefix        <- options[1]
iPairStart    <- as.numeric(options[2])
iPairEnd      <- as.numeric(options[3])
pedfilename   <- paste0("ped_files/",prefix,"_", iPairStart,".ped")
RDatafilename <- paste0("RData/",prefix,"_", iPairStart, ".RData")
bimfile       <- "/QRISdata/Q2970/within-family/ToolsForSimulations/geno-trios-qced/combined.bim"
prune         <- fread("/QRISdata/Q2970/within-family/ToolsForSimulations/Linkage_simulations_map/INDEP_maf_above_0.1_bcf_map.prune.in", h=F)
bim           <- read.table(bimfile,stringsAsFactors=FALSE)
nsnps         <- table(bim[,1])
grps          <- lapply(1:22,function(chr){
  tmp <- bim[which(bim[,1]==chr),]
  A1  <- tmp[,5]
  A2  <- tmp[,6]
  grp <- cbind(paste0(A1,"\t",A1),paste0(A1,"\t",A2),paste0(A2,"\t",A2))
})

grps_pr_index <- lapply(1:22,function(chr){
			  tmp <- bim[which(bim[,1]==chr),]
			    ind <- which(tmp[,2]%in%prune$V1)
			    return(ind)
})
cat("# Simulation of Sib pairs...\n")
cat(paste0("# Output file: ",pedfilename,".\n"))
cat(paste0("# Pair ID range: ",iPairStart," - ",iPairEnd,".\n"))

indexPairs <- iPairStart:iPairEnd
nPairs     <- length(indexPairs)
Ne         <- 2*972 ## Number of haplotypes
pedigrees  <- rep("",2*nPairs)
SEGS       <- NULL
FSTATS     <- NULL

for(index in 1:nPairs){
  pairID <- indexPairs[index]
  cat(paste0("# Pair #",pairID,"\n"))
  simFun  <- function(chr){
    indexParents  <- sample(1:Ne ,4)
    indexParent1a <- indexParents[1]
    indexParent1b <- indexParents[2]
    indexParent2a <- indexParents[3]
    indexParent2b <- indexParents[4]
    simulateGenotype(indexParent1a,indexParent1b,
                     indexParent2a,indexParent2b,
                     Chromosome=chr)
  }

  ResultsAllChroms <- lapply(1:22,simFun)
  
  ## Output results
  iid <- ifelse(pairID<10,paste0("00000",pairID),
                ifelse(pairID<100,paste0("0000",pairID),
                       ifelse(pairID<1000,paste0("000",pairID),
                              ifelse(pairID<10000,paste0("00",pairID),
                                     ifelse(pairID<100000,paste0("0",pairID),
                                            pairID)))))
  FID <- paste0("FAMX",iid)
  PID <- paste0("PIDX",iid)
  MID <- paste0("MIDX",iid)
  ID1 <- paste0("SIB1",iid)
  ID2 <- paste0("SIB2",iid)
  
  segs   <- NULL
  fstats <- NULL
  GEN    <- NULL
  
  for(chr in 1:22){
    ## write plink file
    gcl   <- grps[[chr]]
   
    g1    <-sapply(1:nsnps[chr],function(k) gcl[k,1+ResultsAllChroms[[chr]]$genotype1[k]])
    g2    <-sapply(1:nsnps[chr],function(k) gcl[k,1+ResultsAllChroms[[chr]]$genotype2[k]])
    
    ind   <- grps_pr_index[[chr]]
    gen   <- rbind(g1[ind],g2[ind])
    
    GEN   <- cbind(GEN,gen)
    
    ## compile segments
    segs    <- rbind(segs,ResultsAllChroms[[chr]]$segments)
    fstats  <- rbind(fstats,ResultsAllChroms[[chr]]$fstats)
  }
  SEGS   <- rbind(SEGS,cbind.data.frame(FID=FID,IID1=ID1,IID2=ID2,segs))
  FSTATS <- rbind(FSTATS,cbind.data.frame(FID=FID,IID1=ID1,IID2=ID2,fstats))
  
  ped <- cbind(FID=c(FID,FID),
               IID=c(ID1,ID2),
               PID=c(PID,PID),
               MID=c(MID,MID),
               SEX=sample(1:2,2,replace=TRUE),
               PHENO=c(0,0),
               GEN)
  ped <- c(paste(ped[1,],collapse = "\t"),paste(ped[2,],collapse = "\t"))
  pedigrees[(2*index-1):(2*index)] <- ped
}
cat("# Writing output files...\n")
write(pedigrees,pedfilename)

## reorder SEGS
SEGStmp <- NULL
ufid    <- unique(SEGS[,"FID"])
for(u in ufid){
 tmpfam <- SEGS[which(SEGS[,"FID"]==u),]
 for(chr in 1:22){
   l    <- which(tmpfam[,"CHR"]==chr)
   if(length(l)>0){
     xtmp <- tmpfam[which(tmpfam[,"CHR"]==chr),]
     xtmp <- xtmp[order(xtmp[,"START"]),]
     SEGStmp <- rbind(SEGStmp,xtmp)
   }
 }
}
SEGS <- SEGStmp
save(list=c("SEGS","FSTATS"),file=RDatafilename)


