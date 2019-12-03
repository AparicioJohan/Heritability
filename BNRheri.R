
######################################################################
# TEST HERITABILITY --------------------------------------------------
######################################################################


library(BGLR)
library(tidyverse)
library(sommer)

# geno =  string with  file direction of marker information .in  
# samp =  string with  file direction of samples  .txt
# phen =  string with  file direction of phenotypic Inform `.txt`


BNHerita <- function(geno, samp, phen , genoname="Line"){

samp = read.delim(samp, header = F)[,1]
G = read.table(geno, row.names = as.character(samp), header = F)
G <- G[ order(row.names(G)) , ]

phen <- read.delim(phen)
phen <- arrange(phen, get(genoname))

Results <- matrix(nrow = 3,ncol = ncol(phen)-1, data = NA , 
                  dimnames = list(c("Narrow_Sense","Broad","RKHS"),names(phen)[-1]))

for (i in names(phen)[-1]) {
  
  
  if(i==names(phen)[2]){
    cat("\n[]==============================================================[]")
    cat("\n[]======== Marker based heritability calculation ===============[]")
    cat("\n[]=============== sommer & BGLR package  =======================[]")
    cat("\n[]======= Last update: 2019-10-25  Johan Aparicio ==============[]")
    cat("\n[]==============================================================[]\n")
    cat("\n Trait","\t", "h2_N","\t h2_B","\t h2_RKHS" ,"\tn_Phen","\tn_Gen","\tShare", "\tfinishedAt")
  }
  
  LinesA <- as.character(phen[!is.na(phen[,i]),genoname])  # A 
  LinesB <- as.character(rownames(G))                      # B
  
  
  
  GT <- G[rownames(G)%in%intersect(LinesA,LinesB),]   # Marker
  
  ### look at the 
  A <-A.mat(GT)  # additive relationship matrix
  D <-D.mat(GT)  # dominance relationship matrix
  E <-E.mat(GT)  # epistatic relationship matrix
  
  phen %>% subset(.[,1]%in%rownames(GT)) %>% droplevels()  -> phen2
  
  phen2$idd <-phen2[,genoname]
  phen2$id <-phen2[,genoname]
  
  phen2$var <- phen2[,i]
  
  ans.ADE <-mmer(var~1,
                 random=~vs(id,Gu=A)+ vs(idd,Gu=D),
                 rcov=~units,
                 data=phen2, verbose = F)
  
  hn <- pin(ans.ADE, h2~(V1)/( V1+V3) )$Estimate
  hb <- pin(ans.ADE, h2~(V1+V2)/( V1+V2+V3) )$Estimate
  
  #------ RKHS ---------
  
  y <- phen2$var
  D = as.matrix( dist( GT, method="euclidean")) ^2
  D = D / mean(D)
  h = 0.5
  K = exp(-h * D)
  # K <- A.mat(GT)
  ETA = list(list(model="RKHS",   K=K))
  
  fm = BGLR( y=y, ETA=ETA, 
              nIter=10000, burnIn=1000, thin=5, 
              verbose=FALSE)
  
  hrk <- fm$ETA[[1]]$varU/(fm$ETA[[1]]$varU + fm$varE)
  
  corrk <- cor(fm$y,fm$yHat)
  
  #----------------------
  
  
  cat("\n" , i ,
      "\t", round(hn,2),
      "\t", round(hb,2),
      "\t", round(hrk,2),
      "\t", sum(!is.na(phen[,i])),
      "\t", length(rownames(G)),
      "\t", sum(rownames(G)%in%LinesA),
      "\t", paste0(Sys.time()))
  
  Results[,i] <- c(hn,hb,hrk)
  
}
# Heritability
Herita <- data.frame(trait=rownames(t(Results)), round(t(Results),2),row.names = NULL)
invisible(Herita)
}

# Example
# BNHerita(geno = geno,samp = samp,phen = phen,genoname = "Line")

