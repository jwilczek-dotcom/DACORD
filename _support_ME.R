# SET PATH TO DACORD:
WD <- c("C:/DACORD")
setwd(WD)
source("functions.R")

# SET PATH TO DACORD:
WD <- c("C:/DACORD")

# Mesh error external calculation
# see 'mesh.error' for more information

# Libraries needed for used packages:
library(parallel)   # makeCluster
library(foreach)
library(iterators); library(doParallel)
library(Rvcg)

set.seed(12345)
a <- Sys.time()
mesh <- vcgImport(paste(WD,"/temp/temp.ply",sep=""), clean=F, readcolor=T)
profil <- read.table(paste(WD,"/temp/temp_profil.txt",sep=""), sep=";")
ME <- mesh.error(mesh,method="profile",list(profil=profil))
print(paste("The calculation taken:",Sys.time()-a,"seconds"))
write.table(ME, paste(WD,"/temp/temp_ME.txt",sep=""), sep=";", row.names=F, col.names=F)


