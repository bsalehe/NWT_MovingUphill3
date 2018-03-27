#Loading/saving/packages needed
#This was originally done on R 3.4.1 "Kite-eating tree"

setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill2")
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata")

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis2.Rdata")  #alternate between 1 and 2

load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis2.Rdata")


#for data cleaning

#for installing phyloseq
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

library(phyloseq)
#packageVersion("phyloseq")
library(picante) #for phylogenetic diversity

#for cooccurrence networks
#library(foreach)
#library(doParallel)
library(HMSC)
library(vegan)
library(corrplot)
library(circlize)
library(Hmisc)
library(boral)
library(Matrix)

#for plotting
library(igraph)
#library(fdrtool)
library(ggplot2)
library(grid) #for unit function in ggplot2 for legend 

library(vegan)

#for network stats
library(NetIndices)

#for manipulating datasets for plotting 
library(tidyr)
library(dplyr)
library(plotrix)

detach(package:igraph)
sessionInfo()

#extra not needed
library(reshape)
library(plotrix)
library(Kendall)


library(data.table)
#library(BiodiversityR) #this requires X11 and takes a while to load, you need to close the window that it opens in rcommander

