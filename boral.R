## Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM with 4 latent variables, a random site effect, and the effect of fixed environmental variables

## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

library(boral) #Need version 0.7 or later, available on CRAN.
library(Matrix)


####### Get species and environment data together #####

#microbes relative abundance data, filtered doubletons, singletons, taxa with < .002 summed rel abundance
comm.dataEukS<-datEukS3otu3
comm.dataEukN<-datEukN3otu3
comm.dataBac<-datBacS3otu3
comm.dataITS<-datITSS3otu3

#microbes count data, filtered doubletons, singletons, taxa with < .002 summed rel abundance
#use this
comm.dataEukS<-datEukS3cotu3
comm.dataEukN<-datEukN3cotu3
comm.dataBac<-datBacS3cotu3
comm.dataITS<-datITSS3cotu3

#plants
plantcomp2


#Merge things. all microbe datasets (not plants) should have the same 90 samples
#first merge comm.dataEuk with comm.data16S
#I need to remove all the description columns in one of the files, then merge
comm.dataEukSa<-cbind(Sample_name=comm.dataEukS$Sample_name,comm.dataEukS[,-c(1:31)])
comm.dataALL1<-merge(comm.dataEukN,comm.dataEukSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataBaca<-cbind(Sample_name=comm.dataBac$Sample_name,comm.dataBac[,-c(1:31)])
comm.dataALL2<-merge(comm.dataALL1,comm.dataBaca,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataITSa<-cbind(Sample_name=comm.dataITS$Sample_name,comm.dataITS[,-c(1:31)])
comm.dataALL3<-merge(comm.dataALL2,comm.dataITSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataALL3$Sample_name
dim(comm.dataALL3)[2]-31
1091+142+4853+1122 #matches, good, 7208 microbial taxa total

#then merge plants with microbes
comm.dataALL4<-merge(comm.dataALL3,plantcomp2,"Sample_name",sort=F,all.y=F)
comm.dataALL4$Sample_name

#substitute S for the N, since the mapping file was from the nematode dataset
comm.dataALL4$X.SampleID<-sub("N", "S", comm.dataALL4$X.SampleID) 

#delete mapping file data except for X.SampleID
comm.dataALL5<-comm.dataALL4[,-c(1,3:31)]
comm.dataALL5[1:10,1:10]



# biogeochemistry and plant density/cover data
biogeo6$X.SampleID

#Merge the biogeo6 with comm.dataALL, then split them to make sure the same samples and order are in each dataset
comm.bio<-merge(biogeo6,comm.dataALL5)
comm.bio[1:10,1:60]

#the comm.bio.csv file that is saved here is from the old bioinformatics. I saved it b/c the R environment file takes so long to load. however I might want to reinstate this, if this environment gets really big
#write.csv(comm.bio,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.bio.csv",row.names=F)
#comm.bio<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/comm.bio.csv")





##### Split and subset datasets for modeling #####

dim(comm.bio)

hmscY<-comm.bio[,54:7315] #for count data
#hmscY<-comm.bio[,54:5895] #for relabundance data
#hmscY<-comm.bio[,54:300] #for practice

rownames(hmscY)<-comm.bio$X.SampleID
hmscY[1:10,1:10]

#the intercept is a holdover from using the HMSC package, I will leave it in for how incase I ever want to go back to that
hmscX<-data.frame(inter=rep(1,75),snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture,lomehi=comm.bio$lomehi) #,plantcov=comm.bio$plantcov  ,whc=comm.bio$WHC
rownames(hmscX)<-comm.bio$X.SampleID

rcorr(as.matrix(hmscX[,2:(dim(hmscX)[2]-1)]))
#plot(hmscX$TC,hmscX$moisture)

#select lo/me/hi
ind<-which(hmscX$lomehi=="hi")
hmscXb<-hmscX[ind,]
hmscYb<-hmscY[ind,]

#select species with greater than X (X+1 or more) occurrences and remove lo me hi (since you can't have text in a matrix or the whole matrix becomes character)
ind<-which(colSums(hmscYb>0)>8)
length(ind)
hmscYc<-hmscYb[,ind]
hmscXc<-hmscXb[,1:dim(hmscXb)[2]-1]#
dim(hmscYc)
dim(hmscXc)

#the y data are not normal (the only options I have are normal, binary, poisson, overdispersed poisson), so I could do a sqrt transformation on Y (log(0) is -Inf). log(x+1) doesn't work since the proportions are so low, could do log(x*100+1) but the sqrt actually makes it more normal
#hmscYd<-sqrt(hmscYc*100)
#hmscYd<-hmscYc*100 #this is odd, 9 taxa have negative R2
#hmscYd<-log(hmscYc+1)
#hist(hmscYc[,12])
#hist(hmscYd[,12])

#check if the values are too low that some tolerance is messing up the CI estimates, yes important to scale y, instead of scaling Y I will use the sqrt transform of the percent. I will scale x, since they differ so much in range
hmscXd<-scale(hmscXc)
hmscXd[,1]<-1

#hmscYd2<-scale(hmscYd)

#make them matrices
hmscXe<-as.matrix(hmscXd)

hmscYe<-as.matrix(hmscYc)
#hmscYe<-as.matrix(hmscYd2)

dim(hmscYe)
dim(hmscXe)


## Covariates need to be stored as a matrix for boral. no need for intercept
covX <- hmscXe[,-1]




##### Fit the LVM using boral and calculate residual correlation matrix#####

#List of files produced:
#fit.hilv4occ10exp4
#fit.lolv4occ10exp4
#fit.hilv4occ9exp4
#fit.lolv4occ9exp4
#fit.hilv4occ8exp4, do not have the corresponding rescor for this b/c my computer crashed
#fit.melv4occ9exp4
#fit.melv3occ9exp4
#fit.melv5occ9exp4
#fit.melv6occ9exp4
#fit.melv2occ9exp4
#fit.lolv4occ9exp4f - f means final model fitted with long mcmc chain, start 2:00pm, model finished 4:20pm, rescor finished 
#fit.melv4occ9exp4f - f means final model fitted with long mcmc chain
#fit.hilv4occ9exp4f - f means final model fitted with long mcmc chain
#fit.lolv4occ9exp4nosite - with no site random effect, fit with long chain
#fit.lolv4occ9exp0nosite - with no site random effct and only latent vaiables
#fit.melv4occ9exp4nosite


#Using the default mcmc parameters, the models take about 2 hrs to fit.  mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
#Using shorter chains, it takes about 12 min to fit.  mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 3, seed = 123)
# in the tutorial they add this to the model fitting code, hypparams = c(20,20,20,20), however, now if you wanted to change this you need to put it in a prior.control statement or something. I am just using the default here
fit.hilv4occ9exp4f <- boral(y = hmscYe, X = covX, num.lv = 4, family = "negative.binomial", row.eff = "random", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#calc.ics = F, use ics (information criteria at your own risk) 
#fitting models for comparing percent covariance explained by env, take row.eff out
fit.hilv4occ9exp4nosite <- boral(y = hmscYe, X = covX, num.lv = 4, family = "negative.binomial", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))
fit.lolv4occ9exp0nosite <- boral(y = hmscYe, X = NULL, num.lv = 4, family = "negative.binomial", save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))


#Calculate residual correlation matrix
#Using
#Using shorter chains (see above comment), takes about 11 min to fit
#This will not calculate with a cuttoff of occ6, R studio crashes on my laptop
rescor.lolv4occ9exp4f <- get.residual.cor(fit.lolv4occ9exp4f) 
rescor.melv4occ9exp4f <- get.residual.cor(fit.melv4occ9exp4f) 
rescor.hilv4occ9exp4f <- get.residual.cor(fit.hilv4occ9exp4f) 

rescor.lolv4occ9exp4nosite <- get.residual.cor(fit.lolv4occ9exp4nosite) 
rescor.lolv4occ9exp0nosite <- get.residual.cor(fit.lolv4occ9exp0nosite) 
rescor.melv4occ9exp4nosite <- get.residual.cor(fit.melv4occ9exp4nosite) 
rescor.hilv4occ9exp4nosite <- get.residual.cor(fit.hilv4occ9exp4nosite) 


#Use modified function below to get 90% CI
#res.corshiocc8.90<-get.residual.cor2(fit.lvmhiocc8)
save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis2.Rdata")  


##### Look at results and check convergence/fit #####

#Model fit
fit.melv2occ9exp4$ics[1]
fit.melv3occ9exp4$ics[1]
fit.melv4occ9exp4$ics[1]
fit.melv5occ9exp4$ics[1]
fit.melv6occ9exp4$ics[1]
i<-1
plot(2:6,c(fit.melv2occ9exp4$ics[i],
           fit.melv3occ9exp4$ics[i],
           fit.melv4occ9exp4$ics[i],
           fit.melv5occ9exp4$ics[i],
           fit.melv6occ9exp4$ics[i]),type = "b")


summary(fit.hilv4occ9exp4f) # To look at estimated parameter values
fit.lolv4occ10exp4$hpdintervals # 95% credible intervals for model parameters.

#check information criteria
fit.lolv4occ10exp4$ics

#Check convergence
#Geweke diagnostic - a z test testing whether the first 10% and the last 50% are diffrent (i think those are the fractions, doesn't really matter exactly), if it is significant, then the means are different and it didn't converge
plot(get.mcmcsamples(fit.hilv4occ9exp4)[,1])
plot(get.mcmcsamples(fit.hilv4occ9exp4f)[,1])

#the order is effct of snowdepth for each of th 600 species, then effect of TC, then pH, then moisture 
mcmchi<-get.mcmcsamples(fit.lolv4occ9exp4f)
dim(mcmchi)
colnames(mcmchi)[2500:3000]
mcmchi[1:10,1:5]

#TRUE means these did not converge
gew.pvals <- 2*pnorm(abs(unlist(fit.hilv4occ9exp4f$geweke.diag[[1]])), lower.tail = FALSE)
length(gew.pvals)
gew.pvals[1:5]
gew.pvals[which(gew.pvals<.05)] #technically these did not converge, however, the trace plots look fine to me
p.adjust(gew.pvals, method = "holm")

fit.hilv4occ9exp4f$geweke.diag
str(fit.hilv4occ9exp4f$geweke.diag)

#example of one that did not converge
#(1st species) N6f914ead2160e51670d3dc70c25e107b for snowdepth did not converge, but looking at the trace plot, it seems fine
#geweke diagnostic
fit.hilv4occ9exp4f$geweke.diag$geweke.diag$X.coefs[1:5,]
#trace plot (it is the very first parameter)
plot(get.mcmcsamples(fit.hilv4occ9exp4f)[,1])
#mean of the mcmc chain to make sure I'm looking at the right parameter
mean(get.mcmcsamples(fit.hilv4occ9exp4f)[,1]) #mean is -1.710607
#mean of the extractd model coefficients (to mak sure I'm looking at the right parameter)
fit.hilv4occ9exp4f$X.coefs.mean  #-1.710607072, yes checks



##### Percent covariation explained by env #####
#One approach to quantify how much of the species co-occurrence is explained by covariates (that is, how well the predictor variables describe the species assemblage) is through differences in the trace of the estimated residual covariance matrix induced by the latent variables (Warton et al. 2015). From the above code, this can be obtained as rescors$trace. For the spider data set, when we compared a pure latent variable model (similar to the equation 1 but without site effects) to the correlated response model, the trace decreased from 178.92 to 107.92. This implies that environmental covariates accounted for approximately 40% of the covariation between species.
(107.92-178.92)/178.92

#Need to fit a model with only latent variables
rescor.lolv4occ9exp4f$trace
rescor.lolv4occ9exp0f$trace






##### Environmental correlations #####
envcor.lolv4occ9exp4f<-get.enviro.cor(fit.lolv4occ9exp4f)
corrplot(envcor.lolv4occ9exp4f$cor[1:100,1:100], type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45,tl.pos="n")
corrplot(envcor.lolv4occ9exp4f$sig.cor[1:100,1:100], type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45,tl.pos="n")




##Extract model coefficients
cbind(fit.lolv4occ10exp4$lv.coefs.mean,fit.lolv4occ10exp4$X.coefs.mean)
fit.hilv4occ9exp4f$X.coefs.mean

##Dunn-Smyth residual plots to check model assumption, outliers etc. The first plot should not have a funnel
plot(fit.lolv4occ9exp4f)
plot(fit.melv4occ9exp4f)
plot(fit.hilv4occ9exp4f)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lolv4occ10exp4)



##### Extract number of significant correlations #####

(length(which(rescor.melv3occ9exp4$sig.correlaton!=0))-dim(rescor.melv3occ9exp4$sig.correlaton)[1])/2
(length(which(rescor.melv4occ9exp4$sig.correlaton!=0))-dim(rescor.melv4occ9exp4$sig.correlaton)[1])/2
(length(which(rescor.melv5occ9exp4$sig.correlaton!=0))-dim(rescor.melv5occ9exp4$sig.correlaton)[1])/2

(length(which(rescor.lolv4occ9exp4f$sig.correlaton!=0))-dim(rescor.lolv4occ9exp4f$sig.correlaton)[1])/2
(length(which(rescor.melv4occ9exp4f$sig.correlaton!=0))-dim(rescor.melv4occ9exp4f$sig.correlaton)[1])/2
(length(which(rescor.hilv4occ9exp4f$sig.correlaton!=0))-dim(rescor.hilv4occ9exp4f$sig.correlaton)[1])/2

#strange - using more mcmc iterations changs how many significant interactions there are (ex: lo: 3000 iter - 3398 interactions, 30000 iter - 1806 interactions). the only way I can explain this is that thre is more error in th paramter estimates with the short chain (i.e. th mixing is wider and the density histogram is wider), thus if you are more confident in th location of the species in ordination space, thn the locations will not overlap as much, and you will say that thy are less correlated.


###### Use corrplot package to plot residual correlations between species #####

corrplot(rescor.melv3occ9exp4$sig.correlaton[1:100,1:100], diag = F, type = "lower", title = "Residual correlations from LVM", mar=c(3,0.5,2,1), tl.srt = 45,method = "color")#
corrplot(rescor.melv4occ9exp4$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")
corrplot(rescor.melv5occ9exp4$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")

corrplot(rescor.lolv4occ9exp4f$sig.correlaton[1:200,1:200], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color",tl.pos="n")




##### Plotting with igraph #####

##### colors #####

labelcols<-data.frame(rbind(c("Bacteria","#7879BC"),
                            c("Eukaryota","#94BA3C"),
                            c("Mesofauna","#ff9c34"),
                            c("Fungi","#F6EC32"),
                            c("Plant","#E95275")))
colnames(labelcols)=c("group","color")


# including photosynthetic/not information
labelcols<-data.frame(rbind(c("PhotosyntheticEukaryota","#49874c"),# 466D24
                            c("HeterotrophicEukaryota","#673482"),
                            c("PhotosyntheticBacteria","#94BA3C"),
                            c("HeterotrophicBacteria","#7879BC"),
                            c("ChemoautotrophicBacteria","#6295cd"),
                            c("UnknownEukaryota","gray50"),
                            c("UnknownBacteria","gray70"),
                            c("Mesofauna","#ff9c34"),
                            c("Fungi","#F6EC32"),
                            c("Plant","#E95275")))
colnames(labelcols)=c("group2","color")


head(labelfile)

#labelsall<-merge(labelfile,labelcols,"group",all.x=F,all.y=F) #"labels"
labelsall<-merge(labelfile,labelcols,"group2",all.x=F,all.y=F) #"labels"
labelsall$color<-as.character(labelsall$color)
head(labelsall)
labelsall$group2<-factor(labelsall$group2,levels=c("HeterotrophicBacteria","PhotosyntheticBacteria","ChemoautotrophicBacteria","UnknownBacteria","Fungi","HeterotrophicEukaryota","PhotosyntheticEukaryota","UnknownEukaryota","Mesofauna","Plant"))

#old                            
# c("NonphotosyntheticEukaryota","#673482"),
# c("PhotosyntheticEukaryota","#466D24"),
# c("Fungi","#F6EC32"),
# c("Metazoa","#ff9c34"),
# c("Plant","#E95275"),
# c("PhotosyntheticBacteria","#94BA3C"),
# c("NonphotosyntheticBacteria","#7879BC"),
# c("OtherMetazoa","black"),
# c("Nematoda","gray30"),
# c("Tardigrada","gray55"),
# c("Rotifera","gray80"),
# c("Arthropoda","white"),
# c("AF","red"),
# c("AP","red"),
# c("BF","black"),
# c("FF","gray30"),
# c("OM","gray55"),
# c("PP","gray80"),
# c("RA","white"),
# c("unknown","red")))

#colnames(labelcols)=c("labels","color")

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/legend2.pdf")
plot(c(1,1),c(1,1))
#legend("topright",c("Bacteria","Small Eukaryota","Mesofauna","Fungi","Plants"),pt.bg=c("#7879BC","#94BA3C","#ff9c34","#F6EC32","#E95275"),bty="n",pch=21,cex=1.4)
legend("topright",c("Heterotrophic bacteria","Photosynthetic bacteria","Chemoautotophic bacteria","Unknown bacteria","Heterotrophic Eukaryota","Photosynthetic Eukaryota","UnknownEukaryota","Mesofauna","Fungi","Plants"),pt.bg=c("#7879BC","#94BA3C","#6295cd","gray70","#673482","#466D24","gray50","#ff9c34","#F6EC32","#E95275"),bty="n",pch=21,cex=1.4)
legend("topleft",c("Positive","Negative"),col=c("#ce4d42","#687dcb"),lty=1,bty="n",cex=1.4)
#legend("top",as.character(1:10),col=c("#111110","#660011","#A80013","#118877","#4c3d3e","#118877","#7f783f","#aa8888","#aabbdd","#ff99a4"),lty=1,lwd=3)
dev.off()
 
#colors from nico: 111110,660011,112288,A80013,4c3d3e,118877,7f783f,aa8888,aabbdd,ff99a4, Ffccd1,ddd7d7,d8d3ad,e5001a



##### lo #####
#creating sparse matrix
colMatlo<-rescor.lolv4occ9exp4f$sig.correlaton
colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton>0)]<-1
#colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton>0)]<-0
colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton<0)]<- -1
#colMatlo[which(rescor.lolv4occ9exp4f$sig.correlaton<0)]<- 0

colMatlo<-rescor.lolv4occ9exp4nosite$sig.correlaton
colMatlo[which(rescor.lolv4occ9exp4nosite$sig.correlaton>0)]<-1
colMatlo[which(rescor.lolv4occ9exp4nosite$sig.correlaton<0)]<- -1

#colMatlo<-rescor.lolv4occ9exp4f$sig.correlaton
#colMatlo[which(colMatlo>.6)]<-1
#colMatlo[which(colMatlo<(-.6))]<- -1
#colMatlo[which(colMatlo<.6&colMatlo>(-.6))]<-0

temp<-colMatlo[upper.tri(colMatlo)]
temp2<-temp[temp!=0]
hist(temp2)
sort(temp2)
length(temp2)

graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges

length(which(myedgelistlo$weight==1))
length(which(myedgelistlo$weight==-1))
length(which(myedgelistlo$weight==1))/(length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1)))

graphlo2<-graph.edgelist(as.matrix(myedgelistlo[,1:2]),directed=FALSE)
graphlo2

verticesgraphlo<-data.frame(otu=rownames(as.matrix(V(graphlo2))))
colorgraphlo<-merge(verticesgraphlo,labelsall,"otu",all.y=F,all.x=F,sort=F)
#shapesgraplo<-ifelse(colorgraphlo$group%in%c("Eukaryota"),"csquare",'circle')

##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderlo<-order(colorgraphlo$group2)
#orderlo<-order(verticesgraphlo$otu)
graphlo2$layout <- layout_in_circle(graphlo2,order=orderlo)
#graphlo2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networklocircle.pdf") 
plot(graphlo2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistlo$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphlo$color,edge.width=.7)#,vertex.shape=shapesgraplo
#dev.off()

colorgraphlo[which(colorgraphlo$group=="Mesofauna"),]
colorgraphlo[which(colorgraphlo$group2=="PhotosyntheticBacteria"),]

myedgelistlo[which(myedgelistlo$X1=="B783ef4ce2388b995de6b9b27b0c9209e"|myedgelistlo$X2=="B783ef4ce2388b995de6b9b27b0c9209e"),]

temp<-colorgraphlo[which(colorgraphlo$group=="Mesofauna"),"otu"]
temp2<-myedgelistlo[which(myedgelistlo$X1%in%temp|myedgelistlo$X2%in%temp),]
dim(temp2)


##### me #####
#creating sparse matrix
colMatme<-rescor.melv4occ9exp4f$sig.correlaton
colMatme[which(rescor.melv4occ9exp4f$sig.correlaton>0)]<-1
colMatme[which(rescor.melv4occ9exp4f$sig.correlaton<0)]<- -1

# colMatme<-rescor.melv4occ9exp4f$sig.correlaton
# colMatme[which(colMatme>.6)]<-1
# colMatme[which(colMatme<(-.6))]<- -1
# colMatme[which(colMatme<.6&colMatme>(-.6))]<-0

colMatme<-rescor.melv4occ9exp4nosite$sig.correlaton
colMatme[which(rescor.melv4occ9exp4nosite$sig.correlaton>0)]<-1
colMatme[which(rescor.melv4occ9exp4nosite$sig.correlaton<0)]<- -1

graphme1<-graph_from_adjacency_matrix(colMatme, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistme<-data.frame(as_edgelist(graphme1),weight=E(graphme1)$weight) #just the edges

length(which(myedgelistme$weight==1))
length(which(myedgelistme$weight==-1))
length(which(myedgelistme$weight==1))/(length(which(myedgelistme$weight==1))+length(which(myedgelistme$weight==-1)))

graphme2<-graph.edgelist(as.matrix(myedgelistme[,1:2]),directed=FALSE)
graphme2

verticesgraphme<-data.frame(otu=rownames(as.matrix(V(graphme2))))
colorgraphme<-merge(verticesgraphme,labelsall,"otu",all.y=F,all.x=F,sort=F)

##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderme<-order(colorgraphme$group2)
#orderme<-order(verticesgraphme$otu)
graphme2$layout <- layout_in_circle(graphme2,order=orderme)
#graphme2$layout <- layout_in_circle

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkmecircle.pdf") 
plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
dev.off()

temp<-colorgraphme[which(colorgraphme$group=="Mesofauna"),"otu"]
temp2<-myedgelistme[which(myedgelistme$X1%in%temp|myedgelistme$X2%in%temp),]
dim(temp2)



##### hi #####
#creating sparse matrix
colMathi<-rescor.hilv4occ9exp4f$sig.correlaton
colMathi[which(rescor.hilv4occ9exp4f$sig.correlaton>0)]<-1
colMathi[which(rescor.hilv4occ9exp4f$sig.correlaton<0)]<- -1
#colMathi[which(rescor.hilv4occ9exp4f$sig.correlaton<0)]<- 0

# colMathi<-rescor.hilv4occ9exp4f$sig.correlaton
# colMathi[which(colMathi>.5)]<-1
# colMathi[which(colMathi<(-.5))]<- -1
# colMathi[which(colMathi<.5&colMathi>(-.5))]<-0

colMathi<-rescor.hilv4occ9exp4nosite$sig.correlaton
colMathi[which(rescor.hilv4occ9exp4nosite$sig.correlaton>0)]<-1
colMathi[which(rescor.hilv4occ9exp4nosite$sig.correlaton<0)]<- -1

graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges

length(which(myedgelisthi$weight==1))
length(which(myedgelisthi$weight==-1))
length(which(myedgelisthi$weight==1))/(length(which(myedgelisthi$weight==1))+length(which(myedgelisthi$weight==-1)))

graphhi2<-graph.edgelist(as.matrix(myedgelisthi[,1:2]),directed=FALSE)
graphhi2

verticesgraphhi<-data.frame(otu=rownames(as.matrix(V(graphhi2))))
colorgraphhi<-merge(verticesgraphhi,labelsall,"otu",all.y=F,all.x=F,sort=F)

#order starts at 3:00 and goes counterclockwise
##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderhi<-order(colorgraphhi$group2)
#orderhi<-order(verticesgraphhi$otu)
graphhi2$layout <- layout_in_circle(graphhi2,order=orderhi)
#graphhi2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkhicircle.pdf") 
plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
#dev.off()

colorgraphhi[which(colorgraphhi$group=="Mesofauna"),]

temp<-colorgraphhi[which(colorgraphhi$group=="Mesofauna"),"otu"]
temp2<-myedgelisthi[which(myedgelisthi$X1%in%temp|myedgelisthi$X2%in%temp),]
dim(temp2)

#Creating a subgraph
colorgraphhi2<-colorgraphhi[which(colorgraphhi$group2=="Mesofauna"),]
myedgelisthi2<-myedgelisthi[which(myedgelisthi[,"X1"]%in%colorgraphhi2$otu|myedgelisthi[,"X2"]%in%colorgraphhi2$otu),]
graph3<-subgraph.edges(graphhi2, eids=which(myedgelisthi2[,"X1"]%in%colorgraphhi2$otu|myedgelisthi2[,"X2"]%in%colorgraphhi2$otu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)



###### network statistics #######
graph.density(graphlo2)
graph.density(graphme2)
graph.density(graphhi2)

length(E(graphlo2))/length(V(graphlo2))
length(E(graphme2))/length(V(graphme2))
length(E(graphhi2))/length(V(graphhi2))

temp<-colorgraphlo; temp$ones<-1
aggregate.data.frame(temp$ones,by=list(temp$group2),sum)

temp<-colorgraphme; temp$ones<-1
aggregate.data.frame(temp$ones,by=list(temp$group2),sum)

temp<-colorgraphhi; temp$ones<-1
aggregate.data.frame(temp$ones,by=list(temp$group2),sum)












colMat<-res.corslolv4occ8exp5$sig.correlaton
colMat[which(res.corslolv4occ8exp5$sig.correlaton>0)]<-1
colMat[which(res.corslolv4occ8exp5$sig.correlaton<0)]<- -1

colMat<-res.corshilv4occ8exp4$sig.correlaton
colMat[which(res.corshilv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corshilv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corsmelv4occ8exp4$sig.correlaton
colMat[which(res.corsmelv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corsmelv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corslolv4occ8exp4$sig.correlaton
colMat[which(res.corslolv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corslolv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corshilv4occ8exp4$sig.correlaton
colMat[which(colMat>.6)]<-1
colMat[which(colMat<(-.6))]<- -1
colMat[which(colMat<.6&colMat>(-.6))]<-0

colMat<-res.corslolv4occ8exp4$sig.correlaton
colMat[which(colMat>.6)]<-1
colMat[which(colMat<(-.6))]<- -1
colMat[which(colMat<.6&colMat>(-.6))]<-0

#hist(colMat)
#dim(colMat)

#colMat2<-as(colMat, "dgCMatrix")
graph1<-graph_from_adjacency_matrix(colMat, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelist<-data.frame(as_edgelist(graph1),weight=E(graph1)$weight) #just the edges

length(which(myedgelist$weight==1))
length(which(myedgelist$weight==-1))
length(which(myedgelist$weight==1))/(length(which(myedgelist$weight==1))+length(which(myedgelist$weight==-1)))

#graph2<-simplify(graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE))
graph2<-graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE)
graph2

#(sum(colMat,na.rm=T)-322)/2 #(134 interactions)
#(24*24-24)/2 #every pairwise interaction (276)

graph2$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph2))))
colorgraph<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraph$color)#,layout=l3

unique(tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%colorgraph$oldotu])
temp<-myedgelist[which(myedgelist$X1%in%c("8f0fb68673a8b3c26f3fd210e29b7916")|myedgelist$X2%in%c("8f0fb68673a8b3c26f3fd210e29b7916")),]
8f0fb68673a8b3c26f3fd210e29b7916 #clostridium bowmanii, byproducts of acetic and butyric acid
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

temp<-myedgelist[which(myedgelist$X1%in%c("6472eb8b1e09f892aca2f23182962903")|myedgelist$X2%in%c("6472eb8b1e09f892aca2f23182962903")),]
8f0fb68673a8b3c26f3fd210e29b7916  #bradyrhizobium fixes N
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%c(as.character(temp2$X1),as.character(temp2$X2))] #one chlorophyceae - possible syntrophic relationship
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]


"631c2b290edb6b71d9e6b5577a0613ae"# syntrophobacteriaceae
"8a8895d160b2fda424ad5446e9bf8f4e"#Desulfosporosinus sulfate reducing
temp<-myedgelist[which(myedgelist$X1%in%c("8a8895d160b2fda424ad5446e9bf8f4e")|myedgelist$X2%in%c("8a8895d160b2fda424ad5446e9bf8f4e")),]
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

"b8f64099026ffffde7cd2c6aa4642240" #nitrosomonadaceae, ammonia oxidizers
temp<-myedgelist[which(myedgelist$X1%in%c("b8f64099026ffffde7cd2c6aa4642240")|myedgelist$X2%in%c("b8f64099026ffffde7cd2c6aa4642240")),]
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

myedgelist[myedgelist$X1%in%temp,]
myedgelist[myedgelist$X2%in%temp,]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(temp)]
colorgraph[which(colorgraph$oldotu=="854124b2d8e0badeb138eb9e0e1808c2"),]

#Creating a subgraph
colorgraph3<-colorgraph[which(colorgraph$labels=="Bacteria"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Bacteria"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
tax_table(datEukN)[rownames(tax_table(datEukN))%in%colorgraph3$oldotu]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph3<-subgraph.edges(graph2, eids=which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraph$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)

#investigating nematode relationships
#lo
temp<-"b85db42310af5ddb08354eef2427cc8e" #omnivore nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="b85db42310af5ddb08354eef2427cc8e"|myedgelist$X2=="b85db42310af5ddb08354eef2427cc8e"),]
temp2<-temp1[temp1$weight==-1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#5 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#2 euks
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#0
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#6
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#1
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#6
tax_table(datEukN)[rownames(tax_table(datEukN))%in%temp2$X2]#0

#hi
6f914ead2160e51670d3dc70c25e107b __Aphelenchida fungal feeder
e485c9c4bdda247bc813c36275d99dce __Plectidae bacterial feeder
temp<-"6f914ead2160e51670d3dc70c25e107b" #ff nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="6f914ead2160e51670d3dc70c25e107b"|myedgelist$X2=="6f914ead2160e51670d3dc70c25e107b"),]
#only positive rels
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#2 bacteria
temp<-"e485c9c4bdda247bc813c36275d99dce" #bf nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="e485c9c4bdda247bc813c36275d99dce"|myedgelist$X2=="e485c9c4bdda247bc813c36275d99dce"),]
temp2<-temp1[temp1$weight==-1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#3 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#0 euks
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#0
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#3 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]
# 2 plants 





#trying to keep the layout the same and plot only partial graph
l3<-layout_in_circle(graph2)#layout_with_fr(graph2) #I don't want to run this accidentally
rownames(l3) <- V(graph2)$name
l3b<-layout.norm(l3,xmin=-1,xmax=1,ymin=-1,ymax=1)
#plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph3,layout=l3b,rescale=F,xlim=c(-1,1),ylim=c(-1,1))




##Exploring nematode correlations in hi
#use colorgraph from below
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
colorgraph2
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%colorgraph2$oldotu]
tax_table(datEukN)[rownames(tax_table(datEukN))%in%colorgraph2$oldotu]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%colorgraph2$oldotu]
temp<-c("a64f0e90f25d607e4c34fffa870f19e4","71f3f9709801ee03a20e33c46a8b797d","a7a7353d6231f298e75b69f4161926a1","ee45b87e6ac04257510a6d99644760cd","eafa1d1d561a47042f099e619cd25693","e4ae90d4c05ef00793fc6b93fb6a9af7","f7b389ab203fda4cc3cf1fffe23004a6","9bd3df7f714986d788a2cce65be56d10","7a97813269020725286a63989363ceef","d3898987a65cc50f9226529374935e74","854124b2d8e0badeb138eb9e0e1808c2")#lo 10
temp<-c("f4e7445ced25e7028c593551658c091a","f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8","479660228dd139e20ad8ff42698a123e","65cb8aa63bed1f44c9a576ae2e94e0ab")#hi 8
temp<-c("f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8")#hi 10
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$oldotu%in%temp),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Fungi"),]

#tax_table(datEukN3)[rownames(tax_table(datEukN3))=="6f914ead2160e51670d3dc70c25e107b"]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph3<-graph.edgelist(as.matrix(myedgelist2[,1:2]),directed=FALSE)
graph3
graph3$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph3))))
colorgraph4<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"red","blue"),vertex.color=colorgraph$color,layout=l3b)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)



#using hi occ8 with correlation >.5 or <-.5, I get 10 mesofauna taxa
#using hi occ10 significnat correlations, I get 5 mesofauna and 4 heterotrophic eukaryotes. 67/51
#using hi occ8 significnat correlations, I get 5 mesofauna and 7 heterotrophic eukaryotes. 53/45
#lo 10, significant, 1 mesofuana, 10 heterotrophic euks 128/81
#using hi occ8 with 90%CI, 7 mesofauna (didn't look for euks)









##### Chord diagrams #####
(length(which(res.corslo$sig.correlaton!=0))-dim(res.corslo$sig.correlaton)[1])/2
colMat <- matrix(NA, nrow = nrow(rescor.lolv4occ10exp4$sig.correlaton), ncol = ncol(rescor.lolv4occ10exp4$sig.correlaton))
colMat[which(res.cors$sig.correlaton > 0.63, arr.ind = TRUE)] <- "red"
colMat[which(res.cors$sig.correlaton < -0.63, arr.ind = TRUE)] <- "blue"
chordDiagram(res.cors$sig.correlaton, symmetric = TRUE,annotationTrack = c("name", "grid"),grid.col = "grey",col=colMat)

colMat <- matrix(NA, nrow = nrow(res.corshiocc8$correlation), ncol = ncol(res.corshiocc8$correlation),dimnames = list(rownames(res.corshiocc8$correlation),rownames(res.corshiocc8$correlation)))
colMat[which(res.corshiocc8$correlation > 0.5, arr.ind = TRUE)] <- 1
colMat[which(res.corshiocc8$correlation < -0.5, arr.ind = TRUE)] <- -1
colMat[which(res.corshiocc8$correlation<.5&res.corshiocc8$correlation>(-.5))]<-0

length(which(res.cor90$correlation>.7&res.cor90$correlation<1,arr.ind = T))
length(which(res.corshi8$correlation>.7&res.corshi8$correlation<1,arr.ind = T))
hist(res.cors$correlation)
hist(res.corshi8$correlation)





#Trial runs were stored in MovingUphill3_Workspace_Analysis3b.Rdata
#fit.lvmhiocc6 
#load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis3b.Rdata")
#fit.hilv4occ8exp5
#fit.lolv4occ8exp5
#fit.hilv4occ8exp4 (and b)
#fit.melv4occ8exp4

#res.cors (hi with 4latent)
#res.corslo (lo with 4 latent)
#res.corshi8 (hi with 8 latent) - hardly any correlations significant or strong
#res.corshiocc8 (hi with 8 occurrences)
#res.corshiocc8.90 (hi with 8 occurrences and 90CI)
#res.corshiocc6 (hi 6 occ)
#res.corshilv4occ8exp5 (hi with 5 environmetnal variables)
#res.corslolv4occ8exp5



##### Checking the models with multiple regressions ####
#there is no negative binomial option in glm, so using quasipoisson
#overall, I think I'm ok with th number of nonzero points (9 or greater) and having 4 explanatory variables
hmscXe
hmscYe

#n=10
sum(hmscYe[,17]>0)
m1<-glm(hmscYe[,17]~hmscXe[,2:5],family=quasipoisson)
m1<-glm(hmscYe[,17]~hmscXe[,c(2,3,4)],family=quasipoisson)
m1<-glm(hmscYe[,17]~hmscXe[,c(2,4,5)],family=quasipoisson)
summary(m1)
m1fitted<-predict(m1,type="response")
plot(hmscXe[,3],hmscYe[,17])
points(hmscXe[,3],m1fitted,col=2)
#

###### Checking whether I should include a random plot effect #####




##### funcions #####
#90% CI
get.residual.cor2<-function (object, est = "median", prob = .9) 
{
  if (is.null(object$jags.model)) 
    stop("MCMC samples not found")
  fit.mcmc <- get.mcmcsamples(object)
  y <- object$y
  X <- object$X
  num.lv <- object$num.lv
  if (length(grep("lvs", colnames(fit.mcmc))) == 0) 
    stop("Cannot find MCMC samples corresponding to latent variables")
  n <- nrow(y)
  p <- ncol(y)
  sig_rescor_mat <- rescor_mat <- rescov_mat <- matrix(0, nrow = p, 
                                                       ncol = p)
  sig_respres_mat <- respres_mat <- matrix(0, nrow = p, ncol = p)
  if (is.null(colnames(y))) 
    colnames(y) <- 1:ncol(y)
  rownames(rescor_mat) <- colnames(rescor_mat) <- rownames(sig_rescor_mat) <- colnames(sig_rescor_mat) <- colnames(y)
  rownames(rescov_mat) <- colnames(rescov_mat) <- colnames(y)
  rownames(respres_mat) <- colnames(respres_mat) <- rownames(sig_respres_mat) <- colnames(sig_respres_mat) <- colnames(y)
  all_rescor_mat <- all.rescov_mat <- all.respres_mat <- array(0, 
                                                               dim = c(nrow(fit.mcmc), p, p))
  all_trace_rescor <- numeric(nrow(fit.mcmc))
  for (k0 in 1:nrow(fit.mcmc)) {
    lv.coefs <- matrix(fit.mcmc[k0, grep("lv.coefs", colnames(fit.mcmc))], 
                       nrow = p)
    lambdalambdaT <- tcrossprod(as.matrix(lv.coefs[, 2:(num.lv + 
                                                          1)]))
    all.rescov_mat[k0, , ] <- (lambdalambdaT)
    all_trace_rescor[k0] <- sum(diag(lambdalambdaT))
    all_rescor_mat[k0, , ] <- cov2cor(all.rescov_mat[k0, 
                                                     , ])
    all.respres_mat[k0, , ] <- ginv(all_rescor_mat[k0, , 
                                                   ])
  }
  for (j in 1:p) {
    for (j2 in 1:p) {
      if (est == "median") {
        rescov_mat[j, j2] <- median(all.rescov_mat[, 
                                                   j, j2])
        rescor_mat[j, j2] <- median(all_rescor_mat[, 
                                                   j, j2])
        respres_mat[j, j2] <- median(all.respres_mat[, 
                                                     j, j2])
      }
      if (est == "mean") {
        rescov_mat[j, j2] <- mean(all.rescov_mat[, j, 
                                                 j2])
        rescor_mat[j, j2] <- mean(all_rescor_mat[, j, 
                                                 j2])
        respres_mat[j, j2] <- mean(all.respres_mat[, 
                                                   j, j2])
      }
      sig_rescor_mat[j, j2] <- rescor_mat[j, j2]
      get.hpd.cors <- HPDinterval(as.mcmc(all_rescor_mat[, 
                                                         j, j2]), prob = .9)
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
        sig_rescor_mat[j, j2] <- 0
      sig_respres_mat[j, j2] <- respres_mat[j, j2]
      get.hpd.cors <- HPDinterval(as.mcmc(all.respres_mat[, 
                                                          j, j2]), prob = .9)
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
        sig_respres_mat[j, j2] <- 0
    }
  }
  if (est == "median") 
    final_trace <- median(all_trace_rescor)
  if (est == "mean") 
    final_trace <- mean(all_trace_rescor)
  return(list(correlation = rescor_mat, sig.correlaton = sig_rescor_mat, 
              covariance = rescov_mat, precision = respres_mat, sig.precision = sig_respres_mat, 
              trace = final_trace))
}
