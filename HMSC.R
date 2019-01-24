
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis5hmsc.Rdata")

library(HMSC)
library(vegan)
library(corrplot)
library(circlize)
library(Hmisc)
library(igraph)
library(plotrix)
library(tidyr)
library(dplyr)
library(ggplot2)

####### Get species and environment data together #####
#taken from boral.R and modified - filtered taxa more stringently

#microbes relative abundance data, filtered doubletons, singletons, taxa with < .002 summed rel abundance

dim(comm.dataEukS)
comm.dataEukS[1:10,1:10]

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

max(datEukS3otu3[,32:dim(datEukS3otu3)[2]])
[1] 0.4006889
max(datEukN3otu3[,32:dim(datEukN3otu3)[2]])
[1] 1
max(datBacS3otu3[,32:dim(datBacS3otu3)[2]])
[1] 0.1005081
max(datITSS3otu3[,32:dim(datITSS3otu3)[2]])
[1] 0.3743274

sort(matrix(as.matrix(datITSS3otu3[,32:dim(datITSS3otu3)[2]]),ncol=1),decreasing=T)

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
#hmscY<-comm.bio[,54:5895] #for relabundance data #I think this was from the old bioinformatics (clustering at 97% (?), otherwise the dimensions should not be different for count or rel abun data)
#hmscY<-comm.bio[,54:300] #for practice

rownames(hmscY)<-comm.bio$X.SampleID
hmscY[1:10,1:10]

#take out plants if you want
hmscY<-hmscY[,1:7208]

#the intercept is a holdover from using the HMSC package, I will leave it in for now in case I ever want to go back to that
#hmscX<-data.frame(inter=rep(1,75),snowdepth=comm.bio$snowdepth,TC=comm.bio$TC,pH=comm.bio$pH,moisture=comm.bio$moisture,lomehi=comm.bio$lomehi) #,plantcov=comm.bio$plantcov  ,whc=comm.bio$WHC
hmscX<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,TC=comm.bio$TC,plantcov=comm.bio$plantcov,lomehi=comm.bio$lomehi) #
pclo$vectors[,1:2]
rownames(hmscX)<-comm.bio$X.SampleID

#spatial xy coordinates UTMS
autoxy<-data.frame(plot=comm.bio$X.SampleID,x=comm.bio$X,y=comm.bio$Y) 
rownames(autoxy)<-comm.bio$X.SampleID

#correlations among explanatory variables
#rcorr(as.matrix(hmscX[,2:(dim(hmscX)[2]-1)]))
#plot(hmscX$plantcov,hmscX$TC)

#select lo/me/hi
ind<-which(hmscX$lomehi=="hi")
hmscXb<-hmscX[ind,]
hmscYb<-hmscY[ind,]
autoxyb<-autoxy[ind,]
#autoxyblo<-autoxyb
#autoxybme<-autoxyb
#autoxybhi<-autoxyb
autoxyb$plot<-as.factor(1:25) #if this is not 1:25 (i.e. numbered rather than a plot name), then as.HMSCdata fails

#select species with greater than X (X+1 or more) occurrences and remove lo me hi (since you can't have text in a matrix or the whole matrix becomes character)
ind<-which(colSums(hmscYb>0)>16)#15
length(ind)
hmscYc<-hmscYb[,ind]
hmscXc<-hmscXb[,1:dim(hmscXb)[2]-1]#
dim(hmscYc)
dim(hmscXc)

hmscYc[1:10,1:10]


#the y data are not normal (the only options I have are normal, binary, poisson, overdispersed poisson), so I could do a sqrt transformation on Y (log(0) is -Inf). log(x+1) doesn't work since the proportions are so low, could do log(x*100+1) but the sqrt actually makes it more normal. I will do log(x+1) on read count data
#hmscYd<-sqrt(hmscYc*100)
#hmscYd<-hmscYc*100 #this is odd, 9 taxa have negative R2
hmscYd<-log(hmscYc+1)
#hist(hmscYc[,30])
#hist(hmscYd[,30])

hmscYd[1:10,1:10]
apply(hmscYd,2,max)

#check if the values are too low that some tolerance is messing up the CI estimates, yes important to scale y, instead of scaling Y I will use the sqrt transform of the percent. I will scale x, since they differ so much in range
hmscXd<-scale(hmscXc)
hmscXd<-scale(data.frame(pchi$vectors[,1:2],hmscXc))

#hmscXd[,1]<-1

#hmscYd2<-scale(hmscYd)

#make them matrices
hmscXe<-as.matrix(hmscXd)

#hmscYe<-as.matrix(hmscYc)
hmscYe<-as.matrix(hmscYd)

#for space
#rownames(hmscYe)<-1:50
#rownames(hmscXe)<-1:50
#rownames(autoxyb)<-1:50

dim(hmscYe)
dim(hmscXe)
dim(autoxyb)

#for space
#hmscYe[hmscYe>0]<-1
#ind<-which(colSums(hmscYe)<45)
#hmscYe<-hmscYe[,ind]
#hmscYe<-hmscYe[,1:10]
#dim(hmscYe)






##### modeling #####

#make the random (residual matrix)
pimat<-data.frame(plot=1:dim(hmscYe)[1])
pimat$plot<-as.factor(pimat$plot)
rownames(pimat)<-rownames(hmscYe)

#summary from space trials: I got it to work using, autoxyb, with probit, with presence absence data, and when sample size is 50, I could only model ~10 y species. Even if I took out the x variables completely, it still crashed with 60 species
#formdata <- as.HMSCdata(Y = hmscYe, X = hmscXe, Auto=autoxyb, interceptX = F, scaleX=T)
#model <- hmsc(formdata, family = "probit", niter = 10000, nburn = 1000, thin = 10) #worked with pres/abs with 50 samples, not 25, if you try with lots of sampls cut it down to niter=100



formdata <- as.HMSCdata(Y = hmscYe, X = hmscXe, Random=pimat, interceptX = T, scaleX=T)
model <- hmsc(formdata, family = "gaussian", niter = 40000, nburn = 10000, thin = 20)
#model <- hmsc(formdata, family = "overPoisson", niter = 10000, nburn = 1000, thin = 10)

#if you don't include random, you don't get a correlation matrix
#formdata <- as.HMSCdata(Y = hmscYe, X = hmscXe, interceptX = T, scaleX=T) 
#model <- hmsc(formdata, family = "gaussian", niter = 10000, nburn = 1000, thin = 10)

modelhi<-model
modello<-model
modelme<-model

modelhi2<-model
modello2<-model
modelme2<-model

#formprior <- as.HMSCprior(formdata,family="gaussian",shrinkOverall = NULL, shrinkSpeed = NULL,shrinkLocal = NULL) #not necessary, this just generates flat priors
#formprior <- as.HMSCprior(formdata,family="overPoisson") #not necessary, this just generates flat priors
#formparam <- as.HMSCparam(formdata, formprior,latent = NULL,paramLatent = NULL, shrinkLocal = NULL, paramShrinkGlobal = NULL) #not necessary, this just generates random staring parameters


mixing <- as.mcmc(model,parameters = "paramX")
temp<-as.vector(mixing[1:900,6])
plot(temp,type="l")
#hist(temp)
mixing <- as.mcmc(model,parameters = "paramLatent")
str(mixing)

### Convert the mixing object to a matrix
#mixingDF <- as.data.frame(mixing)
#boxplot(mixingDF[,2], las = 2)

#CI
average <- apply(model$results$estimation$paramX, 1:2, mean)
### 95% confidence intervals
CI.025 <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.025)
CI.975 <- apply(model$results$estimation$paramX, 1:2, quantile, probs = 0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI), type = "n", xlab = "", ylab = "", main="paramX")
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2], code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)

#checking that intercepts make sense, yes more or less
N6f914ead2160e51670d3dc70c25e107b intercept is .47
mean(hmscYd$N6f914ead2160e51670d3dc70c25e107b) is .44
#summary(lm(hmscYg[,27]~decostand(hmscXe[,2],method="standardize"))) #yes the coefficients are very close, bdenovo90996

### Summary table
paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")


###### Variance partitioning######
#this code is essentially going through and calculating the parameters*Xdata for each group. it is stll variance though so I'm not sure how it relates to R2, it is the contribution of each parameter to the predicted value

#lo
variationPartlo <- variPart(modello,c(rep("abiotic",4),rep("biotic",2)))#
variationPartlo2 <- variPart(modello2,c(rep("space",3),rep("abiotic",3),rep("biotic",2)))#
#variationPart <- variPart(model,c(rep("abiotic",5)))
barplot(t(variationPartlo), legend.text=colnames(variationPartlo),
        args.legend=list(y=1.1, x=nrow(variationPartlo)/2, xjust=0.5, horiz=T))
#effect of abiotic:
mean(variationPartlo[,1])
#biotic (plant):
mean(variationPartlo[,2])
#plot
mean(variationPartlo[,3])

vardatlo<-data.frame(succession=rep("Early",139*2),source=rep(c("Abiotic","Plant"),each=139),taxon=c(rep("Eukaryota",12),rep("Bacteria",106),rep("Fungi",21)),var=c(variationPartlo[,1],variationPartlo[,2]))
vardatlo2<-data.frame(succession=rep("Early",139*3),source=rep(c("Space","Abiotic","Plant"),each=139),taxon=c(rep("Eukaryota",12),rep("Bacteria",106),rep("Fungi",21)),var=c(variationPartlo2[,1],variationPartlo2[,2],variationPartlo2[,3]))

vardatlo2<-data.frame(succession=rep("Early",139),taxon=c(rep("Eukaryota",12),rep("Bacteria",106),rep("Fungi",21)),Space=variationPartlo2[,1],Abiotic=variationPartlo2[,2],Plant=variationPartlo2[,3])
vardatlo2$tot<-vardatlo2$Space+vardatlo2$Abiotic+vardatlo2$Plant
vardatlo2$Spacer<-vardatlo2$Space/vardatlo2$tot
vardatlo2$Abioticr<-vardatlo2$Abiotic/vardatlo2$tot
vardatlo2$Plantr<-vardatlo2$Plant/vardatlo2$tot
vardatlo3<-vardatlo2%>%
  select(succession,taxon,Spacer,Abioticr,Plantr)%>%
  gather(source,var,Spacer:Plantr)


#me
variationPartme <- variPart(modelme,c(rep("abiotic",4),rep("biotic",2)))#
variationPartme2 <- variPart(modelme2,c(rep("space",3),rep("abiotic",3),rep("biotic",2)))#
barplot(t(variationPartme), legend.text=colnames(variationPartme),
        args.legend=list(y=1.1, x=nrow(variationPartme)/2, xjust=0.5, horiz=T))
#effect of abiotic:
mean(variationPartme[,1])
#biotic (plant):
mean(variationPartme[,2])
#plot
mean(variationPartme[,3])

vardatme<-data.frame(succession=rep("Mid",133*2),source=rep(c("Abiotic","Plant"),each=133),taxon=c(rep("Eukaryota",6),rep("Bacteria",116),rep("Fungi",11)),var=c(variationPartme[,1],variationPartme[,2]))
vardatme2<-data.frame(succession=rep("Mid",133*3),source=rep(c("Space","Abiotic","Plant"),each=133),taxon=c(rep("Eukaryota",6),rep("Bacteria",116),rep("Fungi",11)),var=c(variationPartme2[,1],variationPartme2[,2],variationPartme2[,3]))

vardatme2<-data.frame(succession=rep("Mid",133),taxon=c(rep("Eukaryota",6),rep("Bacteria",116),rep("Fungi",11)),Space=variationPartme2[,1],Abiotic=variationPartme2[,2],Plant=variationPartme2[,3])
vardatme2$tot<-vardatme2$Space+vardatme2$Abiotic+vardatme2$Plant
vardatme2$Spacer<-vardatme2$Space/vardatme2$tot
vardatme2$Abioticr<-vardatme2$Abiotic/vardatme2$tot
vardatme2$Plantr<-vardatme2$Plant/vardatme2$tot
vardatme3<-vardatme2%>%
  select(succession,taxon,Spacer,Abioticr,Plantr)%>%
  gather(source,var,Spacer:Plantr)


#hi
variationParthi <- variPart(modelhi,c(rep("abiotic",4),rep("biotic",2)))#
variationParthi2 <- variPart(modelhi2,c(rep("space",3),rep("abiotic",3),rep("biotic",2)))#
barplot(t(variationParthi), legend.text=colnames(variationParthi),
        args.legend=list(y=1.1, x=nrow(variationParthi)/2, xjust=0.5, horiz=T))
#effect of abiotic:
mean(variationParthi[,1])
#biotic (plant):
mean(variationParthi[,2])
#plot
mean(variationParthi[,3])

vardathi<-data.frame(succession=rep("Late",120*2),source=rep(c("Abiotic","Plant"),each=120),taxon=c(rep("Mesofauna",4),rep("Eukaryota",10),rep("Bacteria",101),rep("Fungi",5)),var=c(variationParthi[,1],variationParthi[,2]))
vardathi2<-data.frame(succession=rep("Late",120*3),source=rep(c("Space","Abiotic","Plant"),each=120),taxon=c(rep("Mesofauna",4),rep("Eukaryota",10),rep("Bacteria",101),rep("Fungi",5)),var=c(variationParthi2[,1],variationParthi2[,2],variationParthi2[,3]))

vardathi2<-data.frame(succession=rep("Late",120),taxon=c(rep("Mesofauna",4),rep("Eukaryota",10),rep("Bacteria",101),rep("Fungi",5)),Space=variationParthi2[,1],Abiotic=variationParthi2[,2],Plant=variationParthi2[,3])
vardathi2$tot<-vardathi2$Space+vardathi2$Abiotic+vardathi2$Plant
vardathi2$Spacer<-vardathi2$Space/vardathi2$tot
vardathi2$Abioticr<-vardathi2$Abiotic/vardathi2$tot
vardathi2$Plantr<-vardathi2$Plant/vardathi2$tot
vardathi3<-vardathi2%>%
  select(succession,taxon,Spacer,Abioticr,Plantr)%>%
  gather(source,var,Spacer:Plantr)

head(vardathi2)

vardat<-rbind(vardatlo,vardatme,vardathi)
vardat2<-rbind(vardatlo2,vardatme2,vardathi2)
vardat3<-rbind(vardatlo3,vardatme3,vardathi3)

temp<-vardat %>%
  mutate(succession=factor(succession,levels=c("Early","Mid","Late")))%>%
  group_by(succession,source)%>%
  summarise(mean=mean(var),se=std.error(var))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Meetings/ESA2018/varexp.pdf",width=7.5, height=5)#
ggplot(temp,aes(x=succession,y=mean,col=source))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=15),axis.text.x=element_text(size=rel(1.2)),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"),legend.position="none")+
  labs(x = "",y="Variance Explained")+
#  geom_line(stat = "identity", position = "identity",size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=1)+
  facet_wrap(~source,scales = "free")
dev.off()

vardat3$source<-plyr::revalue(vardat3$source,c("Spacer"="Space","Abioticr"="Abiotic","Plantr"="Plant"))
vardat3$source<-factor(vardat3$source,levels=c("Space","Abiotic","Plant"))

temp<-vardat3 %>%
  mutate(succession=factor(succession,levels=c("Early","Mid","Late")))%>%
  group_by(succession,source)%>%
  summarise(mean=mean(var),se=std.error(var))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Meetings/ESA2018/varexpsapnotfree.pdf",width=9.5, height=5)#
ggplot(temp,aes(x=succession,y=mean,col=source))+
  #theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=15),axis.text.x=element_text(size=rel(1.2)),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"),legend.position="none")+
  labs(x = "",y="Variance Explained")+
  #  geom_line(stat = "identity", position = "identity",size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=1)+
  facet_wrap(~source)#,scales = "free"
dev.off()



temp<-vardat3 %>%
  mutate(succession=factor(succession,levels=c("Early","Mid","Late")))%>%
  group_by(succession,source,taxon)%>%
  summarise(mean=mean(var),se=std.error(var))

ggplot(temp,aes(x=succession,y=mean,col=taxon))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=14),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  labs(x = "",y="Variance Explained")+
  #  geom_line(stat = "identity", position = "identity",size=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=1)+
  facet_wrap(~source*taxon,scales = "free")



colnames(modelhi$data$Y)
#lo:  nem:0 euk:1:12 bact:13:118  fun:119:139  
#me: nem:0 euk: 1:6  bact: 7:122 fun 123:133
#hi:  nem:1:4 euk:5:14 bact:15:115  fun:116:120

#effect of abiotic:
mean(variationPartlo[1:12,1])
mean(variationPartlo[13:118,1])
mean(variationPartlo[119:139,1])

mean(variationPartme[1:6,1])
mean(variationPartme[7:122,1])
mean(variationPartme[123:133,1])

mean(variationPart[1:4,1])
mean(variationPart[5:14,1])
mean(variationPart[15:115,1])
mean(variationPart[116:120,1])

#biotic (plant):
mean(variationPartlo[1:12,2])
mean(variationPartlo[13:118,2])
mean(variationPartlo[119:139,2])

mean(variationPartme[1:6,2])
mean(variationPartme[7:122,2])
mean(variationPartme[123:133,2])

mean(variationPart[1:4,2])
mean(variationPart[5:14,2])
mean(variationPart[15:115,2])
mean(variationPart[116:120,2])

#plot
mean(variationPart[,3])



###### R2 ######
#I don't understand the difference b/t variance partitioning and R2, the numbers are not the same. This calcualtes R2 <- 1 - ssRes/ssY. Ok I think the difference is that the total height of the bar (which is 100% in the variance partitining) should be scaled to this R2 value; because above is the contribution of each parameter to the predicted value and here is the difference b/t the predicted value and the observed value
#I think this doesn't work unless gaussian or probit
#Ymean <- apply(model$data$Y,2,mean)
R2 <- Rsquared(modello, averageSp=FALSE)
mean(R2)
mean(R2[1:16])
mean(R2[17:148])
mean(R2[149:172])

#plot(Ymean,R2,pch=19)

#checking that R2 makes sense. this doesn't work b/c I guess R2 is relative to a model with only an intercept. but it is good b/c bdenovo990996 has high correlations with other taxa
bdenovo990996 R2 is 0.49353501
summary(lm(hmscYe$bdenovo90996~decostand(hmscXe[,2],method="standardize")))$r.squared #r2 is .25 but that doesn't include a random plot effect, and I'm not sure if it is relevant to varipart b/c th random+fixed effects R2 adds up to 100, so it is relative??
str(summary(lm(hmscYe$bdenovo90996~1)))
test=data.frame(y=hmscYe$bdenovo90996,x=decostand(hmscXe[,2],method="standardize"),r=1:29)

require(MuMIn)
#The marginal R squared values are those associated with your fixed effects, the conditional ones are those of your fixed effects plus the random effects. 
m1<-lme(y~x,random=~1|r,data=test)
r.squaredGLMM(m1)
#marginal (fixed only) R2: .244, conditional (f+r): .907
#these r2 are not matching up with the ones from hmsc, i can't figure out why, unless it is because the hmsc r2 consider the effect of interactions with other species??? I could try fitting the hmsc model with more or fewer species to see if the R2 changes




###### Species co-occurrances #####

#complexity = linkage density = #interactions/taxon


###### Circle network ######

#lo
corMatlo <- corRandomEff(modello, cor = TRUE)
averageCorlo <- apply(corMatlo[, , , 1], 1:2, mean)
colMatlo <- matrix(NA, nrow = nrow(averageCorlo), ncol = ncol(averageCorlo))
colMatlo[which(averageCorlo > 0.9, arr.ind = TRUE)] <- "blue"
colMatlo[which(averageCorlo < -0.9, arr.ind = TRUE)] <- "red"

rownames(averageCorlo)<-paste(sapply(rownames(averageCorlo),function(x) substr(x,1,1)),1:length(rownames(averageCorlo)),sep="")
colnames(averageCorlo)<-paste(sapply(rownames(averageCorlo),function(x) substr(x,1,1)),1:length(rownames(averageCorlo)),sep="")

#taxon=c(rep("Eukaryota",12),rep("Bacteria",106),rep("Fungi",21))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Meetings/ESA2018/networklo.pdf",width=5, height=5)#
chordDiagram(averageCorlo, symmetric = TRUE,
             annotationTrack = c( "grid"),#"name",
             grid.col = c(rep("#57aa60",12),rep("#5f4686",106),rep("#ba9500",21)),col=colMatlo)
dev.off()

colMat <- matrix(0, nrow = nrow(averageCorlo), ncol = ncol(averageCorlo))
colMat[which(averageCorlo > 0.9, arr.ind = TRUE)] <- 1
colMat[which(averageCorlo < -0.9, arr.ind = TRUE)] <- 1
#sum colmat, subtract the diagonal, divide by 2, divide by number of species
(sum(colMat)-dim(colMat)[1])/2/dim(colMat)[1]


#me
corMatme <- corRandomEff(modelme, cor = TRUE)
averageCorme <- apply(corMatme[, , , 1], 1:2, mean)
colMatme <- matrix(NA, nrow = nrow(averageCorme), ncol = ncol(averageCorme))
colMatme[which(averageCorme > 0.9, arr.ind = TRUE)] <- "blue"
colMatme[which(averageCorme < -0.9, arr.ind = TRUE)] <- "red"

rownames(averageCorme)<-paste(sapply(rownames(averageCorme),function(x) substr(x,1,1)),1:length(rownames(averageCorme)),sep="")
colnames(averageCorme)<-paste(sapply(rownames(averageCorme),function(x) substr(x,1,1)),1:length(rownames(averageCorme)),sep="")

#taxon=c(rep("Eukaryota",6),rep("Bacteria",116),rep("Fungi",11))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Meetings/ESA2018/networkme.pdf",width=5, height=5)#
chordDiagram(averageCorme, symmetric = TRUE,
             annotationTrack = c( "grid"),#"name",
             grid.col = c(rep("#57aa60",6),rep("#5f4686",116),rep("#ba9500",11)),col=colMatme)
dev.off()

colMat <- matrix(0, nrow = nrow(averageCorme), ncol = ncol(averageCorme))
colMat[which(averageCorme > 0.9, arr.ind = TRUE)] <- 1
colMat[which(averageCorme < -0.9, arr.ind = TRUE)] <- 1
#sum colmat, subtract the diagonal, divide by 2, divide by number of species
(sum(colMat)-dim(colMat)[1])/2/dim(colMat)[1]


#hi
corMathi <- corRandomEff(modelhi, cor = TRUE)
averageCorhi <- apply(corMathi[, , , 1], 1:2, mean)
colMathi <- matrix(NA, nrow = nrow(averageCorhi), ncol = ncol(averageCorhi))
colMathi[which(averageCorhi > 0.9, arr.ind = TRUE)] <- "blue"
colMathi[which(averageCorhi < -0.9, arr.ind = TRUE)] <- "red"

rownames(averageCorhi)<-paste(sapply(rownames(averageCorhi),function(x) substr(x,1,1)),1:length(rownames(averageCorhi)),sep="")
colnames(averageCorhi)<-paste(sapply(rownames(averageCorhi),function(x) substr(x,1,1)),1:length(rownames(averageCorhi)),sep="")

#taxon=c(rep("Mesofauna",4),rep("Eukaryota",10),rep("Bacteria",101),rep("Fungi",5))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Meetings/ESA2018/networkhi.pdf",width=5, height=5)#
chordDiagram(averageCorhi, symmetric = TRUE,
             annotationTrack = c( "grid"),#"name",
             grid.col = c(rep("#d60062",4),rep("#57aa60",10),rep("#5f4686",101),rep("#ba9500",5)),col=colMathi)
dev.off()

colMat <- matrix(0, nrow = nrow(averageCorhi), ncol = ncol(averageCorhi))
colMat[which(averageCorhi > 0.9, arr.ind = TRUE)] <- 1
colMat[which(averageCorhi < -0.9, arr.ind = TRUE)] <- 1
#sum colmat, subtract the diagonal, divide by 2, divide by number of species
(sum(colMat)-dim(colMat)[1])/2/dim(colMat)[1]

#Legend
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Meetings/ESA2018/legend2.pdf")
plot(c(1,1),c(1,1))
legend("topright",c("Bacteria","Small eukaryotes","Fungi","Soil mesofauna"),pt.bg=c("#5f4686","#57aa60","#ba9500","#d60062"),bty="n",pch=22,cex=1.5)
legend("topleft",c("Positive","Negative"),col=c("#687dcb","#ce4d42"),lty=1,lwd=2,bty="n",cex=1.5)
dev.off()




###### Distance decay #####
vegdistlo <- vegdist(modello$data$Y) # Bray-Curtis
envdistlo <- vegdist(autoxyblo[,2:3], "euclid")
#sqrt((445132-445139)^2+(4433621-4433623)^2)
mantel(vegdistlo, envdistlo)#,method="kendall"
mantel(vegdistlo, log(envdistlo))
plot(envdistlo,vegdistlo)
abline(lm(vegdistlo~envdistlo))
plot(log(envdistlo),vegdistlo)

vegdistme <- vegdist(modelme$data$Y) # Bray-Curtis
envdistme <- vegdist(autoxybme[,2:3], "euclid")
#sqrt((445132-445139)^2+(4433621-4433623)^2)
mantel(vegdistme, envdistme)
mantel(vegdistme, log(envdistme))
plot(envdistme,vegdistme)
abline(lm(vegdistme~envdistme))
plot(log(envdistme),vegdistme)

vegdisthi <- vegdist(modelhi$data$Y) # Bray-Curtis
envdisthi <- vegdist(autoxybhi[,2:3], "euclid")
#sqrt((445132-445139)^2+(4433621-4433623)^2)
mantel(vegdisthi, envdisthi)
mantel(vegdisthi, log(envdisthi))
plot(envdisthi,vegdisthi)
abline(lm(vegdisthi~envdisthi))
plot(log(envdisthi),vegdisthi)



##### pcnm #####
pclo<-pcnm(envdistlo)
pclo$vectors[,1:2]
ordisurf(autoxyblo[,2:3], scores(pclo, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(autoxyblo[,2:3], scores(pclo, choi=2), bubble = 4, main = "PCNM 2")

pcme<-pcnm(envdistme)
ordisurf(autoxybme[,2:3], scores(pcme, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(autoxybme[,2:3], scores(pcme, choi=2), bubble = 4, main = "PCNM 2")

pchi<-pcnm(envdisthi)
ordisurf(autoxybhi[,2:3], scores(pchi, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(autoxybhi[,2:3], scores(pchi, choi=2), bubble = 4, main = "PCNM 2")






###### Correlations #####
#Confidence intervals for pairwise correlations, these distributions are very odd, they are bimodal with a lot at -1 and a lot at 1, so the mean is somewhere in the middle
corMat <- corRandomEff(model, cor = TRUE)
#when I did lo, I needed to set this
#rownames(corMat)<-colnames(model$data$Y)
#colnames(corMat)<-colnames(model$data$Y)
#hist(corMat[15,14 , , 1])
ltri <- lower.tri(apply(corMat[, , , 1], 1:2, quantile, probs = 0.025),diag=F) #originally diag=T
### Average
averageCor <- as.vector(apply(corMat[, , , 1], 1:2, mean)[ltri])
medianCor <- as.vector(apply(corMat[, , , 1], 1:2, median)[ltri])
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs = 0.025)[ltri])
corMat.975 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs=0.975)[ltri])
CICor <- cbind(corMat.025, corMat.975)
head(CICor)
#which(CI[,1]>0&CI[,2]>0)

#put labels on the pairwise correlations
rownames(corMat)
CorSp<-data.frame(sp1=rep(NA,length(averageCor)),sp2=rep(NA,length(averageCor)))
r=1
for(i in 1:(length(rownames(corMat))-1)){
  for(j in (i+1):length(rownames(corMat))){
    CorSp[r,1]<-rownames(corMat)[i]
    CorSp[r,2]<-rownames(corMat)[j]
    r=r+1
  }
}
head(CorSp)
tail(CorSp)
CorsCI<-cbind(CorSp,averageCor,CICor)
head(CorsCI)
CorsCIhi<-CorsCI

#these are the same as with the Cov matrix, so I guess that's ok, even though the distributions are kind of wonky
ind<-which(CorsCI$corMat.025<0&CorsCI$corMat.975<0|CorsCI$corMat.025>0&CorsCI$corMat.975>0)
length(ind)
CorsCI[ind,]
hist(CorsCI$averageCor)
unique(c(CorsCI$sp1[ind],CorsCI$sp2[ind]))
length(unique(c(CorsCI$sp1[ind],CorsCI$sp2[ind])))

ind<-which(CorsCI$averageCor>.9)
CorsCI[ind,]
ind<-which(CorsCI$averageCor<(-.9))
CorsCI[ind,]

hist(corMat["S6922370212b10cc1f286d15985037751","B9041166ce31053d58c3df9ec9ae26229" , , 1])


#looking at some of the correlations
#CorsCI3<-subset(CorsCI2,qval<.01)
#plot(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo145037"])
#abline(lm(hmscYe[,"bdenovo145037"]~hmscYe[,"bdenovo193772"]))
#summary(lm(hmscYe[,"bdenovo145037"]~hmscYe[,"bdenovo193772"]))


###### Plotting matrix #####
averageCor2 <- apply(corMat[, , , 1], 1:2, mean) #this is the full matrix, not just the lower triangle, for the mean it might not matter, but for CI it might? but no I don't think it does for the CI, because the CI is taking the confidence interval over the mcmc chain, it just does the calcualtion twice if you do it on the whole matrix
averageCor2[1:10,1:10]

rownames(averageCor2)<-paste(sapply(rownames(averageCor2),function(x) substr(x,1,1)),1:length(rownames(averageCor2)),sep="")
colnames(averageCor2)<-paste(sapply(rownames(averageCor2),function(x) substr(x,1,1)),1:length(rownames(averageCor2)),sep="")

corrplot(averageCor2, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200))


###### Circle network ######
corMat2 <- corRandomEff(modello, cor = TRUE)
averageCor <- apply(corMat2[, , , 1], 1:2, mean)
colMat <- matrix(NA, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor > 0.9, arr.ind = TRUE)] <- "red"
colMat[which(averageCor < -0.9, arr.ind = TRUE)] <- "blue"

rownames(averageCor)<-paste(sapply(rownames(averageCor),function(x) substr(x,1,1)),1:length(rownames(averageCor)),sep="")
colnames(averageCor)<-paste(sapply(rownames(averageCor),function(x) substr(x,1,1)),1:length(rownames(averageCor)),sep="")

chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = c(rep("grey",4),rep("blue",13),rep("red",124),rep("yellow",7)),col=colMat)#grid.col = c(rep("grey",4),rep("blue",13),rep("red",124),rep("yellow",7))

chordDiagram(averageCor, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey",col=colMat)

colMat <- matrix(0, nrow = nrow(averageCor), ncol = ncol(averageCor))
colMat[which(averageCor > 0.9, arr.ind = TRUE)] <- 1
colMat[which(averageCor < -0.9, arr.ind = TRUE)] <- 1
#sum colmat, subtract the diagonal, divide by 2, divide by number of species
(sum(colMat)-dim(colMat)[1])/2/dim(colMat)[1]

rownames(corMat2)
chordDiagram(averageCor[1:50,1:50], symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "gray",col=colMat, reduce=.01)#
?chordDiagram


##### Graphing with igraph #####
CorsCIhi$direction<-sign(CorsCIhi$averageCor)
dim(CorsCIhi)
ind<-which(CorsCIhi$corMat.025<0&CorsCIhi$corMat.975<0|CorsCIhi$corMat.025>0&CorsCIhi$corMat.975>0)
length(ind)
inputhi<-CorsCIhi[ind,]
dim(inputhi)
head(inputhi)
#inputhiv<-subset(edge_listsKS32no2b,qval<.05&trt=="hi")#
#vertexsizes1<-unique(data.frame(otu=c(as.character(inputhiv$taxa1),as.character(inputhiv$taxa2)),abun=c(inputhiv$ab1,inputhiv$ab2)))
graph1<-simplify(graph.edgelist(as.matrix(inputhi[,1:2]),directed=FALSE))
graph1$layout <- layout_in_circle
#verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
#colnames(verticesgraph1)<-"otu" #
#colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph1<-ifelse(verticesgraph1$otu%in%hubshi,8,4)
plot(graph1,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(inputhi$direction==-1,"blue","red"))
#plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2  vertex.label=as.character(colorgraph1$orders)  


ltri<-lower.tri(array(NA,dim=c(dim(corMat2)[1],dim(corMat2)[2])),diag=F)
averageCor <- as.vector(rowSums(corMat2[,,,1],dims=2)[ltri])/dim(corMat2)[3] #faster but this is a vector not a matrix, but I can use it to count
length(which(averageCor>.7|averageCor<(-.7)))
length(which(averageCor>.8|averageCor<(-.8)))
length(which(averageCor>.9|averageCor<(-.9)))
which(averageCor<(-.9))
which(averageCor>(.9))
hist(averageCor)

averageCor2<-data.frame(taxon1=rownames(averageCor),averageCor)
averageCor3<-gather(averageCor2,taxon2,cor,bdenovo195709:idenovo18189) #not labeling like I want
head(averageCor3)
subset(averageCor3,averageCor3$cor>.9&averageCor3$cor<1)

#checking, for species 2 and 3 intercept and slope should be significant
plot(formdata$X[,2],formdata$Y[,2])
abline(a=-2.571260,b=0.5421437)
abline(a=-7.1447,b=1.4906,col=2)
summary(lm(formdata$Y[,1]~formdata$X[,2]))
summary(lm(formdata$Y[,4]~formdata$X[,2]))


#trying just calculating intercept on simulated data, everything checks out, confidence intervals are correct and standard error in lm is the same as the calucated SE and you can calculate the CI from the standard error in lm()
hmscYd<-as.matrix(data.frame(y1=rnorm(1000,mean=1,sd=.2),y2=rnorm(1000,mean=1,sd=.2)))
hmscXd<-as.matrix(rep(1,1000))

formdata <- as.HMSCdata(Y = hmscYd, X = hmscXd, interceptX = F, scaleX=T)
model <- hmsc(formdata, family = "gaussian", niter = 10000, nburn = 1000, thin = 10)
0.99505+1.96*.2/sqrt(1000)
summary(lm(hmscYd[,1]~1))
.2/sqrt(1000)
sd(hmscYd)/sqrt(1000)
0.995050-0.006077*1.96
0.995050+0.006077*1.96

#trying to calculate likelihood
-sum(dnorm(formdata$Y[,2],mean=-100+1.4906*formdata$X[,2],sd=1,log=T)) 
-sum(dnorm(formdata$Y[,2],mean=-7.1447+1.4906*formdata$X[,2],sd=1,log=T)) #I don't know the sd, but I think it doesn't matter if I'm only looking at relative differences in loglik
-sum(dnorm(formdata$Y[,2],mean=-2.5712604+0.5421437*formdata$X[,2],sd=1,log=T))
#in terms of likelihood, the lm() model estimates are better




##### Covariances #####
#Confidence intervals for pairwise covariances. Covariances mcmc are much more normally distributed, so I feel more comfortable using them for calculating z statsitics and p values. They aren't perfectly normally distributed (so it would probably be more accurrate to use some kind of on sided z-statistic) but it is probably ok and would be easier than explaining all that.
corMat <- corRandomEff(model,cor=F) #cor=F gives the covariance. These covariances are not actualy fit by the model, they are calcuated post hoc based on the latent variables at each iteration, start 3:05, end 3:15ish

#ltri <- lower.tri(apply(corMat[, , , 1], 1:2, quantile, probs = 0.025),diag=F)#start 3:20, end waited until 4:20 but it hadn't stopped so I stopped it #used to be diat=T
#apply function takes a ridiculousy long time, need to vectorize things
ltri<-lower.tri(array(NA,dim=c(dim(corMat)[1],dim(corMat)[2])),diag=F) #does the same thing and is way way faster

### Average
#averagec <- as.vector(apply(corMat[, , , 1], 1:2, mean)[ltri])
averagec<- as.vector(rowSums(corMat[,,,1],dims=2)[ltri])/dim(corMat)[3] #faster
medianc <- as.vector(apply(corMat[, , , 1], 1:2, median)[ltri])
head(cbind(averagec,medianc))
plot(corMat[1,2,,1],type="l")
hist(corMat[1,2,,1]) #dims 1 and 2 are the matrix species (35) by species, 3 is the iteration, 4 is the plot level (I only have 1 level in that dimension), bdenovo32932 and bdenovo195709
plot(density(corMat[1,2,,1]))
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs = 0.025)[ltri]) 
corMat.975 <- as.vector(apply(corMat[, , , 1], 1:2, quantile,probs=0.975)[ltri])
#try something like this to make it faster???
#a<-array(data=1:90,dim=c(3,3,10,1))
#colSums(matrix(a,nrow=10,ncol=9,byrow=T))/dim(a)[3]

CICov <- cbind(corMat.025, corMat.975)
head(CICov)
# plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI), type = "n", xlab = "", main = "cov(paramLatent[[1, 1]])")
# abline(h = 0, col = "grey")
# arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2], code = 3, angle = 90, length = 0.05)
# points(1:nrow(CI), average, pch = 15,cex = 1.5)

#put labels on the pairwise correlations
rownames(corMat)
CorSp<-data.frame(sp1=rep(NA,length(averagec)),sp2=rep(NA,length(averagec)))
r=1
for(i in 1:(length(rownames(corMat))-1)){
  for(j in (i+1):length(rownames(corMat))){
    CorSp[r,1]<-rownames(corMat)[i]
    CorSp[r,2]<-rownames(corMat)[j]
    r=r+1
  }
}
head(CorSp)
tail(CorSp)
CovsCI<-cbind(CorSp,averagec,CICov)
head(CovsCI)

ind<-which(CovsCI$corMat.025<0&CovsCI$corMat.975<0|CovsCI$corMat.025>0&CovsCI$corMat.975>0)
CovsCI[ind,]#there are 960 (hi) that don't overlap 0!!! need to see if my code below is correct
length(ind)
unique(c(CovsCI$sp1[ind],CovsCI$sp2[ind]))




#upshot
#mcmc is very sensitive to the range in the y variables, if the numbers are too low, the CI will be huge b/c somehow it takes huge jumps in the mcmc - solution: scale the Y
#mcmc is also sensitive to the absolute value of the x variables (but I did this on my original data trial so this was not the problem then). even if it ranges from 6-8 the mcmc gives non-optimal results compared if the range is -1 to 1 


















#####Old code#####

# Covariances
#translate to p values
#P from CI for a difference
#If the upper and lower limits of a 95% CI are u and l respectively:
#  1 calculate the standard error: SE = (u − l)/(2*1.96)
#  2 calculate the test statistic: z = Est/SE
#  3 calculate the P value2: P = exp(−0.717×z − 0.416×z2). #I changed to using pnorm(), more accurate
#Even with a model with only an intercept, the lowest p value is .001, which is a qval of .8, so nothing is significant
#Thoughts - since these are all estimated simulteneously, I don't know if I really need to account for multiple comparions. Also, I could mayb just use the cutoff "significance is when th 95% CI don't cross zero" this is the bayesian credible interval, right? yes, and a 95% CI should correspond to a 0.05 p value
CovsCI$SE<-(CovsCI$corMat.975-CovsCI$corMat.025)/(2*1.96)
CovsCI$z<-CovsCI$averagec/CovsCI$SE
CovsCI$absz<-abs(CovsCI$z)
#CovsCI$Pone<-pnorm(-CovsCI$absz) #one tailed test, don't use
CovsCI$P<-2*pnorm(-CovsCI$absz) #two tailed test. it should be a two tailed test b/c it could be higher or lower than 0 (and that's what lme does for example)
#CovsCI$Pold<-exp(-.717*CovsCI$absz-.416*CovsCI$absz^2)
which(CovsCI$Ptwo<.1)
CovsCI$qval<-p.adjust(CovsCI$P,method="fdr")#
sort(CovsCI$qval)
sort(CovsCI$P) 
CovsCI2<-subset(CovsCI,averagec>0)
which(CovsCI2$P<.1)
CovsCI2[37,]#8412,3112
CovsCI2[135,]#8412,3112

CovsCI[which(CovsCI$sp1=="bdenovo193772"&CovsCI$sp2=="bdenovo20595"),]
hist(corMat["bdenovo193772","bdenovo20595",,1])
hist(corMat["bdenovo117183","bdenovo85595",,1])
#so I think my p values are high b/c I'm dealing with slighly non normal data

#looking into one species pair with strong positive covariance
bdenovo88234 bdenovo137544
bdenovo184998 bdenovo117183
bdenovo198139 bdenovo51656
bdenovo193772 bdenovo20595
hmscXe
hmscYe[,"bdenovo193772"]
hmscYe[,"bdenovo20595"]
cov(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"])
cor(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"])
#library(Hmisc)
#rcorr(cbind(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"]))$P #same P value as regression model
plot(hmscYe[,"bdenovo193772"],hmscYe[,"bdenovo20595"])
summary(lm(hmscYe[,"bdenovo193772"]~hmscYe[,"bdenovo20595"]))
abline(lm(hmscYe[,"bdenovo193772"]~hmscYe[,"bdenovo20595"]))
lm1<-lm(hmscYe[,"bdenovo193772"]~hmscXe[,"pH"])
lm2<-lm(hmscYe[,"bdenovo20595"]~hmscXe[,"pH"])
plot(resid(lm1),resid(lm2))
abline(lm(resid(lm1)~resid(lm2)))
summary(lm(resid(lm1)~resid(lm2)))





# Correlations
#translate to p values
#P from CI for a difference
#If the upper and lower limits of a 95% CI are u and l respectively:
#  1 calculate the standard error: SE = (u − l)/(2*1.96)
#  2 calculate the test statistic: z = Est/SE
#  3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).
#This only works when the CIs are symmetric around the estimate, as shown below, many of my CIs are not asymmetric b/c correlations can't be larger than 1 so the CI is squinched at that side
CorsCI$SE<-(CorsCI$corMat.975-CorsCI$corMat.025)/(2*1.96)
CorsCI$z<-CorsCI$averageCor/CorsCI$SE
CorsCI$absz<-abs(CorsCI$z)
CorsCI$P<-exp(-.717*CorsCI$absz-.416*CorsCI$absz^2)
which(CorsCI$P<.001)
CorsCI$qval<-p.adjust(CorsCI$P,method="fdr")
CorsCI2<-subset(CorsCI,averageCor>0)

#trying it by taking the "one-sided" p value
CorsCI$SE2<-ifelse(CorsCI$averageCor>0,(CorsCI$averageCor-CorsCI$corMat.025)/(1*1.96),(CorsCI$corMat.975-CorsCI$averageCor)/(1*1.96)) #this way is still doing a "two tailed test" b/c I'm using 1.96 sd away from the mean, if I were doing a one tailed test I would use 1.645 sd. In other words this way is like doing a one-tailed test with an alpha=0.025, and I'm using a one-tailed test simply b/c I don't have the proper data (i.e. asymmetrical CIs) to calcualted a true two tailed test. however, it seems that the z statistic calculation is based on the type of confidence interval you have, so if you have 95% CI, then you should use 1.96, it doesnt make sense to use a different number. I think the alpha cutoff is for translating the z value to a pvalue and giving you a cut off. Another thought is that I maybe should be using the median instead of the mean here 
CorsCI$z2<-CorsCI$averageCor/CorsCI$SE2
CorsCI$absz<-abs(CorsCI$z)
CorsCI$P<-exp(-.717*CorsCI$absz-.416*CorsCI$absz^2)
which(CorsCI$P<.001)
CorsCI$qval<-p.adjust(CorsCI$P,method="fdr")
CorsCI2<-subset(CorsCI,averageCor>0)


length(which(CorsCI2$P<.05))
length(which(CorsCI2$qval<.01))
hist(CorsCI2$averageCor[which(CorsCI2$qval<.01)])







# Mixing object, I don't think this is useful. I don't totally understand latent variables- they somehow help with fitting a model with so many parameters (the covoariance matrix has 630 parameters in this example case accounting for all the species pairwise correaltions.) when you run a model with no env factors, the latent variables are essentially ordination axes
mixing <- as.mcmc(model, parameters = "paramLatent")
### Draw trace and density plots for all combination of paramters
dim(mixing) #900 x 140, I don't undrsand what this is, what the 140 refers to, it looks like there are about 3-5 latent vaiables per species
plot(mixing[,140])
mean(mixing[,1])
hist(mixing[,1])
cbind(corMat[1,2,,1],mixing[,1],corMat[1,2,,1]/mixing[,1])
### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing[[2]])
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)

