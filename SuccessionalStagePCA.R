#Successional stage PCA


##### biogeochemistry #####
biogeo<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Biogeochemistry/CN_enzymes_pH_moisture_whc_Niwot2015t.csv")
head(biogeo)
dim(biogeo)
#note on dataset - I updated this because for the data from microbial biomass (gravimetric moisture, IN, DOC, etc), sample 33 and 34 were mixed up. origianlly sample 33 was missing from the first dataset dorota gave me (CN_enzymes_pH_moisture_whc_Niwot2015.xlsx), and the nubmers fom 34 were really the values for 33. Then I got the All Data file from Dorota and this cleared it up b/c it had both samples 33 adn 34 in it, I also double checked with the Biogeochemistry_Niwot_2015.xlsx file where the calculations were done.

#merge with one of the mapping files
biogeo2<-merge(biogeo,datEukS3otu[,1:31])
#cbind(biogeo2$moisture,biogeo2$WHC)
head(biogeo2)

#I could replace all negative numbers with 0, not sure if I should do this, since they are relative, but it doesn't really make sense to have negative microbial biomass or inorganic N
#biogeo2[biogeo2<0]<-0



##### plant cover #####
#get plant cover data to use as "light" and merge with biogeo
plantcov<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Plants/Niwot_MovingUpHill_plots2015.csv")
names(plantcov)[1]<-"Sample_name"
plantcov$plantcover<-(plantcov$MOSS+plantcov$VEG)/plantcov$TOTAL
plantcov$plantcover[which(plantcov$Sample_name==64)]<-0 #this I'm 100% sure has zero plants, however bare and rock were recorded as 0s. I looked up the sample data sheet and it looks like Sam and I did not sample it because a soil smaple was not taken in 2007, however since it had 0 plants we added that plot and forgot to do cover on it.
head(plantcov)

#take out some columns that I don't want to confuse with mapping file columns
plantcov2<-plantcov%>%select(Sample_name,X,Y,plantcover)

head(biogeo2)
biogeo3<-merge(biogeo2,plantcov2,"Sample_name") #with merge, they got put in the right order
rownames(biogeo3)<-biogeo3$X.SampleID
head(biogeo3)
#write.csv(biogeo3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/biogeo3.csv",row.names=F)



##### PCA #####
#Use these plant variables: they include mosses/liverworts but not lichen
#Plant_Dens
#Plant_Div

#datITSS3otu3$X

biogeo4<-biogeo3%>%
  select(Plant_Dens,Plant_Div,TC,TN,NH4,NO3,MicC,MicN,pH,WHC,moisture,snowdepth,elevation,plantcover)#,IN,DOC,DON
ind<-which(is.na(rowSums(biogeo4))==F)
biogeo5<-biogeo4[ind,]
dim(biogeo5)

mypca<-rda(biogeo5,scale=T,na.action=na.omit)
plot(mypca)
summary(mypca)

succession<-(-scores(mypca)$sites[,1])

succession2<-data.frame(pca1=succession,X.SampleID=names(succession))
succession3<-succession2[order(succession2$pca1),]
succession3$lomehi<-rep(c('lo','me','hi'),each=25)

plot(biogeo2$snowdepth,biogeo2$Plant_Dens)
plot(biogeo2$snowdepth,biogeo2$pH)

summary(lm(biogeo3$Plant_Dens~biogeo3$WHC))
summary(lm(biogeo3$Plant_Dens~biogeo3$moisture))


###### Merge back with biogeo3 #####
biogeo6<-merge(biogeo3[,-which(names(biogeo3)=="lomehi")],succession3)
head(biogeo6)
cbind(biogeo6$plantcover,biogeo6$Plant_Dens,biogeo6$lomehi)
#it is kind of strange that some of the plots with 0 plants were classified as medium, but there is not much I can do if that's how the pca shakes out. i'll just have to see what the networks look like

