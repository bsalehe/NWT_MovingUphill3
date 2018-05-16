#Functional redundancy

#smaller exponent a means larger functional redundancy


#All genes

#open the txt file in excel and delete the # next to #OTU ID and put OTUID together as one word, and delete the first line
metagenome<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/closedRef_forPICRUSt_97db_97sim/PICRUSt/otu_table21.txt",header=T)
row.names(metagenome)<-metagenome$OTUID
metagenome$OTUID<-NULL
head(metagenome)

metagenomet<-as.data.frame(t(as.matrix(metagenome)))
metagenomet[1:9,1:9]
metagenomet$X.SampleID<-row.names(metagenomet)

metagenomet2<-merge(metagenomet,data.frame(X.SampleID=biogeo6$X.SampleID,lomehi=biogeo6$lomehi))
metagenomet2[1:9,1:9]
metagenomet2$lomehi
metagenomet2$orthorich<-rowSums(metagenomet2[,2:6910]>0)
dim(metagenomet2)
metagenomet2$orthodiv<-vegan::diversity(metagenomet2[,2:6910],index="shannon")

m1<-metagenomet2%>%
  group_by(lomehi)%>%
  summarise(mean=mean(orthorich),se=std.error(orthorich))

ggplot(m1,aes(x=lomehi,y=mean,group=lomehi))+
  labs(x = "",y="Ortho shannon diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)




otutablem<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/closedRef_forPICRUSt_97db_97sim/exported-feature-table/otu_table.txt",header=T)
head(otutablem)

#I rarefied to 10825 so this is actually calculation of how many data sequences were matched (rather than repset) = 79.1%
mean(colSums(otutablem[,2:95]))/10825

row.names(otutablem)<-otutablem$OTUID
otutablem$OTUID<-NULL
head(otutablem)

otutablemt<-as.data.frame(t(as.matrix(otutablem)))
otutablemt[1:9,1:9]
otutablemt$X.SampleID<-row.names(otutablemt)

otutablemt2<-merge(otutablemt,data.frame(X.SampleID=biogeo6$X.SampleID,lomehi=biogeo6$lomehi))
otutablemt2[1:9,1:9]
otutablemt2$lomehi
otutablemt2$sprich<-rowSums(otutablemt2[,2:4632]>0)
dim(otutablemt2)
otutablemt2$spdiv<-vegan::diversity(otutablemt2[,2:4632],index="shannon")

colSums(otutablem[,2:95]>0)

plot(otutablemt2$sprich,metagenomet2$orthorich,col=otutablemt2$lomehi)
curve(exp(8.20449)*x^0.04866,add=T,col=2)
curve(exp(7.94693)*x^0.08837,add=T,col=3)
curve(exp(8.28269)*x^0.03516,add=T,col=1)

temp<-data.frame(X.SampleID=otutablemt2$X.SampleID,lomehi=otutablemt2$lomehi,sprich=otutablemt2$sprich,orthorich=metagenomet2$orthorich)

summary(lm(log(orthorich)~log(sprich),data=subset(temp,lomehi=="lo")))
summary(lm(log(orthorich)~log(sprich),data=subset(temp,lomehi=="me")))
summary(lm(log(orthorich)~log(sprich),data=subset(temp,lomehi=="hi")))




###### Rarefaction curve #####
#Do resampling to get rarefaction curve (Ortholog richness vs species richness)

#take out the taxonomy columns b/c they give an error if they are blank
awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' Galaxy3-[Metagenome_Contributions_by_Higher_Category_on_data_2]Metabolism.txt > otufile3.txt

metaotu<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/closedRef_forPICRUSt_97db_97sim/PICRUSt/otufile3.txt",header=T)
head(metaotu)

otus<-as.numeric(rownames(otutablem))#  4631 taxa
genesmet<-unique(metaotu$Gene) #2933 genes

length(unique(metaotu$OTU)) #they are all represented in the metaotu file

otutablemt2lo<-otutablemt2%>%
  filter(lomehi=="lo")%>%
  select(-lomehi,-sprich,-spdiv,-X.SampleID)%>%
  select(which(colSums(.) > 0))
otuslo<-as.numeric(colnames(otutablemt2lo))

#otutablemt2lo[1:10,1:10]
#colSums(otutablemt2lo)

temp<-data.frame(X.SampleID=biogeo6$X.SampleID,lomehi=biogeo6$lomehi)
templo<-temp$X.SampleID[which(temp$lomehi=="lo")]

metaotu2<-metaotu%>%
  filter(Sample%in%templo)

metaotu3<-unique(metaotu2[,c(1,3)])

length(unique(metaotu3[,c(1)]))
#lo 2699, me 2831, hi 2775

dim(metaotu3)
length(otuslo)

outlo<-data.frame(speciesrich=seq(from=1, to=2241,by=50),orthologrich=rep(NA,length(seq(from=1, to=2241,by=50))))

# for(i in 1:length(seq(from=1, to=2241,by=100))){
#   l<-outlo[i,1]
#   inds<-sample(otuslo,l)
#   ind<-which(metaotu3$OTU%in%inds)
#   temp<-metaotu3[ind,]
#   outlo[i,2]<-length(unique(temp$Gene)) 
# }

for(i in 1:length(seq(from=1, to=2241,by=50))){
  l<-outlo[i,1]
  if(l==1){
    temp1<-data.frame(n=rep(l,20))
    temp2<-data.frame(apply(temp1,1,FUN=function(x){sample(otuslo,x)}))
    ind<-apply(temp2,1,FUN=function(x){which(metaotu3$OTU%in%x)})
    temp<-lapply(ind,function(x){metaotu3[x,]})
    temp3<-lapply(temp,function(x){length(unique(x$Gene))})
    outlo[i,2]<-median(as.numeric(temp3))
  }
  else{
    temp1<-data.frame(n=rep(l,20))
    temp2<-apply(temp1,1,FUN=function(x){sample(otuslo,x)})
    ind<-apply(temp2,2,FUN=function(x){which(metaotu3$OTU%in%x)})
    temp<-lapply(ind,function(x){metaotu3[x,]})
    temp3<-lapply(temp,function(x){length(unique(x$Gene))})
    outlo[i,2]<-median(as.numeric(temp3))
  }
}

outall<-outlo
outall<-rbind(outall,outlo)
outall<-rbind(outall,outlo)

outall$lomehi<-rep(c("lo","me","hi"),each=length(outall$speciesrich)/3)
outall$lomehi<-factor(outall$lomehi,levels=c("lo","me","hi"))
outall$col<-ifelse(outall$lomehi=="lo",1,ifelse(outall$lomehi=="me",2,3))


ggplot(outall, aes(x=speciesrich,y=orthologrich,col=lomehi))+
  labs(x = "OTU richness",y="Gene richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_point() +
  geom_line() 

plot(outall$speciesrich,outall$orthologrich,col=outall$col,xlim=c(0,400))
curve(exp(6.83370)*x^0.14873,add=T,col=2)
curve(exp(6.852773)*x^0.149920,add=T,col=3)
curve(exp(6.83381)*x^0.15011,add=T,col=1)


summary(lm(log(orthologrich)~log(speciesrich),data=subset(outall,lomehi=="lo")))
summary(lm(log(orthologrich)~log(speciesrich),data=subset(outall,lomehi=="me")))
summary(lm(log(orthologrich)~log(speciesrich),data=subset(outall,lomehi=="hi")))






##### Selecing only metabolism genes for ortholog richness and diversity #####

metagenomemet<-metagenome%>%
  mutate(gene=rownames(.))%>%
  filter(gene%in%genesmet)
rownames(metagenomemet)<-metagenomemet$gene
metagenomemet$gene<-NULL

head(metagenomemet)

metagenomemett<-as.data.frame(t(as.matrix(metagenomemet)))
metagenomemett[1:9,1:9]
metagenomemett$X.SampleID<-row.names(metagenomemett)

metagenomemett2<-merge(metagenomemett,data.frame(X.SampleID=biogeo6$X.SampleID,lomehi=biogeo6$lomehi))
metagenomemett2[1:9,1:9]
metagenomemett2$lomehi
metagenomemett2$orthorich<-rowSums(metagenomemett2[,2:2934]>0)
dim(metagenomemett2)
metagenomemett2$orthodiv<-vegan::diversity(metagenomemett2[,2:2934],index="shannon")

m1<-metagenomemett2%>%
  group_by(lomehi)%>%
  summarise(mean=mean(orthorich),se=std.error(orthorich))

m1$lomehi<-factor(m1$lomehi,levels=c("lo",'me','hi'))

ggplot(m1,aes(x=lomehi,y=mean,group=lomehi))+
  labs(x = "",y="Gene richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)

m1<-gls(orthorich~lomehi,data=metagenomemett2)
anova(m1,type="marginal")
