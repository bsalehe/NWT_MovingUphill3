##Diversity

#pd16S<-pd(as.matrix(datBacr3fotu[,-c(1:31)]),phy_tree(datBac3f),include.root=TRUE)
pdBac$X.SampleID<-rownames(pdBac)
pdBac2<-merge(pdBac,biogeo6,"X.SampleID")
pdBac2$type<-"1Bacteria"

#colnames(richBac)<-"PD"
#richBac$SR<-richBac$PD
#richBac$X.SampleID<-rownames(richBac)
#richBac2<-merge(richBac,biogeo6,"X.SampleID")
#richBac2<-merge(richBac,datBacS3otu[,1:31],"X.SampleID")
#richBac2$type<-"1Bacteria"

colnames(richITS)<-"PD"
richITS$SR<-richITS$PD
richITS$X.SampleID<-rownames(richITS)
richITS2<-merge(richITS,biogeo6,"X.SampleID")
richITS2$type<-"2Fungi"

pdEukS$X.SampleID<-rownames(pdEukS)
pdEukS2<-merge(pdEukS,biogeo6,"X.SampleID")
pdEukS2$type<-"3Small Eukaryotes"

pdEukN$X.SampleID<-rownames(pdEukN)
pdEukN$X.SampleID<-gsub(pattern = "N", replace = "S", x = pdEukN$X.SampleID)
pdEukN2<-merge(pdEukN,biogeo6,"X.SampleID")
pdEukN2$type<-"4Soil Mesofauna"


richdata<-rbind(pdBac2,richITS2,pdEukS2,pdEukN2)

ggplot(richdata,aes(x=log10(Plant_Dens+1),y=SR))+# as.numeric(fert),color=species
  labs(x="Plant density",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")

#could log transform if I wanted to, it looks less skewed
ggplot(richdata,aes(x=log10(pca1+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")


richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(PD),se=std.error(PD))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/diversitybysuccessionalstage.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()

anova(lm(PD~lomehi,data=pdBac2))
anova(lm(PD~lomehi,data=richITS2))
anova(lm(PD~lomehi,data=pdEukS2))
anova(lm(PD~lomehi,data=pdEukN2))



ggplot(pdBac2,aes(x=log10(Plant_Dens+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Plant density",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black")

summary(lm(PD~log10(Plant_Dens+1),data=pdBac2))

ggplot(pdBac2,aes(x=log(pca1+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black")

summary(lm(PD~log(pca1+1),data=pdBac2))


richmeans<-pdBac2%>%
  group_by(lomehi)%>%
  summarise(mean=mean(PD),se=std.error(PD))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

anova(lm(PD~lomehi,data=richBac2))

ggplot(richmeans,aes(x=lomehi,y=mean,group=1))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)







