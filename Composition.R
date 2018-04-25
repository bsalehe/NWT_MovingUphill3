###### Change in relative abundance ######
datEukS3k2
datEukN3k2
datBacS3k2
datITSS3k2

#change a colname, it can't have a dash in it
colnames(datBacS3k2)[78]<-"WPS.2"
names(which(colSums(datBacS3k2[,32:81])>2))
relBac<-datBacS3k2 %>% 
  dplyr::select(Sample_name,Acidobacteria,Actinobacteria,Bacteroidetes,Cyanobacteria,Gemmatimonadetes,Heterotrophic_Chloroflexi,Heterotrophic_Planctomycetes,Heterotrophic_Proteobacteria,Verrucomicrobia,WPS.2) %>%
  gather(Taxa,abun,Acidobacteria:WPS.2) %>%
  mutate(type="A. Bacteria")

names(which(colSums(datITSS3k2[,32:42])>.5))
relITS<-datITSS3k2 %>% 
  dplyr::select(Sample_name,Ascomycota,Basidiomycota,Glomeromycota,Mortierellomycota) %>%
  gather(Taxa,abun,Ascomycota,Basidiomycota,Glomeromycota,Mortierellomycota) %>%
  mutate(type="B. Fungi")

names(which(colSums(datEukS3k2[,32:46])>4))#,Nonphotosynthetic_Excavata,Photosynthetic_Stramenopiles,
relEukS<-datEukS3k2 %>% 
  dplyr::select(Sample_name,Alveolata,Archaeplastida,Photosynthetic_Stramenopiles,Rhizaria) %>%
  gather(Taxa,abun,Alveolata,Archaeplastida,Photosynthetic_Stramenopiles,Rhizaria) %>%
  mutate(type="C. Small Eukaryotes")

names(which(colSums(datEukN3k2[,32:39])>1))
relEukN<-datEukN3k2 %>% 
  dplyr::select(Sample_name,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  gather(Taxa,abun,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  mutate(type="D. Soil mesofauna")

relALL1<-rbind(relBac,relITS,relEukS,relEukN)#
head(relALL1)

#merge with biogeo6 to get pca1
relALL<-merge(relALL1,biogeo6,"Sample_name")
head(relALL)

#plotdata<-relALL %>%
#  mutate(typeTaxa=paste(type,Taxa)) %>%
#  group_by(Taxa,lomehi,type,typeTaxa) %>%
#  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
#  #%>%filter(mean_abun>.04)

#this was weird, maybe something changed in ggplot or dplyr because the colors were messing up and it was listing the legend in alfabetical order by taxa rather than the order in the "plotdata" dataframe. the workaroudn was to set the levels of plotdata$Taxa so they were correct
plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  group_by(typeTaxa,Taxa,lomehi,type) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

as.data.frame(plotdata)
plotdata$lomehi<-factor(plotdata$lomehi,levels=c("lo","me","hi"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroupsR.pdf",width=6.5,height=4.3)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()

#10 bacteria, 4 fungi, 4 small euks, 4 large euks
mycols<-c("#4BC366",#light green
          "#D9A125",#yellow
          "#6F94DE",#light blue
          "#B4405E",#red
          "#D185E0",
          "#ff99a4",#light pink,
          "#659125",#green
          "#cf6f23",#orange
          "#5C426C",#dark purple
          "#6768A3",#medium blue last bact
          "#D9A125",#yellow
          "#B4405E",#red
          "#659125",
          "#6768A3",
          "#6F94DE",
          "#5C426C",#dark purple
          "#D185E0",#light purple
          "#cf6f23",#orange
          "#B4405E", #red
          "#4BC366", #light green
          "#5C426C", #dark purple
          "#D9A125") #yellow

                    "#659125", #darker green
          "#D185E0", #light purple
          "#6768A3", #medium purple
          "#D9A125", #yellow
          "#6F94DE") #blue)

#scatter plots
head(relALL)

plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa))
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroupsBFSLENscatter.pdf",width=6.5,height=6)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=log10(Plant_Dens+1),y=abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  scale_color_manual(values=mycols) +
  geom_point(size=.2)+
  geom_smooth(method=lm,se=F,size=.8) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()


#Doing anova on all of the above taxa groups
length(unique(relALL$Taxa))
anovaoutput<-data.frame(Taxa=rep(NA,22),F=rep(NA,22),P=rep(NA,22))
for(i in 1:22){
  current.taxa<-unique(relALL$Taxa)[i]
  temp<-relALL %>%
    filter(Taxa==current.taxa)
  mod<-anova(lm(abun~lomehi,data=temp))
  anovaoutput[i,1]<-current.taxa
  anovaoutput[i,2]<-mod$`F value`[1]
  anovaoutput[i,3]<-mod$`Pr(>F)`[1]
}
anovaoutput$qval<-p.adjust(anovaoutput$P,method="fdr")
anovaoutput$qval<-format(anovaoutput$qval,scientific=F)

anovaoutput$Taxa<-factor(anovaoutput$Taxa,levels=unique(plotdata$Taxa))
anovaoutput[order(anovaoutput$Taxa),]






#90 samples, relative abundance otu table, grouped by "kingdom" aka "labels"
datBacS3k2
#dim(comm.dataBac)
#comm.dataBac[1:10,1:32]

#colnames(datBac3fk2)[78]<-"WCHB1.60"
names(which(colSums(datBacS3k2[,32:74])>3))
relBac<-datBacS3k2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Acidobacteria,Actinobacteria,Bacteroidetes,Heterotrophic_Chloroflexi,Cyanobacteria,Planctomycetes,Proteobacteria,Verrucomicrobia) %>%
  gather(Taxa,abun,Acidobacteria:Verrucomicrobia) %>%
  mutate(type="A. Bacteria")

plotdata<-relBac %>%
  group_by(Taxa,lomehi) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
plotdata$lomehi<-factor(plotdata$lomehi,levels=c("lo","me","hi"))
data.frame(plotdata)

ggplot(plotdata,aes(x=lomehi,y=mean_abun,color=Taxa)+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  #theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)

mycols<-c("#4BC366",
          "#D9A125",
          "#659125",
          "#6768A3",
          "#5C426C",
          "#D185E0",
          "#6F94DE",
          "#B4405E")
            
ggplot(plotdata,aes(x=lomehi,y=mean_abun,color=Taxa,group=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))















mynmds<-metaMDS(comm.dataBac[,32:3430],dist="bray",trymax = 1000)
#old lomehi
col=ifelse(comm.dataBac$lomehi=="lo","lightblue",NA)
col[which(comm.dataBac$lomehi=="me")]<-"dodgerblue"
col[which(comm.dataBac$lomehi=="hi")]<-"darkblue"

plot(scores(mynmds),col=col,pch=21,bg=col)#-scores(mynmds)[,1],scores(mynmds)[,2]
ordiellipse(mynmds,groups=col,col=c("darkblue","dodgerblue","lightblue"),conf=.99999,kind="se",lwd=2)#
legend("bottomright",c("Early","Mid","Late"),col=c("#ab3b57","#5268cb","#639e51"),pch=21,pt.bg=c("#ab3b57","#5268cb","#639e51"),lty=1,bty="n")


#contains all data plus biogeochemistry and plant cover, 75 samples
#biogeo info 53
#N 143
#S 1124
#Bact 3399
#ITS 1122
ordi.bact<-cbind(comm.bio[,1:53],comm.bio[,1321:4719])

mynmds<-metaMDS(ordi.bact[,54:3452],dist="bray",trymax = 1000)
#new lomehi
col=ifelse(ordi.bact$lomehi=="lo","lightblue",NA)
col[which(ordi.bact$lomehi=="me")]<-"dodgerblue"
col[which(ordi.bact$lomehi=="hi")]<-"darkblue"

plot(scores(mynmds),col=col,pch=21,bg=col)#-scores(mynmds)[,1],scores(mynmds)[,2]
ordiellipse(mynmds,groups=col,col=c("darkblue","dodgerblue","lightblue"),conf=.99999,kind="se",lwd=2)#




comm.dataBac[1:5,32:33]
ordi.bact[1:5,54:55]
