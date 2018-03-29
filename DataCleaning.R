#Euks: joined paired reads, not truncated (b/c short read and primer should be successfully removed b/c primer comes before poor quality bps), used SILVA all taxa nontruncated database
#Bact: joined paired reads, truncated and biased toward the smaller side to remove poor quality reads b/c read length is not very variable, used greengenes b/c there were fewer bacteria;__ and unassigned reads
#ITS: used forward read (b/c reverse read has poor quality early and joining reads removes lots of low quality reads b/c I can't truncate due to variable lengths and b/c forward read resuls in many more sequences retained and more features/taxa compared to joined reads, however, it also has more sequences classified as fungi;__), used UNITE database b/c SILVA resultd in nearly all sequnces unassigned


#######Read in OTU data#######

##### Read in euk files, not filtered #####

#first read in otu table and taxonomy file as txt, then merge (because the OTU table does not have taxonomy in it)

otufileEuk<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table2.txt",header=T)

head(otufileEuk)

#6 ranks
#taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy/taxonomy2.csv",header=T)

#all ranks, consensus, not truncated, all taxa
taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy5/taxonomy2.csv",header=T)

#taxonomyfileEuk[which(taxonomyfileEuk$Confidence>.7&taxonomyfileEuk$Confidence<.72),]
#taxonomyfileEuk[which(taxonomyfileEuk$Confidence<.7),"Taxon"] #all taxa with <.7 confidence are classified as unassigned 
#taxonomyfileEuk[4,]
#hist(taxonomyfileEuk$Confidence)

#“Confidence” is the fraction of top hits that match the consensus taxonomy (at whatever level is provided)

head(taxonomyfileEuk)
colnames(taxonomyfileEuk)[1]<-"OTUID"
colnames(taxonomyfileEuk)[2]<-"taxonomy"

otufileEuk2<-merge(otufileEuk,taxonomyfileEuk,"OTUID")
otufileEuk2<-otufileEuk2[,-dim(otufileEuk2)[2]] #delete the last column which is confidence in taxonomy
head(otufileEuk2)

dim(otufileEuk)
dim(otufileEuk2)
dim(taxonomyfileEuk)

#min(otufileEuk2$Confidence)

#Write it to a text file so that you can convert it to biom, which can then be import back into R using import_biom()
#6 ranks only euks
#write.table(otufileEuk2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table3.txt",sep="\t",row.names = F)

#all ranks and all taxa
#write.table(otufileEuk2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table4.txt",sep="\t",row.names = F)

#all ranks and all taxa, consensus, not truncated
#write.table(otufileEuk2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table5.txt",sep="\t",row.names = F)

#open otu_table3 in excel and add '#OTU ID' as first cell name
#biom convert -i otu_table3.txt -o otu_table3.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#biom convert -i otu_table4.txt -o otu_table4.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

biom convert -i otu_table5.txt -o otu_table5.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table3.biom",parseFunction = parse_taxonomy_greengenes)
#otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table5.biom",parseFunction = parse_taxonomy_greengenes)
otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table4.biom",parseFunction = parse_taxonomy_default)
#I'm not sure what the differences is between parse tax greensgenes and default, the warnigns toldme to use default

head(otu_table(otuEuk))
head(tax_table(otuEuk))


#Import mapping and tree file
mapEuk<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/EukBr_Niwot_20072015_All_MapFilenewlomehi.txt")

treeEuk<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-rooted-tree/tree.nwk")

datEuk<-merge_phyloseq(otuEuk,mapEuk,treeEuk)


##### Euk Soil #####
##Filter, filter only soil samples, 2015, take out bacteria, unassigned, plants, fungi, metazoa, take out singletons
#I think I should try blasting against the database with all levels not just 6. and I could try blasting against the entire database, not just eukaryotes
#in the previous bioinformatics pipeline all taxa had at least three levels (euk;class;"uncultured euk"). Many of the taxa here are only "eukaryote". So I might want to try downloading the older silva database and using that to see if I can replicate, or I can just delete all the "Euk" only assignments (there are rows that are labeled "unassigned" too that I should remove). I think the problem here was the BLAST assign taxonomy difference, not the database

#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement
datEukS<-datEuk%>%
  subset_samples(SampleType=="soil"&year==2015)%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="D_0__Bacteria")%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Unassigned")%>%
  #subset_taxa(is.na(Rank2)==F)%>% #this takes out the Eukaryota;__ taxa, I could or not take them out. they are at least real in the sense that they are not blasting to bacteria, but they could be artifacts (i.e. primer errors). I'll leave them in for now
  subset_taxa(is.na(Rank4)==T|Rank4!="D_3__Fungi")%>%
  subset_taxa(is.na(Rank4)==T|Rank4!="D_3__Metazoa(Animalia)")%>%
  subset_taxa(is.na(Rank7)==T|Rank7!="D_6__Embryophyta")%>%
  filter_taxa(function(x) sum(x) > (1), prune=T) #there are no singletons, odd, but I think that is what DADA2 does, it works with repeated samples to identify errors
datEukS

#using the euk only database, there is roughtly the same number of euk;__ taxa as from the all database euk;__+unassigned taxa
#4168/12475 are euk;__
#5789/12475
#1522 (unassigned)+4168 euk=5690

head(sample_data(datEukS))
head(tax_table(datEukS))
sort(unique(tax_table(datEukS)[,"Rank3"]))

length(which(is.na(tax_table(datEukS)[,"Rank2"])))
#1629 of the 4068 taxa are just assigned euk;__

#sample 61 has 681 reads, next sample76 has 946 reads, so could delete 61
min(sample_sums(datEukS))
sort(sample_sums(datEukS))
which(sample_data(datEukS)$X.SampleID=="S.81.2015")
#sample 81 is not included



###### Euk Nematode Samples #####
##Filter, filter only nematode samples, filter sample 2A (duplicate), only metazoa, take out craniata (chordata here, chordata is a higher level than craniata but the only chordata here are vertebrates), take out singletons (there are none)
datEukN<-datEuk%>%
  subset_samples(SampleType=="nematode"&X.SampleID!="N.2A.2015")%>%
  subset_taxa(Rank4=="D_3__Metazoa(Animalia)")%>%
  subset_taxa(is.na(Rank7)==T|Rank7!="D_6__Chordata")%>%
  filter_taxa(function(x) sum(x) > 1, prune=T) 
datEukN

head(sample_data(datEukS))
head(tax_table(datEukN))
sort(unique(tax_table(datEukN)[,"Rank7"]))

tax_table(datEukN)[which(tax_table(datEukN)[,"Rank7"]=="D_6__Chordata"),]

min(sample_sums(datEukS))
sort(sample_sums(datEukS))

#sample 78 has only 70 reads so delete; however this could be true b/c 78 is a zero plant plot so there are probably not a lot of nematodes, next sample 3 has 700 (3 is also a zero plant plot), next 83 has 2038 (also zero plants)
min(sample_sums(datEukN))
sort(sample_sums(datEukN))
which(sample_data(datEukN)$X.SampleID=="S.33.2015")
which(sample_data(datEukN)$X.SampleID=="S.56.2015")
#33 and 56 are missing

nemsforfunctionalgroup<-tax_table(datEukN)[which(tax_table(datEukN)[,"Rank7"]=="D_6__Nematoda"),]
unique(nemsforfunctionalgroup[,"Rank9"])
unique(nemsforfunctionalgroup[,"Rank11"])
nemsforfunctionalgroup[,"Rank9"]
unique(tax_table(datEukN)[,"Rank7"])






###### Bacteria #######
otufileBac<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-table/otu_table2.txt",header=T)

head(otufileBac)

#taxonomy from greengeenes
taxonomyfileBac<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-taxonomy_gg/taxonomy2.csv",header=T)

#taxonomyfileBac[which(taxonomyfileBac$Confidence<.7),"Taxon"] #all taxa have confidence of >.7
#hist(taxonomyfileBac$Confidence)

#“Confidence” is the fraction of top hits that match the consensus taxonomy (at whatever level is provided)

head(taxonomyfileBac)
colnames(taxonomyfileBac)[1]<-"OTUID"
colnames(taxonomyfileBac)[2]<-"taxonomy"

otufileBac2<-merge(otufileBac,taxonomyfileBac,"OTUID")
otufileBac2<-otufileBac2[,-dim(otufileBac2)[2]] #delete the last column which is confidence in taxonomy
head(otufileBac2)

dim(otufileBac)
dim(otufileBac2)
dim(taxonomyfileBac)

#min(otufileEuk2$Confidence)

#Write it to a text file so that you can convert it to biom, which can then be import back into R using import_biom()
write.table(otufileBac2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-table/otu_table3.txt",sep="\t",row.names = F)

#open otu_table3 in excel and add '#OTU ID' as first cell name
biom convert -i otu_table3.txt -o otu_table3.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

otuBac <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-table/otu_table3.biom",parseFunction = parse_taxonomy_default)
 
#I'm not sure what the differences is between parse tax greensgenes and default, the warnings told me to use default

head(otu_table(otuBac))
head(tax_table(otuBac))


#Import mapping and tree file
mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.txt")

treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-rooted-tree/tree.nwk")

datBac<-merge_phyloseq(otuBac,mapBac,treeBac)


##### Filter bacteria #####
##Filter, filter only soil samples, 2015, keep only bacteria and archaea (this filters unassigned), filter mitochondira and chloroplasts, take out singletons

#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement
datBacS<-datBac%>%
  subset_samples(SampleType=="soil"&year==2015)%>%
  subset_taxa(Rank1=="k__Archaea"|Rank1=="k__Bacteria")%>%
  #subset_taxa(is.na(Rank2)==F&Rank2!="p__")%>% #this takes out the Bacteria;__ taxa, I could or not take them out. they could be real or they could be artifacts (i.e. primer errors). I'll leave them in for now. there are 1652 taxa like this
  subset_taxa(is.na(Rank3)==T|Rank3!="c__Chloroplast")%>%
  subset_taxa(is.na(Rank5)==T|Rank5!="f__mitochondria")%>%
  filter_taxa(function(x) sum(x) > (1), prune=T) #there are no singletons, odd, but I think that is what DADA2 does, it works with repeated samples to identify errors
datBacS

head(sample_data(datBacS))
head(tax_table(datBacS))
sort(unique(tax_table(datBacS)[,"Rank3"]))

length(which(is.na(tax_table(datBacS)[,"Rank2"])))

#sample 5, 34, 126 have low # reads, 2, 46, 117 respectively
min(sample_sums(datBacS))
sort(sample_sums(datBacS))




###### ITS #######
otufileITS<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITSsingle/exported-table/otu_table2.txt",header=T)

head(otufileITS)

#taxonomy from unite
taxonomyfileITS<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITSsingle/exported-taxonomy_unite/taxonomy2.csv",header=T)

#taxonomyfileITS[which(taxonomyfileITS$Confidence<.7),"Taxon"] #all taxa except 5 unassigned have confidence of >.7
#hist(taxonomyfileITS$Confidence)

#“Confidence” is the fraction of top hits that match the consensus taxonomy (at whatever level is provided)

head(taxonomyfileITS)
colnames(taxonomyfileITS)[1]<-"OTUID"
colnames(taxonomyfileITS)[2]<-"taxonomy"

otufileITS2<-merge(otufileITS,taxonomyfileITS,"OTUID")
otufileITS2<-otufileITS2[,-dim(otufileITS2)[2]] #delete the last column which is confidence in taxonomy
head(otufileITS2)

dim(otufileITS)
dim(otufileITS2)
dim(taxonomyfileITS)

#sort(taxonomyfileITS$Confidence,decreasing=T)

#Write it to a text file so that you can convert it to biom, which can then be import back into R using import_biom()
write.table(otufileITS2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITSsingle/exported-table/otu_table3.txt",sep="\t",row.names = F)

#open otu_table3 in excel and add '#OTU ID' as first cell name
biom convert -i otu_table3.txt -o otu_table3.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

otuITS <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITSsingle/exported-table/otu_table3.biom",parseFunction = parse_taxonomy_default)

head(otu_table(otuITS))
head(tax_table(otuITS))


#Import mapping and tree file
mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/ITSsingle/ITS_Niwot_20072015_All_MapFilenewlomehi.txt")

datITS<-merge_phyloseq(otuITS,mapITS)


##### Filter ITS #####
##Filter, filter only soil samples, 2015, keep only Fungi (this filters unassigned and rhizaria), filter mitochondira and chloroplasts, take out singletons

#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement
datITSS<-datITS%>%
  subset_samples(SampleType=="soil"&year==2015)%>%
  subset_taxa(Rank1=="k__Fungi")%>%
  #subset_taxa(is.na(Rank2)==F)%>% #this takes out the Fungi;__ taxa, I could or not take them out. they could be real or they could be artifacts (i.e. primer errors). I'll leave them in for now. there are 5976 taxa like this - that's a lot!
  filter_taxa(function(x) sum(x) > (1), prune=T) #there are no singletons, odd, but I think that is what DADA2 does, it works with repeated samples to identify errors
datITSS

head(sample_data(datITSS))
head(tax_table(datITSS))
sort(unique(tax_table(datITSS)[,"Rank2"]))

length(which(is.na(tax_table(datITSS)[,"Rank2"])))

#sample 112 has kind of low # reads, 1301, the next lowest sample is 56 at 5117 reads but I think I'll keep it since I'm deleting too many samples from the other runs
min(sample_sums(datITSS))
sort(sample_sums(datITSS))
setdiff(sample_data(datBacS)$Sample_name,sample_data(datITSS)$Sample_name)
#ITS is missing 126, but we are already deleting that b/c low reads for bacteria


##### Remove samples across all datasets #####
#(from before) For euksS sample 81 did not amplify and 61 had low # reads. for euksN sample 33 and 56 did not have enough soil so were not done, and 78 had low # reads. for ITS 126 didnt amplify. for bacteria, samples 5,34,126 did not amplify. should have 90 samples left

#From EukS 61, 81
#From EukN 78, 33, 56
#From Bac 5, 34, 126
#From ITS 126

datEukS2<-datEukS%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)
sum(otu_table(datEukS2))

datEukN2<-datEukN%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

datBacS2<-datBacS%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

datITSS2<-datITSS%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)



##### Rarefy and transform to relative abundance #####
min(sample_sums(datEukS2))#rarefy to 946
min(sample_sums(datEukN2))#rarefy to 700
min(sample_sums(datBacS2))#rarefy to 5868
min(sample_sums(datITSS2))#rarefy to 1301

datEukS3<-datEukS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukS2)),rngseed=10,replace=F) %>%
  transform_sample_counts(function(x) x/sum(x) )
#609 OTUs were removed because they are no longer present in any sample after random subsampling, 

datEukN3<-datEukN2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukN2)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
  #91 OTUs were removed because they are no longer present in any sample after random subsampling, 

datBacS3<-datBacS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBacS2)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#1159 OTUs were removed because they are no longer present in any sample after random subsampling, 

datITSS3<-datITSS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS2)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#3276 OTUs were removed because they are no longer present in any sample after random subsampling, 


#rarefying but not calculating relative abundance
datEukS3c<-datEukS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukS2)),rngseed=10,replace=F)
#609 OTUs were removed because they are no longer present in any sample after random subsampling, 

datEukN3c<-datEukN2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukN2)),rngseed=10,replace=F)
#91 OTUs were removed because they are no longer present in any sample after random subsampling, 

datBacS3c<-datBacS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBacS2)),rngseed=10,replace=F)
#1159 OTUs were removed because they are no longer present in any sample after random subsampling, 

datITSS3c<-datITSS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS2)),rngseed=10,replace=F)
#3276 OTUs were removed because they are no longer present in any sample after random subsampling, 






##### Make a label column with kingdom/phylum level labels #####

#Making groups for labeling graphs, I will just use "euks" "mesofauna" "fungi" "bacteria" for now, I can come back and label more precisely if needed

##### Euks Soil #####

labelsEukS<-data.frame(rep("Eukaryota",dim(tax_table(datEukS3))[1]))

# unique(tax_table(datEukS3)[,"Rank3"])
# tax_table(datEukS3)[which(tax_table(datEukS3)[,"Rank2"]=="D_1__Archaeplastida"),]
# labelsEukS<-tax_table(datEukS3)[,"Rank3"]#
# 
# ind<-which(is.na(tax_table(datEukS3)[,"Rank3"]))
# labelsEukS[ind]<-tax_table(datEukS3)[ind,"Rank2"]
# 
# ind<-which(is.na(labelsEukS))
# labelsEukS[ind]<-tax_table(datEukS3)[ind,"Rank1"]
# 
# ind<-which(labelsEukS=="D_2__IncertaeSedis")
# labelsEukS[ind]<-tax_table(datEukS3)[ind,"Rank2"]

colnames(labelsEukS)<-"labels"

#replace tax table
#tax_table(datEukS3)<-cbind(tax_table(datEukS3),labelsEukS)


##### Euks N #####

labelsEukN<-data.frame(rep("Mesofauna",dim(tax_table(datEukN3))[1]))

# unique(tax_table(datEukN3)[,"Rank7"])
# labelsEukN<-tax_table(datEukN3)[,"Rank7"]#
# 
# ind<-which(is.na(labelsEukN))
# labelsEukN[ind]<-tax_table(datEukN3)[ind,"Rank6"]
# 
# ind<-which(is.na(labelsEukN))
# labelsEukN[ind]<-tax_table(datEukN3)[ind,"Rank5"]

colnames(labelsEukN)<-"labels"

#replace tax table
#tax_table(labelsEukN)<-cbind(tax_table(labelsEukN),labelsEukN)


##### Bacteria #####

labelsBac<-data.frame(rep("Bacteria",dim(tax_table(datBacS3))[1]))
colnames(labelsBac)<-"labels"

##### Fungi #####

labelsITS<-data.frame(rep("Fungi",dim(tax_table(datITSS3))[1]))
colnames(labelsITS)<-"labels"


##### Make otu tables (bact takes ~10 min) #####
datEukS3otu<-cbind(sample_data(datEukS3),t(otu_table(datEukS3)))
datEukS3otu$Sample_name<-as.numeric(as.character(datEukS3otu$Sample_name))

datEukN3otu<-cbind(sample_data(datEukN3),t(otu_table(datEukN3)))
datEukN3otu$Sample_name<-as.numeric(as.character(datEukN3otu$Sample_name))

datBacS3otu<-cbind(sample_data(datBacS3),t(otu_table(datBacS3)))
datBacS3otu$Sample_name<-as.numeric(as.character(datBacS3otu$Sample_name))

datITSS3otu<-cbind(sample_data(datITSS3),t(otu_table(datITSS3)))
datITSS3otu$Sample_name<-as.numeric(as.character(datITSS3otu$Sample_name))

#otu for counts (not relative abundance)
datEukS3cotu<-cbind(sample_data(datEukS3c),t(otu_table(datEukS3c)))
datEukS3cotu$Sample_name<-as.numeric(as.character(datEukS3cotu$Sample_name))

datEukN3cotu<-cbind(sample_data(datEukN3c),t(otu_table(datEukN3c)))
datEukN3cotu$Sample_name<-as.numeric(as.character(datEukN3cotu$Sample_name))

datBacS3cotu<-cbind(sample_data(datBacS3c),t(otu_table(datBacS3c)))
datBacS3cotu$Sample_name<-as.numeric(as.character(datBacS3cotu$Sample_name))

datITSS3cotu<-cbind(sample_data(datITSS3c),t(otu_table(datITSS3c)))
datITSS3cotu$Sample_name<-as.numeric(as.character(datITSS3cotu$Sample_name))

datITSS3cotu[1:10,20:50]


###### Filter data sets for network analysis ######
#Follwing Widder et al 2014 PNAS

#take out doubletons and singletons
datEukS3otu2<-cbind(datEukS3otu[,1:31],datEukS3otu[,((which(colSums(datEukS3otu[,32:dim(datEukS3otu)[2]]>0)>2))+31)])

datEukN3otu2<-cbind(datEukN3otu[,1:31],datEukN3otu[,((which(colSums(datEukN3otu[,32:dim(datEukN3otu)[2]]>0)>2))+31)])

datBacS3otu2<-cbind(datBacS3otu[,1:31],datBacS3otu[,((which(colSums(datBacS3otu[,32:dim(datBacS3otu)[2]]>0)>2))+31)])

datITSS3otu2<-cbind(datITSS3otu[,1:31],datITSS3otu[,((which(colSums(datITSS3otu[,32:dim(datITSS3otu)[2]]>0)>2))+31)])

#filter out taxa that have a summed relative abundance of <.002 (.2%)
datEukS3otu3<-cbind(datEukS3otu2[,1:31],datEukS3otu2[,((which(colSums(datEukS3otu2[,32:dim(datEukS3otu2)[2]])>0.002))+31)]) #1124 otu

datEukN3otu3<-cbind(datEukN3otu2[,1:31],datEukN3otu2[,((which(colSums(datEukN3otu2[,32:dim(datEukN3otu2)[2]])>0.002))+31)]) #143 otu

datBacS3otu3<-cbind(datBacS3otu2[,1:31],datBacS3otu2[,((which(colSums(datBacS3otu2[,32:dim(datBacS3otu2)[2]])>0.002))+31)]) #3399 otu

datITSS3otu3<-cbind(datITSS3otu2[,1:31],datITSS3otu2[,((which(colSums(datITSS3otu2[,32:dim(datITSS3otu2)[2]])>0.002))+31)]) #1122 otu


#order of doing things: filtered out unwanted taxa (chloroplasts, spiders) and samples (S.2015), rarefied, relativized, took out doubletons and singletons, took out samples <2% summed abundance. I am not going to rarefy or relativize here again b/c the doubletons/singletons/.2% otus that I removed are real, I could have included them in the network analyis if I wanted.


##### Plants #####

plantcomp<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Plants/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots; there are some plots that have plant data but not microbe data. 69 70 71 77 81 108 117 118 147 148 149 151. This is because when we started doing the surveys we were going to all plots for plants and only sample some for microbes, then we realized that that was insane!
dim(plantcomp2)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantcomp2$LICHEN<-NULL
labelsPlant<-as.data.frame(cbind(labels="Plant",otu=colnames(plantcomp2)[2:55],oldotu=colnames(plantcomp2)[2:55]))
head(labelsPlant)



###### Make labelfile ######

#Test if there is name overlap in full otu tables, no names overlap, no need to rename OTU columns, however it is nice for figuring out what is going on in the modeling, so I will add B, N, S, I before all the otu names
namesEukS<-names(datEukS3otu[,-c(1:31)])
namesEukN<-names(datEukN3otu[,-c(1:31)])
namesBac<-names(datBacS3otu[,-c(1:31)])
namesITS<-names(datITSS3otu[,-c(1:31)])
length(namesBac)+length(namesEukS)
length(union(namesBac,namesEukS))
intersect(namesEukS,namesEukN)
intersect(namesEukS,namesBac)
intersect(namesEukS,namesITS)
intersect(namesEukN,namesBac)
intersect(namesEukN,namesITS)
intersect(namesBac,namesITS)

namesEukS2 <- sub("^", "S", namesEukS)
namesEukN2 <- sub("^", "N", namesEukN)
namesBac2 <- sub("^", "B", namesBac)
namesITS2 <- sub("^", "I", namesITS)

names(datEukS3otu)[-c(1:31)]<-namesEukS2
names(datEukN3otu)[-c(1:31)]<-namesEukN2
names(datBacS3otu)[-c(1:31)]<-namesBac2
names(datITSS3otu)[-c(1:31)]<-namesITS2


#Change colnames on the reduced datasets for the networks
namesEukS<-names(datEukS3otu3[,-c(1:31)])
namesEukN<-names(datEukN3otu3[,-c(1:31)])
namesBac<-names(datBacS3otu3[,-c(1:31)])
namesITS<-names(datITSS3otu3[,-c(1:31)])

namesEukS2 <- sub("^", "S", namesEukS)
namesEukN2 <- sub("^", "N", namesEukN)
namesBac2 <- sub("^", "B", namesBac)
namesITS2 <- sub("^", "I", namesITS)

names(datEukS3otu3)[-c(1:31)]<-namesEukS2
names(datEukN3otu3)[-c(1:31)]<-namesEukN2
names(datBacS3otu3)[-c(1:31)]<-namesBac2
names(datITSS3otu3)[-c(1:31)]<-namesITS2


#Set rownames on the label file
labelsEukS2<-labelsEukS
rownames(labelsEukS2)<-sub("^", "S",rownames(labelsEukS2))
labelsEukN2<-labelsEukN
rownames(labelsEukN2)<-sub("^", "N",rownames(labelsEukN2))
labelsBac2<-labelsBac
rownames(labelsBac2)<-sub("^", "B",rownames(labelsBac2))
labelsITS2<-labelsITS
rownames(labelsITS2)<-sub("^", "I",rownames(labelsITS2))

# rownames(labelsEukS)<-namesEukS
# rownames(labelsEukN)<-namesEukN
# rownames(labelsBac)<-namesBac
# rownames(labelsITS)<-namesITS

labelfile1<-rbind(labelsEukS2,labelsEukN2,labelsBac2,labelsITS2)
head(labelfile1)
labelfile1$otu<-rownames(labelfile1)
labelfile1$oldotu<-substring(rownames(labelfile1), 2)
labelfile<-rbind(labelfile1,labelsPlant)
tail(labelfile)


#Correct rownames for counts data that was not relativized, note I'm overwriting the file names above
namesEukS<-names(datEukS3cotu[,-c(1:31)])
namesEukN<-names(datEukN3cotu[,-c(1:31)])
namesBac<-names(datBacS3cotu[,-c(1:31)])
namesITS<-names(datITSS3cotu[,-c(1:31)])

namesEukS2 <- sub("^", "S", namesEukS)
namesEukN2 <- sub("^", "N", namesEukN)
namesBac2 <- sub("^", "B", namesBac)
namesITS2 <- sub("^", "I", namesITS)

names(datEukS3cotu)[-c(1:31)]<-namesEukS2
names(datEukN3cotu)[-c(1:31)]<-namesEukN2
names(datBacS3cotu)[-c(1:31)]<-namesBac2
names(datITSS3cotu)[-c(1:31)]<-namesITS2







##### Phylogenetic diversity ######
#not done yet
pdBac<-pd(as.matrix(datBacr3fotu[,-c(1:31)]),phy_tree(datBac3f),include.root=TRUE) #takes a few hours, started 3:45pm, finished ~midnight
pdEukN<-pd(as.matrix(datEukN5fotu[,-c(1:31)]),phy_tree(datEukN5f),include.root=TRUE) #took 3 minutes
pdEukS<-pd(as.matrix(datEukS4fotu[,-c(1:31)]),phy_tree(datEukS4f),include.root=TRUE) #took 5 minutes
richITS<-as.data.frame(rowSums(datITS3fotu[,-c(1:31)]>0))
pdMet<-pd(as.matrix(dat99Met2fotu[,-c(1:31)]),phy_tree(dat99Met2f),include.root=TRUE) #run code from metazoa section below first
pdNem<-pd(as.matrix(dat99Met2fNemotu[,-c(1:31)]),phy_tree(dat99Met2fNem),include.root=TRUE) # I need to use root=T, if I use root=F it cannot calculate the diversity in a sample with only one taxon












###### old code #####

#Amoebozoa (kingdom)
#photosynthetic Alveolata (kingdom) (phylum Dinoflagellata: mostly photosynthetic; nonphotosynthetic Protoperidinium, both SL163A10 and Pfiesteria (can if it eats an alga), unknown D244), 
#nonphotosynthetic Alveolata (phyla Ciliophora(predators), protalveolata, apicomplexa (parasites)), BOLA553 is a near apicomplexan, not sure what SCM37C52 is so I put it here
#Archaeplastida (major_clade), combine the kingdoms Chloroplastida, Rhodophyceae, Haptophyta, 
#Rhizaria (kingdom: unicellular amoeboid euks), 
#Holozoa(kingdom: all animals, not fungi), includes metazoa (animals)
#photosynthetic Excavata major clade (was Discoba kingdom) (class Euglenida: mostly photosynthetic), 
#nonphotosnthetic Excavata (was Discoba) (in discoba, phylum Heterolobosea: parasites, free living, symbiotic, amoeba-like, and phylum Jakoba heterotrophic), superphylum metamonada, doesn't have kingdom, parasites not Ps.
#photosynthetic Cryptophyceae is a kingdom, they have one or two chloroplasts, except for Chilomonas, which has leucoplasts (not pigmented not for Ps) and Goniomonas (formerly Cyathomonas) which lacks plastids entirely.P1-31 (i can't find much info on this, some sites called them picoplankton and most of the cryptophyceae are photosynthetic). Cryptomonadales can't find much about it but assume Ps. __SA1-3C06 is at rank2 but notes in silva say "cryptophyte"
#nonphotosynthetic cryptophyceae - Chilomonas and Goniomonas (formerly Cyathomonas)
#Haptophyta kingdom, I think all are Ps
#Centrohelida - encyclopedia of life said some species have photosynthetic symbionts but I can't find any info on Ps taxa.Heterophryidae, Mb-5C (can't find any info on, I'll put it with nonPs), M1-18D08(can't find info),Acanthocystidae(can't find info),Pterocystis(can't find info), so I will leave them all in a group called "centrohelida"
#NonPs Stramenopiles: MAST-3, MAST-12, MAST-7, Labyrinthulomycetes,Bicosoecida,Peronosporomycetes,Hyphochytriales,__Incertae_Sedis(were all __Pirsonia which is a parasite of algae)
#Ps Stramenopiles: Pelagophycea, Diatomea, Eustigmatales, Xanthophyceae, Phaeothamniophyceae, Chrysophyceae,Ochromonas,CCI40 (not much info, some places call it is chrysophyte),Synurales, Raphidophyceae, NA1-2A5 (can't find any info so putting it with Ps), TKR07M.92

#heterotrophic eukaryota (things without a kingdom): Picozoa (phylum nonphotosynthetic unicellular eukaryotes, only 1 OTU); RT5iin14; Apusomonadidae (nonphotosynthetic protist, unknown kingdom); Kathablepharidae (nonPs protist, unkonwn kingdo);Telonema (nonPs protist genus);Breviatea (amoeba like);DH147-EKD23;__Ancyromonadida;GoC2-B10 (no info);Zeuk77 (no info)

labelsEukS<-tax_table(datEukS3)[,"Rank3"]#

ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Amoebozoa")
labelsEukS[ind]<-"Amoebozoa"
ind<-which(tax_table(datEukS4)[,"Rank4"]=="__Dinoflagellata")#__Protoperidinium is included but I then change it below becuase some things have an NA for rank8
labelsEukS[ind]<-"Photosynthetic_Alveolata"
ind<-which(tax_table(datEukS4)[,"Rank8"]=="__Protoperidinium")
labelsEukS[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(tax_table(datEukS4)[,"Rank4"]!="__Dinoflagellata"&tax_table(datEukS4)[,"Rank3"]=="__Alveolata")
labelsEukS[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Archaeplastida")
labelsEukS[ind]<-"Archaeplastida"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Rhizaria")
labelsEukS[ind]<-"Rhizaria"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Holozoa")
labelsEukS[ind]<-"Holozoa"
ind<-which(tax_table(datEukS4)[,"Rank6"]=="__Euglenida")
labelsEukS[ind]<-"Photosynthetic_Excavata"
ind<-which(tax_table(datEukS4)[,"Rank6"]!="__Euglenida"&tax_table(datEukS4)[,"Rank2"]=="__Excavata")
labelsEukS[ind]<-"Nonphotosynthetic_Excavata"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Goniomonas")
labelsEukS[ind]<-"Nonphotosynthetic_Chryptophyceae"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__P1-31"|tax_table(datEukS4)[,"Rank3"]=="__Cryptomonadales"|tax_table(datEukS4)[,"Rank2"]=="__SA1-3C06")
labelsEukS[ind]<-"Photosynthetic_Chryptophyceae"
ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Haptophyta")
labelsEukS[ind]<-"Haptophyta"
ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Centrohelida"&tax_table(datEukS4)[,"Rank3"]!="__uncultured_Sarcosomataceae")
labelsEukS[ind]<-"Centrohelida"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__uncultured_Sarcosomataceae")#i'm putting the wierd centrohelid/fungus with centrohilid (since I supposedly took out al fungi)
labelsEukS[ind]<-"Centrohelida"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS4)[,"Rank4"]%in%c("__MAST-3","__MAST-12","__MAST-7","__Labyrinthulomycetes","__Bicosoecida","__Peronosporomycetes","__Hyphochytriales","__Incertae_Sedis"))
labelsEukS[ind]<-"Nonphotosynthetic_Stramenopiles"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS4)[,"Rank4"]%in%c("__Pelagophycea","__Diatomea","__Eustigmatales","__Xanthophyceae","__Phaeothamniophyceae","__Chrysophyceae","__Ochromonas","__CCI40","__Synurales","__Raphidophyceae","__NA1-2A5","__TKR07M.92","__Pelagophyceae"))
labelsEukS[ind]<-"Photosynthetic_Stramenopiles"

ind<-which(tax_table(datEukS4)[,"Rank3"]=="__RT5iin14"|tax_table(datEukS4)[,"Rank2"]=="__Picozoa"|tax_table(datEukS4)[,"Rank2"]=="__RT5iin25"|tax_table(datEukS4)[,"Rank3"]=="__Apusomonadidae"|tax_table(datEukS4)[,"Rank2"]=="__Kathablepharidae"|tax_table(datEukS4)[,"Rank3"]=="__Telonema"|tax_table(datEukS4)[,"Rank3"]=="__Breviatea"|tax_table(datEukS4)[,"Rank2"]=="__DH147-EKD23"|tax_table(datEukS4)[,"Rank3"]=="__Ancyromonadida"|tax_table(datEukS4)[,"Rank2"]=="__GoC2-B10"|tax_table(datEukS4)[,"Rank2"]=="__Zeuk77")
labelsEukS[ind]<-"Heterotrophic_Eukarya"


#I could separate excavata into kingdoms/super phyla and also separate archaeplastida into its 3 kingdoms
unique(labelsEukS)
colnames(labelsEukS)<-"labels"

#replace tax table for datEukS4 (relative abun) and datEukS3 (not relativized)
tax_table(datEukS4)<-cbind(tax_table(datEukS4),labelsEukS)
tax_table(datEukS3)<-cbind(tax_table(datEukS3),labelsEukS)

#looking where the weird centrohelid fungi comes out in a tree
#ex1<-subset_taxa(datEukr2,Rank3=="__Fungi")
#myTaxa<-names(taxa_sums(ex1))#[1:46]
# myTaxa<-c(myTaxa,"denovo65528")
# ex3<-subset_taxa(datEukr2,Rank2=="__Centrohelida")
# myTaxa<-c(myTaxa,names(taxa_sums(ex3)))
# ex3<-subset_taxa(datEukr2,Rank3=="__Rhizaria")
# myTaxa<-c(myTaxa,names(taxa_sums(ex3)))
# ex2 = prune_taxa(myTaxa, datEukr2)
# plot_tree(ex2, label.tips = "taxa_names",color="Rank3")

#look at tree of abundant taxa
#myTaxa<-c(names(sort(taxa_sums(datEukr2),decreasing=T)))[1:100]
#ex2 = prune_taxa(myTaxa, datEukr2)
#plot_tree(ex2, label.tips = "taxa_names",color="Rank3")







#remove the p__ with substring
labelsITS<-substring(tax_table(datITS3)[,"Rank2"],4)

#remove the __ with substring
labelsBac<-substring(tax_table(datBac3)[,"Rank2"],3)


#Bacteria
#make the class the label for anything in the Phylum Chloroflexi, Only the class Chloroflexi are photosynthetic, the other classes (Ktedontobacteria) are not photosynthetic. There many OTU listed as phylum Chloroflexi but that do not have a Class. The majority of the taxa in Chloroflexi that do have a class are Ktedonobacteria. Thus I will make the default color for heterotrophs (Purple), and just change the Chloroflexales to photosynthetic (green)
ind<-which(labelsBac=="Chloroflexi")
temp<-data.frame(tax_table(datBac3))
length(which(temp$Rank2=="__Chloroflexi"))
length(which(temp$Rank2=="__Chloroflexi"&temp$Rank3=="__Chloroflexi"))#none are known to be in class chloroflexi
length(which(temp$Rank3=="__Chloroflexales"|temp$Rank4=="__Chloroflexaceae")) #81 are Chloroflexales which I imagine are in the class Chloroflexi
length(which(temp$Rank2=="__Chloroflexi"&temp$Rank4=="__Chloroflexaceae")) #61 are Chloroflexaceae, these are included in the 81 above
length(which(temp$Rank2=="__Chloroflexi"&temp$Rank3=="__Ktedonobacteria"))#1497 are Ktedonobacteria
temp<-subset(temp,Rank2=="__Chloroflexi")
sort(unique(temp$Rank3))
sort(unique(temp$Rank4))
#to label the chloroflexales as phototrophs
temp<-data.frame(tax_table(datBac3))
ind<-which(temp$Rank3=="__Chloroflexales"|temp$Rank4=="__Chloroflexaceae")
labelsBacchloro<-labelsBac
labelsBacchloro[ind]<-"Chloroflexiphoto"
#to label the ktedonobacteria as heterotrophs
# temp<-data.frame(tax_table(datBac3))
# ind<-which(temp$Rank3=="__Ktedonobacteria")
# labelsBacchloro<-labelsBac
# labelsBacchloro[ind]<-"Chloroflexihetero"

colnames(labelsBac)<-"labels"
colnames(labelsBacchloro)<-"labels"






#notes: it looks like both prune_samples and subset_samples do NOT remove taxa that have an abundance of 0 after the samples are removed.




###### Grouping by kingdom/phylum #####

datBac3fk<-aggregate.data.frame(otu_table(datBac3f),by=list(labels=tax_table(datBac3f)[,"labels"]),sum)
rownames(datBac3fk)<-datBac3fk$labels
datBac3fk$labels<-NULL
datBac3fk2<-cbind(sample_data(datBac3f),t(datBac3fk))
head(datBac3fk2)














