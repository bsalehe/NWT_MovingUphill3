#Euks: joined paired reads, not truncated (b/c short read and primer should be successfully removed b/c primer comes before poor quality bps), used SILVA all taxa nontruncated database
#Bact: used forward read because it yeilded 1.6 times as many reads as when I tried joining paired reads. truncated and biased toward the smaller side to remove poor quality reads b/c read length is not very variable, used greengenes b/c there were fewer bacteria;__ and unassigned reads
#ITS: used forward read (b/c reverse read has poor quality early and joining reads removes lots of low quality reads b/c I can't truncate due to variable lengths and b/c forward read resuls in many more sequences retained and more features/taxa compared to joined reads, however, it also has more sequences classified as fungi;__), used UNITE database b/c SILVA resultd in nearly all sequnces unassigned


#######Read in OTU data#######

##### Read in euk files, not filtered #####

#first read in otu table and taxonomy file as txt, then merge (because the OTU table does not have taxonomy in it)

otufileEuk<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table2.txt",header=T)

head(otufileEuk)

#6 ranks
#taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy/taxonomy2.csv",header=T)

#all ranks, consensus, truncated, all taxa, 128 release. I added this after the fact because I think I had to have imported it to create otu_table4 below
#taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy2/taxonomy2.csv",header=T)

#all ranks, consensus, not truncated, all taxa, 128 release
#taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy5/taxonomy2.csv",header=T)

#all ranks, consensus, not truncated, all taxa, 111 release
taxonomyfileEuk<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-taxonomy6/taxonomy2.csv",header=T)

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

#all ranks and all taxa, consensus, not truncated, release 128
#write.table(otufileEuk2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table5.txt",sep="\t",row.names = F)

#all ranks and all taxa, consensus, not truncated, release 111
#write.table(otufileEuk2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table6.txt",sep="\t",row.names = F)

#open otu_table3 in excel and add '#OTU ID' as first cell name
#biom convert -i otu_table3.txt -o otu_table3.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
#biom convert -i otu_table4.txt -o otu_table4.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
#biom convert -i otu_table5.txt -o otu_table5.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
biom convert -i otu_table6.txt -o otu_table6.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table3.biom",parseFunction = parse_taxonomy_greengenes)
#otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table5.biom",parseFunction = parse_taxonomy_greengenes)
#otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table4.biom",parseFunction = parse_taxonomy_default) #when I switched to release 111, this was the file that was not commented out. The input is from exported_taxonomy2/ (all taxa database that was truncated to 150bp). I confirmed that otu_table4 was the one I worked with in all the preliminary modeling. NOT otu_table5. however they are very very similar in their taxonomy assigments. except for example for 2e58f27fc29d069b48193979058524ab which is classified as Eukaryota in otu_table5 and Eukaryota;SAR in otu_table4
otuEuk <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-table/otu_table6.biom",parseFunction = parse_taxonomy_default)
#I'm not sure what the differences is between parse tax greensgenes and default, the warnings told me to use default

head(otu_table(otuEuk))
head(tax_table(otuEuk))
tax_table(otuEuk)[2245,] #the row for taxa 2e58f27fc29d069b48193979058524ab

#Import mapping and tree file
mapEuk<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/EukBr_Niwot_20072015_All_MapFilenewlomehi.txt")

treeEuk<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Euks/exported-rooted-tree/tree.nwk")

datEuk<-merge_phyloseq(otuEuk,mapEuk,treeEuk)


##### Euk Soil for 128 release#####
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



###### Euk Nematode Samples 128 release #####
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




##### Euk Soil for 111 release#####
##Filter, filter only soil samples, 2015, take out bacteria, unassigned, archaea, plants, fungi, metazoa, take out singletons
#I think I should try blasting against the database with all levels not just 6. and I could try blasting against the entire database, not just eukaryotes
#in the previous bioinformatics pipeline all taxa had at least three levels (euk;class;"uncultured euk"). Many of the taxa here are only "eukaryote". So I might want to try downloading the older silva database and using that to see if I can replicate, or I can just delete all the "Euk" only assignments (there are rows that are labeled "unassigned" too that I should remove). I think the problem here was the BLAST assign taxonomy difference, not the database

#note!!!! if you filter with subset_taxa and a !=, it will NOT return any rows that are NA, so you always have to do an "or NA" in the statement
datEukS<-datEuk%>%
  subset_samples(SampleType=="soil"&year==2015)%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Bacteria")%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Unassigned")%>%
  subset_taxa(is.na(Rank1)==T|Rank1!="Archaea")%>%
  #subset_taxa(is.na(Rank2)==F)%>% #this takes out the Eukaryota;__ taxa, I could or not take them out. they are at least real in the sense that they are not blasting to bacteria, but they could be artifacts (i.e. primer errors). I'll leave them in for now
  subset_taxa(is.na(Rank3)==T|Rank3!="__Fungi")%>%
  subset_taxa(is.na(Rank3)==T|Rank3!="__Metazoa")%>%
  subset_taxa(is.na(Rank7)==T|Rank7!="__Embryophyta")%>%
  filter_taxa(function(x) sum(x) > (1), prune=T) #there are no singletons, odd, but I think that is what DADA2 does, it works with repeated samples to identify errors
datEukS

head(sample_data(datEukS))
head(tax_table(datEukS))
sort(unique(tax_table(datEukS)[,"Rank7"]))

length(which(is.na(tax_table(datEukS)[,"Rank2"])))
#1680 of the 3932 taxa are just assigned euk;__

#sample 61 has 680 reads, next sample76 has 871 reads, so could delete 61
min(sample_sums(datEukS))
sort(sample_sums(datEukS))
which(sample_data(datEukS)$X.SampleID=="S.81.2015")
#sample 81 is not included



###### Euk Nematode Samples 111 release #####
##Filter, filter only nematode samples, filter sample 2A (duplicate), only metazoa, take out craniata (chordata here, chordata is a higher level than craniata but the only chordata here are vertebrates), take out singletons (there are none)
datEukN<-datEuk%>%
  subset_samples(SampleType=="nematode"&X.SampleID!="N.2A.2015")%>%
  subset_taxa(Rank3=="__Metazoa")%>%
  subset_taxa(is.na(Rank4)==T|Rank4!="__Craniata")%>%
  filter_taxa(function(x) sum(x) > 1, prune=T) 
datEukN

head(sample_data(datEukS))
head(tax_table(datEukN))
sort(unique(tax_table(datEukN)[,"Rank4"]))

tax_table(datEukN)[which(tax_table(datEukN)[,"Rank4"]=="__Craniata"),]

#sample 78 has only 70 reads so delete; however this could be true b/c 78 is a zero plant plot so there are probably not a lot of nematodes, next sample 3 has 700 (3 is also a zero plant plot), next 83 has 2038 (also zero plants)
min(sample_sums(datEukN))
sort(sample_sums(datEukN))
which(sample_data(datEukN)$X.SampleID=="S.33.2015")
which(sample_data(datEukN)$X.SampleID=="S.56.2015")
#33 and 56 are missing

#nemsforfunctionalgroup<-tax_table(datEukN)[which(tax_table(datEukN)[,"Rank4"]=="__Nematoda"),]
#unique(nemsforfunctionalgroup[,"Rank6"])
#unique(nemsforfunctionalgroup[,"Rank8"])
#nemsforfunctionalgroup[,"Rank9"]
#unique(tax_table(datEukN)[,"Rank7"])








###### Bacteria #######
#paired ends
otufileBac<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-table/otu_table2.txt",header=T)

#single (forward) read
otufileBac<-read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/exported-table/otu_table2.txt",header=T)

head(otufileBac)

#taxonomy from greengeenes
#paired ends
taxonomyfileBac<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-taxonomy_gg/taxonomy2.csv",header=T)
#forward read
taxonomyfileBac<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/exported-taxonomy_gg/taxonomy2.csv",header=T)

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
#paired end 
write.table(otufileBac2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-table/otu_table3.txt",sep="\t",row.names = F)
#forward read
write.table(otufileBac2,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/exported-table/otu_table3.txt",sep="\t",row.names = F)

#open otu_table3 in excel and add '#OTU ID' as first cell name
biom convert -i otu_table3.txt -o otu_table3.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

#paired
otuBac <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-table/otu_table3.biom",parseFunction = parse_taxonomy_default)
#forward
otuBac <- import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/exported-table/otu_table3.biom",parseFunction = parse_taxonomy_default)

#I'm not sure what the differences is between parse tax greensgenes and default, the warnings told me to use default

head(otu_table(otuBac))
head(tax_table(otuBac))


#Import mapping and tree file
mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.txt")

#paired
treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bact/exported-rooted-tree/tree.nwk")
#single
treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/QIIME2/Bactsingle/exported-rooted-tree/tree.nwk")

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

#sample 5, 34, 126 have low # reads, 2, 46, 17 respectively
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
min(sample_sums(datEukS2))#rarefy to 871 (was 946 with release 128)
min(sample_sums(datEukN2))#rarefy to 700
min(sample_sums(datBacS2))#rarefy to 10825
min(sample_sums(datITSS2))#rarefy to 1301

datEukS3<-datEukS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukS2)),rngseed=10,replace=F) %>%
  transform_sample_counts(function(x) x/sum(x) )
#609 OTUs were removed because they are no longer present in any sample after random subsampling, 

datEukN3<-datEukN2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datEukN2)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#81 OTUs were removed because they are no longer present in any sample after random subsampling, 

datBacS3<-datBacS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBacS2)),rngseed=10,replace=F)%>%
  transform_sample_counts(function(x) x/sum(x) )
#1983 OTUs were removed because they are no longer present in any sample after random subsampling, 
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
#81 OTUs were removed because they are no longer present in any sample after random subsampling, 

datBacS3c<-datBacS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datBacS2)),rngseed=10,replace=F)
#1983 OTUs were removed because they are no longer present in any sample after random subsampling, 
datITSS3c<-datITSS2%>%
  rarefy_even_depth(sample.size=min(sample_sums(datITSS2)),rngseed=10,replace=F)
#3276 OTUs were removed because they are no longer present in any sample after random subsampling, 






##### Make a label column with kingdom/phylum level labels #####

##### Euks Soil #####

#labelsEukS<-data.frame(rep("Eukaryota",dim(tax_table(datEukS3))[1]))

#Amoebozoa (kingdom)
#Alveolata - right now I'm not separating out Ps vs non-Ps b/c there are only 4 dinoflagallata and I can't tell whether they are Ps or not. there are a lot of things classified only to "Alveolata" so I can't tell what function they are.
  #photosynthetic Alveolata (kingdom) (phylum Dinoflagellata: mostly photosynthetic; nonphotosynthetic Protoperidinium, both SL163A10 and Pfiesteria (can if it eats an alga), unknown D244),) only 4 taxa, not sure if they are Ps but I'll leave them here. 
  #nonphotosynthetic Alveolata (phyla Ciliophora(predators), protalveolata, apicomplexa (parasites)), BOLA553 is a near apicomplexan, not sure what SCM37C52 is so I put it here
#Archaeplastida (major_clade), combine the kingdoms Chloroplastida, Rhodophyceae, Haptophyta, 
#Rhizaria (kingdom: unicellular amoeboid euks), 
#Holozoa(kingdom: all animals, not fungi), includes metazoa (animals)
#photosynthetic Excavata major clade (was Discoba kingdom) (class Euglenida: mostly photosynthetic), 
#nonphotosnthetic Excavata (was Discoba) (in discoba, phylum Heterolobosea: parasites, free living, symbiotic, amoeba-like, and phylum Jakoba heterotrophic), superphylum metamonada, doesn't have kingdom, parasites not Ps. Kinetoplastea: not Ps, Tetramitia: most likely not Ps
#we dont have any cryptophyceae
  #photosynthetic Cryptophyceae is a kingdom, they have one or two chloroplasts, except for Chilomonas, which has leucoplasts (not pigmented not for Ps) and Goniomonas (formerly Cyathomonas) which lacks plastids entirely.P1-31 (i can't find much info on this, some sites called them picoplankton and most of the cryptophyceae are photosynthetic). Cryptomonadales can't find much about it but assume Ps. __SA1-3C06 is at rank2 but notes in silva say "cryptophyte"
  #nonphotosynthetic cryptophyceae - Chilomonas and Goniomonas (formerly Cyathomonas)
#Haptophyta kingdom, I think all are Ps
#Centrohelida class - encyclopedia of life said some species have photosynthetic symbionts but I can't find any info on Ps taxa. I don't have many only like 8 taxa so I will leave them as centrohelida and assume they are nonPs. from before: Heterophryidae, Mb-5C (can't find any info on, I'll put it with nonPs), M1-18D08(can't find info),Acanthocystidae(can't find info),Pterocystis(can't find info), so I will leave them all in a group called "centrohelida"
#Kathablepharidae class - nonPs
#NonPs Stramenopiles: MAST-12, Labyrinthulomycetes,Bicosoecida,Peronosporomycetes,Hyphochytriales
#Ps Stramenopiles: Diatomea, Eustigmatales, Xanthophyceae, Chrysophyceae, Raphidophyceae
#heterotrophic eukaryota (things without a kingdom):Telonema (nonPs protist genus);Breviatea (amoeba like);Apusomonadidae sister to the breviatea nonPs

unique(tax_table(datEukS3)[,"Rank2"])
unique(tax_table(datEukS3)[,"Rank3"])

labelsEukS<-tax_table(datEukS3)[,"Rank3"]#

ind<-which(tax_table(datEukS3)[,"Rank2"]=="__Amoebozoa")
labelsEukS[ind]<-"Amoebozoa"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Alveolata")
labelsEukS[ind]<-"Alveolata"
#ind<-which(tax_table(datEukS3)[,"Rank4"]=="__Dinoflagellata")#
#labelsEukS[ind]<-"Photosynthetic_Alveolata"
#ind<-which(tax_table(datEukS3)[,"Rank4"]!="__Dinoflagellata"&tax_table(datEukS3)[,"Rank3"]=="__Alveolata")
#labelsEukS[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(tax_table(datEukS3)[,"Rank2"]=="__Archaeplastida")
labelsEukS[ind]<-"Archaeplastida"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Rhizaria")
labelsEukS[ind]<-"Rhizaria"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Holozoa")
labelsEukS[ind]<-"Holozoa"
ind<-which(tax_table(datEukS3)[,"Rank6"]=="__Euglenida")
labelsEukS[ind]<-"Photosynthetic_Excavata"
ind<-which(tax_table(datEukS3)[,"Rank6"]!="__Euglenida"&tax_table(datEukS3)[,"Rank2"]=="__Excavata")
labelsEukS[ind]<-"Nonphotosynthetic_Excavata"
ind<-which(tax_table(datEukS3)[,"Rank2"]=="__Haptophyta")
labelsEukS[ind]<-"Haptophyta"
ind<-which(tax_table(datEukS3)[,"Rank2"]=="__Centrohelida")
labelsEukS[ind]<-"Centrohelida"
ind<-which(tax_table(datEukS3)[,"Rank2"]=="__Kathablepharidae")
labelsEukS[ind]<-"Kathablepharidae"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS3)[,"Rank4"]%in%c("__MAST-12","__Labyrinthulomycetes","__Bicosoecida","__Peronosporomycetes","__Hyphochytriales"))
labelsEukS[ind]<-"Nonphotosynthetic_Stramenopiles"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS3)[,"Rank4"]%in%c("__Diatomea","__Eustigmatales","__Xanthophyceae","__Chrysophyceae","__Raphidophyceae"))
labelsEukS[ind]<-"Photosynthetic_Stramenopiles"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Stramenopiles"&is.na(tax_table(datEukS3)[,"Rank4"])==T)
labelsEukS[ind]<-"Unknown_Stramenopiles"
ind<-which(tax_table(datEukS3)[,"Rank2"]=="__SAR"&is.na(tax_table(datEukS3)[,"Rank3"])==T)
labelsEukS[ind]<-"Unknown_Eukaryota"
ind<-which(is.na(tax_table(datEukS3)[,"Rank2"])==T)
labelsEukS[ind]<-"Unknown_Eukaryota"
ind<-which(tax_table(datEukS3)[,"Rank3"]=="__Breviatea"|tax_table(datEukS3)[,"Rank3"]=="__Telonema"|tax_table(datEukS3)[,"Rank3"]=="__Apusomonadidae")
labelsEukS[ind]<-"Heterotrophic_Eukarya"
ind<-which(tax_table(datEukS3)[,"Rank2"]=="__Opisthokonta"&is.na(tax_table(datEukS3)[,"Rank3"])==T)
labelsEukS[ind]<-"Heterotrophic_Eukarya"

#unique(tax_table(datEukS3)[which(tax_table(datEukS3)[,"Rank3"]=="__Discoba"),"Rank2"])
#tax_table(datEukS3)[which(tax_table(datEukS3)[,"Rank3"]=="__Apusomonadidae"),]
#ind<-which(is.na(labelsEukS)==T)
#tax_table(datEukS3)[ind,]
unique(labelsEukS)


#I could separate excavata into kingdoms/super phyla and also separate archaeplastida into its 3 kingdoms
unique(labelsEukS)
colnames(labelsEukS)<-"labels"
labelsEukS2<-as.data.frame(labelsEukS)
unique(labelsEukS2$labels)
labelsEukS2$group<-"Eukaryota"

labelsEukS2$group2<-NA
ind<-which(labelsEukS2$labels=="Archaeplastida"|labelsEukS2$labels=="Photosynthetic_Stramenopiles"|labelsEukS2$labels=="Haptophyta"|labelsEukS2$labels=="Photosynthetic_Excavata")
labelsEukS2$group2[ind]<-"PhotosyntheticEukaryota"
ind<-which(labelsEukS2$labels=="Unknown_Eukaryota"|labelsEukS2$labels=="Unknown_Stramenopiles")
labelsEukS2$group2[ind]<-"UnknownEukaryota"
labelsEukS2$group2[which(is.na(labelsEukS2$group2)==T)]<-"HeterotrophicEukaryota"

head(labelsEukS2)

dim(tax_table(datEukS3))

#columns 12 and 13 are all NAs
unique(tax_table(datEukS3)[,11])
labelsEukS2$taxstring<-paste(tax_table(datEukS3)[,1],tax_table(datEukS3)[,2],tax_table(datEukS3)[,3],tax_table(datEukS3)[,4],tax_table(datEukS3)[,5],tax_table(datEukS3)[,6],tax_table(datEukS3)[,7],tax_table(datEukS3)[,8],tax_table(datEukS3)[,9],tax_table(datEukS3)[,10],tax_table(datEukS3)[,11],sep=";")

#replace tax table
tax_table(datEukS3)<-cbind(tax_table(datEukS3),labelsEukS)

#look at tree of abundant taxa
#myTaxa<-c(names(sort(taxa_sums(datEukr2),decreasing=T)))[1:100]
#ex2 = prune_taxa(myTaxa, datEukr2)
#plot_tree(ex2, label.tips = "taxa_names",color="Rank3")






##### Euks N #####

#for lablefiles: column labels is for bargraphs, group is for a simple network, group2 is for photosynthetic/nonphotosynthetic network

#labelsEukN<-data.frame(rep("Mesofauna",dim(tax_table(datEukN3))[1]))

# ind<-which(is.na(labelsEukN))
# labelsEukN[ind]<-tax_table(datEukN3)[ind,"Rank5"]

unique(tax_table(datEukN3)[,"Rank4"])

labelsEukN<-substring(tax_table(datEukN3)[,"Rank4"],3)
ind<-which(is.na(labelsEukN))
labelsEukN[ind]<-"UnknownMetazoa"
 
colnames(labelsEukN)<-"labels"

labelsEukN2<-as.data.frame(labelsEukN)
unique(labelsEukN2$labels)
labelsEukN2$group<-"Mesofauna"

labelsEukN2$group2<-"Mesofauna"

head(labelsEukN2)

#9-13 are all NAs
unique(tax_table(datEukN3)[,8])
labelsEukN2$taxstring<-paste(tax_table(datEukN3)[,1],tax_table(datEukN3)[,2],tax_table(datEukN3)[,3],tax_table(datEukN3)[,4],tax_table(datEukN3)[,5],tax_table(datEukN3)[,6],tax_table(datEukN3)[,7],tax_table(datEukN3)[,8],sep=";")

#replace tax table
tax_table(datEukN3)<-cbind(tax_table(datEukN3),labelsEukN)




##### Bacteria #####

#labelsBac<-data.frame(rep("Bacteria",dim(tax_table(datBacS3))[1]))
#I check the dataset for chemoautotrophs, no methanogens, yes ammonia oxidizers, 
#the crenarchaeota in archaea: the o__Nitrososphaerales is an ammonia oxidizer(all the abundant crenarchaeota are these), the other less abundant orders (like NRP-J and crenarcheales) I'm not sure about but there are chemoorganotrophs in crenarchaeota so I can't assume they are autotrophs
#Ammonia oxidizing bacteria: 
#f__Nitrosomonadaceae (in Proteobacteria)
#p__Nitrospirae in bacteria (in Nitrospiraceae one genus in this group is chemohterotrophic but it is Thermodesulfovibrio but it is thermophillic so it is most likely that the ones present here are the chemoautotrophic ones), there are other candidate(?) families in the Nitrospirales that I can't find much about, I think I will keep them all labeled the same
https://en.wikipedia.org/wiki/Microbial_metabolism#Sulfur_oxidation)
https://link.springer.com/referenceworkentry/10.1007%2F978-3-642-38954-2_126
#The bradyrhizobiaceae (in proteobacteria) is complicated. the only species listed are bradyrhizobium (heterotroph that fixes N, although in papers some strains can be chemoautotrophs) and other species that I dont know what they do or there is no genus/species listed. Some in this family however are photosynthetic. I will leave them all as heterotrophs for now.
#The Planctomycetes are also complicated. there are five genera that do anammox (chemoautotrophic), none of those genera are found in this dataset, however many of the genus/species info are blanks. the class Pla3, Pla4 (and group VI) are anammox according to the following paper (http://aem.asm.org/content/73/15/4707.full). also according to that paper the genera in the o__Gemmatales and Pirellulales orders are all heterotrophs so I will label those as such

unique(tax_table(datBacS3)[,"Rank2"])

labelsBac<-substring(tax_table(datBacS3)[,"Rank2"],4)

ind<-which(tax_table(datBacS3)[,"Rank2"]=="p__")
labelsBac[ind]<-NA
ind<-which(is.na(labelsBac))
labelsBac[ind]<-"Unknown_Bacteria"
#Only the class Chloroflexi are photosynthetic, the other classes (Ktedontobacteria) are not photosynthetic.
ind<-which(tax_table(datBacS3)[,"Rank2"]=="p__Chloroflexi"&tax_table(datBacS3)[,"Rank3"]=="c__")
labelsBac[ind]<-"Unknown_Chloroflexi"
ind<-which(tax_table(datBacS3)[,"Rank2"]=="p__Chloroflexi"&is.na(tax_table(datBacS3)[,"Rank3"])==T)
labelsBac[ind]<-"Unknown_Chloroflexi"
ind<-which(tax_table(datBacS3)[,"Rank3"]=="c__Chloroflexi")
labelsBac[ind]<-"Photosynthetic_Chloroflexi"
ind<-which(labelsBac=="Chloroflexi")
labelsBac[ind]<-"Heterotrophic_Chloroflexi"
#the o__Nitrososphaerales in the Crenarchaeota
ind<-which(labelsBac=="Crenarchaeota")
labelsBac[ind]<-"Unknown_Crenarchaeota"
ind<-which(tax_table(datBacS3)[,"Rank4"]=="o__Nitrososphaerales")
labelsBac[ind]<-"Chemoautotrophic_Crenarchaeota"
#the nitrosomonadaceae are chemoautotrophs and the bradyrhizobiaceae are lots of things
ind<-which(labelsBac=="Proteobacteria")
labelsBac[ind]<-"Heterotrophic_Proteobacteria"
ind<-which(tax_table(datBacS3)[,"Rank5"]=="f__Nitrosomonadaceae")
labelsBac[ind]<-"Chemoautotrophic_Proteobacteria"
#ind<-which(tax_table(datBacS3)[,"Rank5"]=="f__Bradyrhizobiaceae")
#labelsBac[ind]<-"Unknown_Proteobacteria"
ind<-which(labelsBac=="Planctomycetes")
labelsBac[ind]<-"Unknown_Planctomycetes"
ind<-which(tax_table(datBacS3)[,"Rank3"]=="c__Pla4"|tax_table(datBacS3)[,"Rank3"]=="c__Pla3")
labelsBac[ind]<-"Chemoautotrophic_Planctomycetes"
ind<-which(tax_table(datBacS3)[,"Rank4"]=="o__Gemmatales"|tax_table(datBacS3)[,"Rank4"]=="o__Pirellulales")
labelsBac[ind]<-"Heterotrophic_Planctomycetes"


unique(tax_table(datBacS3)[ind,])

unique(labelsBac)
colnames(labelsBac)<-"labels"

labelsBac2<-as.data.frame(labelsBac)
labelsBac2$group<-"Bacteria"

labelsBac2$group2<-NA
ind<-which(labelsBac2$labels=="Photosynthetic_Chloroflexi"|labelsBac2$labels=="Cyanobacteria"|labelsBac2$labels=="Chlorobi") 
labelsBac2$group2[ind]<-"PhotosyntheticBacteria"
ind<-which(labelsBac2$labels=="Unknown_Chloroflexi"|labelsBac2$labels=="Unknown_Bacteria"|labelsBac2$labels=="Unknown_Planctomycetes"|labelsBac2$labels=="Unknown_Crenarchaeota") 
labelsBac2$group2[ind]<-"UnknownBacteria"
ind<-which(labelsBac2$labels=="Chemoautotrophic_Crenarchaeota"|labelsBac2$labels=="Chemoautotrophic_Proteobacteria"|labelsBac2$labels=="Nitrospirae"|labelsBac2$labels=="Chemoautotrophic_Planctomycetes") 
labelsBac2$group2[ind]<-"ChemoautotrophicBacteria"
ind<-which(is.na(labelsBac2$group2)==T)
labelsBac2$group2[ind]<-"HeterotrophicBacteria"
head(labelsBac2)

tax_table(datBacS3)<-cbind(tax_table(datBacS3),labelsBac)

dim(tax_table(datBacS3))
#8 is my column for "bargraph groups"
unique(tax_table(datBacS3)[,7])
labelsBac2$taxstring<-paste(tax_table(datBacS3)[,1],tax_table(datBacS3)[,2],tax_table(datBacS3)[,3],tax_table(datBacS3)[,4],tax_table(datBacS3)[,5],tax_table(datBacS3)[,6],tax_table(datBacS3)[,7],sep=";")



##### Fungi #####

#labelsITS<-data.frame(rep("Fungi",dim(tax_table(datITSS3))[1]))

unique(tax_table(datITSS3)[,"Rank2"])

labelsITS<-substring(tax_table(datITSS3)[,"Rank2"],4)
ind<-which(is.na(labelsITS)==T)
labelsITS[ind]<-"Unknown_Fungi"
unique(labelsITS)
colnames(labelsITS)<-"labels"

labelsITS2<-as.data.frame(labelsITS)
labelsITS2$group<-"Fungi"
labelsITS2$group2<-"Fungi"

head(labelsITS2)

tax_table(datBacS3)<-cbind(tax_table(datBacS3),labelsBac)

unique(tax_table(datITSS3)[,7])
labelsITS2$taxstring<-paste(tax_table(datITSS3)[,1],tax_table(datITSS3)[,2],tax_table(datITSS3)[,3],tax_table(datITSS3)[,4],tax_table(datITSS3)[,5],tax_table(datITSS3)[,6],tax_table(datITSS3)[,7],sep=";")

#replace tax table
tax_table(datITSS3)<-cbind(tax_table(datITSS3),labelsITS)




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

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis1.Rdata") 




##### Phylogenetic diversity ######
# This needs to be done before any label changes below b/c it uses the tree data.
# I need to use root=T, if I use root=F it cannot calculate the diversity in a sample with only one taxon
pdEukS<-pd(as.matrix(datEukS3cotu[,-c(1:31)]),phy_tree(datEukS3c),include.root=TRUE) #took 5 minutes
pdEukN<-pd(as.matrix(datEukN3cotu[,-c(1:31)]),phy_tree(datEukN3c),include.root=TRUE) #took 3 minutes
pdBac<-pd(as.matrix(datBacS3cotu[,-c(1:31)]),phy_tree(datBacS3c),include.root=TRUE) #takes an hour, started 9:48pm, finished 10:50pm
richITS<-as.data.frame(rowSums(datITSS3cotu[,-c(1:31)]>0))

#Just to test
#richBac<-as.data.frame(rowSums(datBacS3cotu[,-c(1:31)]>0))



###### Grouping by kingdom/phylum #####
#for the bar graphs

datEukS3k<-aggregate.data.frame(otu_table(datEukS3),by=list(labels=tax_table(datEukS3)[,"labels"]),sum)
rownames(datEukS3k)<-datEukS3k$labels
datEukS3k$labels<-NULL
datEukS3k2<-cbind(sample_data(datEukS3),t(datEukS3k))
head(datEukS3k2)

datEukN3k<-aggregate.data.frame(otu_table(datEukN3),by=list(labels=tax_table(datEukN3)[,"labels"]),sum)
rownames(datEukN3k)<-datEukN3k$labels
datEukN3k$labels<-NULL
datEukN3k2<-cbind(sample_data(datEukN3),t(datEukN3k))
head(datEukN3k2)

datBacS3k<-aggregate.data.frame(otu_table(datBacS3),by=list(labels=tax_table(datBacS3)[,"labels"]),sum)
rownames(datBacS3k)<-datBacS3k$labels
datBacS3k$labels<-NULL
datBacS3k2<-cbind(sample_data(datBacS3),t(datBacS3k))
head(datBacS3k2)

datITSS3k<-aggregate.data.frame(otu_table(datITSS3),by=list(labels=tax_table(datITSS3)[,"labels"]),sum)
rownames(datITSS3k)<-datITSS3k$labels
datITSS3k$labels<-NULL
datITSS3k2<-cbind(sample_data(datITSS3),t(datITSS3k))
head(datITSS3k2)
  
  
  
  


###### Filter data sets for network analysis ######
#Follwing Widder et al 2014 PNAS

#take out doubletons and singletons. this doesn't really matter b/c I will take out taxa with 7 or fewer occurrences later on, Im just doing it now to make the dataframe a little smaller and more manageable
ind<-(which(colSums(datEukS3otu[,32:dim(datEukS3otu)[2]]>0)>2))+31
datEukS3otu2<-cbind(datEukS3otu[,1:31],datEukS3otu[,ind])
datEukS3cotu2<-cbind(datEukS3cotu[,1:31],datEukS3cotu[,ind])

ind<-(which(colSums(datEukN3otu[,32:dim(datEukN3otu)[2]]>0)>2))+31
datEukN3otu2<-cbind(datEukN3otu[,1:31],datEukN3otu[,ind])
datEukN3cotu2<-cbind(datEukN3cotu[,1:31],datEukN3cotu[,ind])

ind<-(which(colSums(datBacS3otu[,32:dim(datBacS3otu)[2]]>0)>2))+31
datBacS3otu2<-cbind(datBacS3otu[,1:31],datBacS3otu[,ind])
datBacS3cotu2<-cbind(datBacS3cotu[,1:31],datBacS3cotu[,ind])

ind<-(which(colSums(datITSS3otu[,32:dim(datITSS3otu)[2]]>0)>2))+31
datITSS3otu2<-cbind(datITSS3otu[,1:31],datITSS3otu[,ind])
datITSS3cotu2<-cbind(datITSS3cotu[,1:31],datITSS3cotu[,ind])

#filter out taxa that have a summed relative abundance of <.002 (.2%)
#I tested this and found that there were about 35 bacteria who had 8 or more occurrences but still a summed rel abun of <.002, so this step does remove some taxa from the network independent of the occurrence removal
ind<-(which(colSums(datEukS3otu2[,32:dim(datEukS3otu2)[2]])>0.002))+31
datEukS3otu3<-cbind(datEukS3otu2[,1:31],datEukS3otu2[,ind]) #
datEukS3cotu3<-cbind(datEukS3cotu2[,1:31],datEukS3cotu2[,ind]) #1091 otu

ind<-(which(colSums(datEukN3otu2[,32:dim(datEukN3otu2)[2]])>0.002))+31
datEukN3otu3<-cbind(datEukN3otu2[,1:31],datEukN3otu2[,ind]) #142 otu
datEukN3cotu3<-cbind(datEukN3cotu2[,1:31],datEukN3cotu2[,ind]) #

ind<-(which(colSums(datBacS3otu2[,32:dim(datBacS3otu2)[2]])>0.002))+31
datBacS3otu3<-cbind(datBacS3otu2[,1:31],datBacS3otu2[,ind]) #4853 otu
datBacS3cotu3<-cbind(datBacS3cotu2[,1:31],datBacS3cotu2[,ind]) #

ind<-(which(colSums(datITSS3otu2[,32:dim(datITSS3otu2)[2]])>0.002))+31
datITSS3otu3<-cbind(datITSS3otu2[,1:31],datITSS3otu2[,ind]) #1122 otu
datITSS3cotu3<-cbind(datITSS3cotu2[,1:31],datITSS3cotu2[,ind]) #

#order of doing things: filtered out unwanted taxa (chloroplasts, spiders) and samples (S.2015), rarefied, relativized, took out doubletons and singletons, took out samples <2% summed abundance. I am not going to rarefy or relativize here again b/c the doubletons/singletons/.2% otus that I removed are real, I could have included them in the network analyis if I wanted.


##### Plants #####

plantcomp<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Plants/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots; there are some plots that have plant data but not microbe data. 69 70 71 77 81 108 117 118 147 148 149 151. This is because when we started doing the surveys we were going to all plots for plants and only sample some for microbes, then we realized that that was insane!
dim(plantcomp)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantcomp2$LICHEN<-NULL


###### Make labelfile ######

#Test if there is name overlap in full otu tables, no names overlap, no need to rename OTU columnson the full otu table files, however it is nice for figuring out what is going on in the modeling, so I will add B, N, S, I before all the otu names for the reduced otu table files
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

# namesEukS2 <- sub("^", "S", namesEukS)
# namesEukN2 <- sub("^", "N", namesEukN)
# namesBac2 <- sub("^", "B", namesBac)
# namesITS2 <- sub("^", "I", namesITS)
# 
# names(datEukS3otu)[-c(1:31)]<-namesEukS2
# names(datEukN3otu)[-c(1:31)]<-namesEukN2
# names(datBacS3otu)[-c(1:31)]<-namesBac2
# names(datITSS3otu)[-c(1:31)]<-namesITS2


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
names(datEukS3cotu3)[-c(1:31)]<-namesEukS2
names(datEukN3otu3)[-c(1:31)]<-namesEukN2
names(datEukN3cotu3)[-c(1:31)]<-namesEukN2
names(datBacS3otu3)[-c(1:31)]<-namesBac2
names(datBacS3cotu3)[-c(1:31)]<-namesBac2
names(datITSS3otu3)[-c(1:31)]<-namesITS2
names(datITSS3cotu3)[-c(1:31)]<-namesITS2


#Make combined labelfile
labelsEukS2$otu<-sub("^", "S",rownames(labelsEukS2))
labelsEukN2$otu<-sub("^", "N",rownames(labelsEukN2))
labelsBac2$otu<-sub("^", "B",rownames(labelsBac2))
labelsITS2$otu<-sub("^", "I",rownames(labelsITS2))

labelfile1<-rbind(labelsEukS2,labelsEukN2,labelsBac2,labelsITS2)
head(labelfile1)
labelfile1$oldotu<-rownames(labelfile1)
#labelfile1$oldotu<-substring(rownames(labelfile1), 2)

#combine with plant labelfile
labelsPlant<-as.data.frame(cbind(labels="Plant",group="Plant",group2="Plant",taxstring=colnames(plantcomp2)[2:55],otu=colnames(plantcomp2)[2:55],oldotu=colnames(plantcomp2)[2:55]))
head(labelsPlant)

labelfile<-rbind(labelfile1,labelsPlant)
tail(labelfile)
















#notes: it looks like both prune_samples and subset_samples do NOT remove taxa that have an abundance of 0 after the samples are removed.



