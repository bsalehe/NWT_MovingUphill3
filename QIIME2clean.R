##QIIME2

#notes: this is for installing an old version of qiime2, be sure to find the current release and modify the code accordingly

##First I needed to uninstall and then reinstall miniconda3. it would not update by itself with the update code
https://conda.io/docs/user-guide/install/macos.html

#The second time I tried (when reinstalling for the patch) it worked
conda update conda
conda install wget

##How to install qiime2
https://docs.qiime2.org/2017.12/install/native/

wget https://data.qiime2.org/distro/core/qiime2-2018.2-py35-osx-conda.yml
conda env create -n qiime2-2018.2 --file qiime2-2018.2-py35-osx-conda.yml

# OPTIONAL CLEANUP
rm qiime2-2018.2-py35-osx-conda.yml


#How to uninstll qiime2
#I needed to do this to download the new patch for the ITS dada2 error
conda env remove -n qiime2-2017.12

#Activate Qiime env, must do this in every new terminal tab that is opened
source activate qiime2-2018.2

#Test the new environment
qiime --help

#Paired end tutorial
https://docs.qiime2.org/2017.12/tutorials/atacama-soils/
#General tutorial
https://docs.qiime2.org/2017.12/tutorials/moving-pictures/


#discussion about whether it is worth it to join paired end reads when the reads completely overlap
https://forum.qiime2.org/t/question-about-dada2-denoise-paired-analysis/464
#antoher discussion of relaxing the maxEE filtering parameter when using full reads with poor quality scores near the ends (in qiime dada2 denoise-paired). however the default does delete everyting after the first instance of a bp with quality score of 2 (this can be changed as well)
https://github.com/qiime2/q2-dada2/issues/48






##### Taxonomic databases #####

##### Silva #####
#notes about the silva release: https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_notes.txt
#workflow from this thread: https://forum.qiime2.org/t/18s-classifier-using-silva-database-and-emb-primers/361

#With all taxa and the 111 release, start 1:31pm, end 1:33
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path Silva_111_post/rep_set/Silva_111_full_unique.fasta \
--output-path all_SILVA111_unique_otus.qza


#Import taxonomy, only with 18S.
# There are lots of options here, there is a file that says "no ambiguous" but I dont know if they mean, I assume no ambiguous taxa. there is also a file with consistent 6 ranks (the help said that RDP requires all taxa to have consistent number of ranks. I can try this if the full taxonomy doesn't work). there are also 99 clustered files, but I will use full database
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path Silva_111_post/taxonomy/Silva_111_taxa_map_full.tsv \
--output-path all_SILVA111_unique_taxonomy.qza


##### Extract EMBP variable region for 18S #####
#with single end data you should set a length. When I was playing around I used 150 b/c I was following the person writing the workflow. I could use the max length from my euk data below (which is generally less than 150, except there are about 2000 reads per sample that are 301)
#with paired end data (what I ended up useing), trunc-len should not be used (set at default 0). post saying that you shouldn't truncate with paired end data: https://forum.qiime2.org/t/how-can-i-train-classifier-for-paired-end-reads/1512/7, I think truncation is for when you have only a forward read and it is a particular length, start 1:32pm, end 2:14pm

#Without the 150 truncation, start 2:30pm, end 3:17pm
qiime feature-classifier extract-reads \
--i-sequences all_SILVA111_unique_otus.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads ref-seqs_all_unique_SILVA111.qza

#export to see what it trimmed. 
#Most lengths are less than 140, there are some longer though
qiime tools export \
ref-seqs_all_unique_SILVA111.qza \
--output-dir exported-ref-seqs_all_unique_SILVA111

awk '{print length}' dna-sequences.fasta | sort | uniq -c


#Train classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_unique_SILVA111.qza \
--i-reference-taxonomy all_SILVA111_unique_taxonomy.qza \
--o-classifier all_EMB_SILVA111_classifier.qza


#talking about what the different outputs (truncated taxonomy vs blank genus/species levels means: 
https://forum.qiime2.org/t/consensus-blast-taxonomy-strings/586/2




##### Greengenes database ######
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path gg_13_8_otus/rep_set/99_otus.fasta \
--output-path gg_13_8_otus_99_otus.qza

#Import taxonomy
#Need to duplicate and save the txt file as tsv
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.tsv \
--output-path gg_13_8_otus_99_taxonomy.qza

#start 3:07pm end 3:24pm
#note: when I did this, I was assuming I was going to use paired end reads, so I did not set --p-trunc-len, however then I ended up using just the forward read, since lots of reads were being discarded b/c of poor read quality on the reverse read. I guess I never went back and re-extracted the database with a length. it seemed to work ok though.
qiime feature-classifier extract-reads \
--i-sequences gg_13_8_otus_99_otus.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--o-reads ref-seqs_all_99_gg_13_8.qza

#export to see what it trimmed
qiime tools export \
ref-seqs_all_99_gg_13_8.qza \
--output-dir exported-ref-seqs_all_99_gg_13_8

cd exported-ref-seqs_all_99_gg_13_8
awk '{print length}' dna-sequences.fasta | sort | uniq -c
#Most lengths are 253, there are some longer and shorter from about 251-257

##Train classifier
#start 3:35pm, end 3:42pm
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_gg_13_8.qza \
--i-reference-taxonomy gg_13_8_otus_99_taxonomy.qza \
--o-classifier all_EMB_gg_13_8_classifier.qza




##### UNITE database #####
#I downloaded the file 12_11 alpha release from the qiime page here: http://qiime.org/home_static/dataFiles.html I didn't see it on the unite website...
#tutorial https://github.com/gregcaporaso/2017.06.23-q2-fungal-tutorial
#I tried downloading the UNITE release and processing it but I kept getting errors, I think there were a bunch of formatting errors in the file, so I downloaded an already trained classifir from unite here (UNITE version 7.2): https://forum.qiime2.org/t/unite-ver-7-2-2017-12-01-classifiers-for-qiime2-ver-2017-12-available-here/3020
#She followed the emp protocol for training it

unite-ver7-99-classifier-01.12.2017.qza



##### Euk #####

##### Demultiplexing #####

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/Euks

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip reverse.fastq 
gzip forward.fastq

#start 4:08pm, end 4:15
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences \
--output-path emp-paired-end-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv

#barcode is on the forward read (see word doc showing how the R1 has the reverse primers in the read, while the R2 has the forward primers, this is b/c the sequence is so short that it kept reading the other half of the primer. However the mapping file is correct: "reverse primer" matches to the R2). looking at the barcodes and the barcode file, it looks like it is the reverse complement
#start 5:20pm, end 6:30
qiime demux emp-paired \
--m-barcodes-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize, start 8:44pm, end 9:03pm
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#Export it so you can look inside at the files
#started 8:58pm, end 9:00pm
qiime tools export \
demux.qza \
--output-dir exported-demux
  
##### Trim primers/adapters #####
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 9:44pm, end 10:02
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f GTAGGTGAACCTGCAGAAGGATCA \
--p-adapter-r GACGGGCGGTGTGTAC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose
#note that it said that for the reverse read, the primer was preceeded by C extremely often, which might mean it is part of the primer. it shouldn't be though, so I'm leaving it in.
#note I could have used --p-error-rate 0 \  meaning no errors in the adapter, the default used above is 10%, not sure what difference it would make

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs
#look at this tomorrow and make sure it trimmed things correctly!!! Yes they are, looks good.
awk '{print length}' N.0.2015_49_L001_R1_001copy.fastq | sort | uniq -c

# note that for some of the forward reads start with N for a basepair, this is actually "true" when you look at the reverse and forward reads, the N represents a certain bp that was apparently unknown in the sequencing. (also when it starts with a g rather than an n (gctac) the g blasts to something so the g is correct). The DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html states that no N are allowable in DADA2, so I need to get rid of that first basepair.

grep --color -n "^N"  N.0.2015_49_L001_R1_001copy.fastq # ^ means at the beginning of the line. there actually aren't that many, only like 30 reads start with N
grep --color -n "^N"  N.0.2015_49_L001_R2_001copy.fastq #not any Ns at the beginning of reads like in R1

##### Denoising with DADA2 ##### 
#220 and 210 respectively are when the median quality scores hit 25, however reading the documentation it says reads that are shorter than the truncation numbers are discarded (this seems silly but it is what it is), thus I should keep everything (I tried doing truncation prior to dada2 but it didn't work). Also, I don't think it matters as much for euks b/c the read is so short, if the primer is successfully removed, there is no issue with quality. reads are ~130bp, primer is 24 or 16bp, so the primers are well within the good quality read.
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#the p-trim-left-f is set at 1 to get rid of the N basepair
#start 8:11pm, end 9:22
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 6 \
--p-trim-left-f 1 \
--p-trim-left-r 0 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0

#export rep seqs just to take a look
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat filename|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 244bp, and blasts to a fungus

#export table to get otu table
qiime tools export \
table.qza \
--output-dir exported-table

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
sed '2s/[ ]//' otu_table.txt | sed '2s/.//' > otu_table2.txt

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv


###### Create a phylogenetic tree #####

#do the alignment
#start 9:56pm, end  10:00. could add --p-n-threads to the code
qiime alignment mafft \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza

#Mask (or filter) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree.
#start 10:30pm, end 10:38
qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

#generate tree from masked alignment
#start 10:39, end 10:45
qiime phylogeny fasttree \
--i-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza

#apply midpoint rooting to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree
#start 2:04pm, and 2:05pm 
qiime phylogeny midpoint-root \
--i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#export
qiime tools export \
rooted-tree.qza \
--output-dir exported-rooted-tree

###### Assign taxonomy #####

#note: I tried using blast just to see what I would get. I accepted a lot of defaults here because I couldn't figure out what they were in the last 97% analysis. These defaults might actually be wrong b/c thre is one about percent identity whose default is 80%. There is a post saying blast is terrible and will classify anything, no matter how bad the fit is: https://github.com/benjjneb/dada2/issues/323
#I don't think I'll use blast or uclust, I feel like they work with OTU clustering methods, not with amplicon sequence variant analysis like DADA2.

#start 3:25, end 3:36
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA111_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy6.qza


#explanation of sklearn and training 
#https://forum.qiime2.org/t/classifier-training-questions/1162/3

qiime tools export \
taxonomy6.qza \
--output-dir exported-taxonomy6

#navigate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy6.qza \
--o-visualization taxonomy6.qzv

qiime tools view taxonomy6.qzv


#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy6.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots6.qzv

qiime tools view taxa-bar-plots6.qzv
#you can download a csv file from this visualization, might be useful

#I did a bunch of trials with assigning taxonomy using different databases and classifiers (code not shown above) here are my notes: there are unassigned taxa in the file that was blasted against the all database, there are no unassigned taxa in the file that blasted agains the Euk only database. this might be because when you blast against euks it just leavs out anything that doesnot have a hit (b/c they are also likely bacteria so it doesn't make sense to call them unassigned ??)
#looking only at eukaroyte classified reads, many of the numbers in the groups (level2) are identical, a few are a little different, with the Euk only database always having more reads (if they are not the same). the euk dataset also has many many more reads for the Eukaryota;__ classification. That is where the main difference is. This is because it looks like any read that is unclassified or classified as bacteria in the dataset (from the all database) is calssified as a Eukaryota;__ in the euk only database. (i.e. the total number of reads is the same). Looking at some of the taxonomic groups at Level 2, some of the numbers of reads are higher in the all data base, some higher in the euk-only database. I can't explain that, except it is a different training so there is variability
#comparing concensus and majority at level 2 they are nearly identical, does not affect the number of reads in Eukaryota;__
#comparing truncated vs not, they look almost identical
#comparing release 128 to 111, 111 has fewer unassigned, about the same number of Eukaryota;__





##### Bact - processing only forward reads #####

##### Demultiplexing #####

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/Bactsingle

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip sequences.fastq 

#start 10:17, end 10:26
qiime tools import \
--type EMPSingleEndSequences \
--input-path emp-single-end-forward-sequences \
--output-path emp-single-end-forward-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv

#barcode is on the forward read, it is not the reverse complement
#start 10:28 pm, end 11:36 pm
#the first time I did this, I got an error:
Plugin error from demux:
  Mismatched sequence ids: M01918:229:000000000-ALYY5:1:2102:16350:8485 and ACAGAGGAGAG7:CAACTGEGGG<E*CBT8GGGGGGGGTAGGTCCG,:@EGGGDEB;TGCGFFGGGGGEGGGGGGGGG4:00-ALYY5:1:2102:1GGGGGGG18:229:0000C*CGCTGGCTGACG3EFTAAAGGGFGF?,CB+C8*GGGGG0CGGGGGGCGFGGGGGFGGGGGGGCACCTATCCTTGCGCAGGGGGCGCACCTGDGGGGGFGF?,GGGGGGGGGGGG7GF+CGTGG25:1:2TGTAGCGGTGGAATG8>5/CCFBGGFGGGGGGGGGGG2/CCFBG::::GGGGFGGGGGGGGGAAGTGGGGGGGGGGG)7GGGGGGGGGG2
#then I realized that I couldn't gunzip the sequences file (it was corrupted or something) so I re-created the gzip file from the raw .fastq file and then it ran fine

#Second time I tried this, I did demultiplexing and then trimming. during the trimming, I got a crazy error that I couldn't find any info about, but it mentioned sample N.47.2015 so I went back to the mapping file and deleted N.47.2015, then I redid demultiplexing and trimming and it was fine.

qiime demux emp-single \
--m-barcodes-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--m-barcodes-column BarcodeSequence \
--i-seqs emp-single-end-forward-sequences.qza \
--o-per-sample-sequences demux 

#summarize, start 6:21pm, end
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#export it so you can look inside at the files
#started 6:22pm, end
qiime tools export \
demux.qza \
--output-dir exported-demux


##### Trim primers/adapters ##### 
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 1:24am, end 1:37pm, ran fine, no error
qiime cutadapt trim-single \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter ATTAGAWACCCBNGTAGTCC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv
#nice, the samples with low # sequences here matches that in the DataCleaning file (samples 5, 34, 126)

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs

#Getting a histogram of line lengths. at least for S.99.2015, the "histogram" below is exactly the same as for paired reads
awk '{print length}' S.0.2015_5_L001_R1_001copy.fastq | sort | uniq -c
#read 1 80036 at 253, 572 at 252, 84 at 251, 2 at 250

awk '{print length}' S.99.2015_101_L001_R1_001copy.fastq | sort | uniq -c
#read 1 96184 at 253, 1066 at 252, 224 at 251, 34 at 250, 3862 at 301


##### Denoising with DADA2 #####
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 275bp is when the median quality scores hit 25 - based on previous r script. using visualizaiton above (after primer trimming - the cutoffs for quality 25 is 254)
#the reads I'm seeing are mostly 253bp, if the primers are 20 and 19 bp, then for R1, the primer should be within the high quality region, so should be taken off but there were some reads still at 301bp. thus trimming at 253 should take the primer off even if there were lots of bp errors in it. However, I will truncate at 251 to delete the potentially bad quality bp without losing too many small fragments
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#start 4:41am, end 1:23pm
qiime dada2 denoise-single \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 0 \
--p-trim-left 0 \
--p-trunc-len 251 \

#export rep seqs just to take a look at. 46,673 rep sequences
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat dna-sequences.fasta|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 460bp, and blasts to a fungus, odd

#26277 are 253 bp
awk '{print length}' dna-sequences.fasta | sort | uniq -c

#export table to get otu table
qiime tools export \
table.qza \
--output-dir exported-table

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
sed '2s/[ ]//' otu_table.txt | sed '2s/.//' > otu_table2.txt

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv


##### Create a phylogenetic tree #####

#do the alignment
#start 4:22pm, end 5:52pm. --p-n-threads -1 means use all cores
qiime alignment mafft \
--i-sequences rep-seqs.qza \
--p-n-threads -1 \
--o-alignment aligned-rep-seqs.qza

#Mask (or filter) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree.
#start 6:09pm, end 6:56
qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

#generate tree from masked alignment
#start 8:12, end 9:10
qiime phylogeny fasttree \
--i-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza

#apply midpoint rooting to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree
#start 1:04pm, and 1:04pm 
qiime phylogeny midpoint-root \
--i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#export
qiime tools export \
rooted-tree.qza \
--output-dir exported-rooted-tree


##### Assign taxonomy #####

#using greengenes
#start 1:50pm, end 2:13pm
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_gg_13_8_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy_gg.qza

qiime tools export \
taxonomy_gg.qza \
--output-dir exported-taxonomy_gg


#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy_gg.qza \
--o-visualization taxonomy_gg.qzv

qiime tools view taxonomy_gg.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_gg.qza \
--m-metadata-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--o-visualization taxa-bar-plots_gg.qzv

qiime tools view taxa-bar-plots_gg.qzv












###### ITS - Processing only the forward reads #####

#in directory QIIME2/ITSsingle/
#In the previous analysis with QIIME1, we just used the reverse read,I think b/c that is where the barcode was. Here i am using the forward read b/c it was better quality, there was one paper (in qiime2 directory) who also used that excuse to choose the forward read
#importing only the forward read 
#the sequences need to be called "sequences.fastq.gz"
qiime tools import \
--type EMPSingleEndSequences \
--input-path emp-single-end-forward-sequences \
--output-path emp-single-end-forward-sequences.qza

###### Demultiplexing single reads ######
#start 11:43pm, end 12:27am
qiime demux emp-single \
--m-barcodes-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-column BarcodeSequence \
--i-seqs emp-single-end-forward-sequences.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize, start 8:45am, end 8:
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#export it so you can look inside at the files
#started 8:54, end 9:08
qiime tools export \
demux.qza \
--output-dir exported-demux

#The amplicon size should be 230 bp (according to the earth microbiome website). So there will be primers in many reads.

###### Trim primers/adapters ######
#To get these adapters I looked in the raw sequence data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 9:50am, end 9:58am
qiime cutadapt trim-single \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter GCATCGATGAAGAACGCAGC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs

#Getting a histogram of line lengths
awk '{print length}' S.99.2015_141_L001_R1_001_copy.fastq | sort | uniq -c
#read 1 lots of variability, but most around 213 (4000 reads) and 235 (3000) and 301 (13834 reads)

##### Denoising with DADA2 #####
#Trimming to remove poor quality reads - if I wanted to use fastx_trimmer - When I used it before, I did it on the initial file from the sequencing company, prior to demultiplexing (a fastq file). I might be able to do it on multiple files at a later stage (after demultiplexing), but then I'd have to figure out how to do it on many files (each sample individually) and then convert those file back to a qza file. Not impossible but maybe not worth it. I tried trimming on the initial fastq file with fastx_trimmer and got through demultiplexing, but then got an error when I tried to trim adapters using cut-adapt. it could have been that there was a particular sample that was messed up or it could be genearlly that cut-adapt is not able to deal with input reads of different lengths. the error was about the length of the read and length of the quality scores not being the same.
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 290bp is when the median quality scores hit 25 - based on previous r script. using visualizaiton above, however the reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step.
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#Originally I was getting a bunch of errors in denoising, but then I found that they released a patch: https://forum.qiime2.org/t/qiime-2-2017-12-release-is-now-live/2308/10?u=thermokarst
#I just needed to uninstall then reinstall qiime2. and it was fine, no more errors
#start 12:04pm, end 12:54pm
qiime dada2 denoise-single \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 0 \
--p-trim-left 0 \
--p-trunc-len 0

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat dna-sequences.fasta|awk '{print length, $0}'|sort -nr|head -1

#most are between 214-251 bp, but a lot at 301 too
awk '{print length}' dna-sequences.fasta | sort | uniq -c

#export table to get otu table
qiime tools export \
table.qza \
--output-dir exported-table

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
sed '2s/[ ]//' otu_table.txt | sed '2s/.//' > otu_table2.txt

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv


##### Assign taxonomy #####

#start 11:47pm, end 12:13 pm 
qiime feature-classifier classify-sklearn \
--i-classifier unite-ver7-99-classifier-01.12.2017.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy_unite.qza

qiime tools export \
taxonomy_unite.qza \
--output-dir exported-taxonomy_unite

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy_unite.qza \
--o-visualization taxonomy_unite.qzv

qiime tools view taxonomy_unite.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_unite.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots_unite.qzv

qiime tools view taxa-bar-plots_unite.qzv







