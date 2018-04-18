##QIIME2

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

#OTUs
#first time
SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta
#second and third time - all taxa (bacteria and euks)
SILVA_128_QIIME_release/rep_set/rep_set_all/99/99_otus.fasta

#Taxonomy
#Explanation for the different taxonomy files in SILVA https://groups.google.com/forum/#!topic/qiime-forum/TqztahqRFSk
#first time
SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_all_levels.txt
#second time - 
SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_all_levels.txt
#third time - 
SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/majority_taxonomy_all_levels.txt

#Import OTUs, try it with the 18S only first
#start 6:23pm, end 6:24
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta \
--output-path 18S_SILVA128_99_otus.qza

#second and third time, with all taxa, start 7:33pm, end 7:36
##USE THIS if you want 128 release
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_128_QIIME_release/rep_set/rep_set_all/99/99_otus.fasta \
--output-path all_SILVA128_99_otus.qza

#6th time, with all taxa and the 111 release, start 1:31pm, end 1:33
##USE THIS
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path Silva_111_post/rep_set/Silva_111_full_unique.fasta \
--output-path all_SILVA111_unique_otus.qza


#Import taxonomy, only with 18S. I'm using consensus taxonomy, which means all of the sequences need to have the same exact species taxonomy to be identified as a species, otherwise it is identified as "ambiguous taxa". This is very conservative. the alternative is using majority_taxonomy which 90% of the sequences in the database have to match. One comment in a forum is that the majority should be conservative enough for most uses, but if you are looking into a particular taxon, you should probably take the sequence and blast it yourself separately to make sure you are referencing the correct taxonomy. B/c I might do this when looking at hubs/connectors, I'll leave it a the highly conservative consensus, since I don't want to reblast things and I don't want to mis-report.
#Need to duplicate and save the txt file as tsv
#start 6:50, end 6:50
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_7_levels.tsv \
--output-path 18S_SILVA128_99_taxonomy.qza

#second time, all taxa. this is odd, there is also a file called just "taxonomy_all_levels.txt" not sure how that is different from consnsus_taxonomy_all_levels.txt  In consensus taxonomy, identical sequences with different names are renamed using consensus
#USE THIS if you want 128 release
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_all_levels.tsv \
--output-path all_SILVA128_99_taxonomy.qza

#third time, all taxa majority or I could also try the vsearch classifier as that is comparable to uclust https://forum.qiime2.org/t/import-i-reference-taxonomy-taxonomy-tsv-to-qza/781/14. i used blast though. or I could do a different train classifier thing below
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/majority_taxonomy_all_levels.tsv \
--output-path all_SILVA128_99_taxonomy_majority.qza

#fourth time, the file called just "taxonomy_all_levels.txt" 
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/taxonomy_all_levels.tsv \
--output-path all_SILVA128_99_taxonomy_plain.qza

#sixth time (not sure what happened with 5). THere are lots of options here, there is a file that says "no ambiguous" but I dont know if they mean ambiguous taxa or bases, I assume taxa since this is a taxonomy file not base pairs. there is also a file with consistent 6 ranks (the help said that RDP requires all taxa to have consistent number of ranks. I can try this if the full taxonomy doesn't work). there are also 99 clustered files, but I will use full database
#USE THIS if you want 111 release
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path Silva_111_post/taxonomy/Silva_111_taxa_map_full.tsv \
--output-path all_SILVA111_unique_taxonomy.qza


##### Extract EMBP variable region for 18S #####
#I'm not entirely sure what length to use, probably should use the max length from my euk data below (which is generally less than 150, except there are about 2000 reads per sample that are 301) or 150, since that is exactly what the person writing this workflow used. update: with paired end data, trunc-len should not be used (set at default 0), see below
#start 7:12, end 7:17
qiime feature-classifier extract-reads \
--i-sequences 18S_SILVA128_99_otus.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--p-trunc-len 150 \
--o-reads ref-seqs_18S_99_SILVA128_RL150.qza

#start 7:45, end 8:29 (way slower! so many more taxa!)
qiime feature-classifier extract-reads \
--i-sequences all_SILVA128_99_otus.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--p-trunc-len 150 \
--o-reads ref-seqs_all_99_SILVA128_RL150.qza

#Try it wihtout the 150 truncation (default is 0). post saying that you shouldn't truncate with paired end data: https://forum.qiime2.org/t/how-can-i-train-classifier-for-paired-end-reads/1512/7, I think truncation is for when you have only a forward read and it is a particular length, start 1:32pm, end 2:14pm
##USE THIS for 128 release
qiime feature-classifier extract-reads \
--i-sequences all_SILVA128_99_otus.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads ref-seqs_all_99_SILVA128.qza


#Without the 150 truncation, start 2:30pm, end 3:17pm
##USE THIS for 111 release
qiime feature-classifier extract-reads \
--i-sequences all_SILVA111_unique_otus.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads ref-seqs_all_unique_SILVA111.qza


#export to see what it trimmed. there are definitely sequences that are 101 bp, thus lower than 150, so it's not like it deleted short sequences. 150 is the largest length
qiime tools export \
ref-seqs_18S_99_SILVA128_RL150.qza \
--output-dir exported-ref-seqs_18S_99_SILVA128_RL150

awk '{print length}' dna-sequences.fasta | sort | uniq -c

#Most lengths are less than 140, there are some longer though
#USE THIS
qiime tools export \
ref-seqs_all_unique_SILVA111.qza \
--output-dir exported-ref-seqs_all_unique_SILVA111.qza

awk '{print length}' dna-sequences.fasta | sort | uniq -c

#Train classifier
#start 8:02, end 8:03
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_18S_99_SILVA128_RL150.qza \
--i-reference-taxonomy 18S_SILVA128_99_taxonomy.qza \
--o-classifier 18S_EMB_SILVA128_classifier.qza

#start 10:25, end 10:26, renamed added RL150 at the end
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_SILVA128_RL150.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy.qza \
--o-classifier all_EMB_SILVA128_classifierRL150.qza

#start
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_SILVA128_RL150.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy_majority.qza \
--o-classifier all_EMB_SILVA128_classifier_majority.qza

#start
###USE THIS for 128
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_SILVA128.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy.qza \
--o-classifier all_EMB_SILVA128_classifier.qza

#start
###USE THIS for 111
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_unique_SILVA111.qza \
--i-reference-taxonomy all_SILVA111_unique_taxonomy.qza \
--o-classifier all_EMB_SILVA111_classifier.qza


#talking about what the different outputs (truncated taxonomy vs blank genus/species levels means: 
https://forum.qiime2.org/t/consensus-blast-taxonomy-strings/586/2



##### Extract region for 16S #####
#With paired end data, trunc-len should not be used (set at default 0), see below

#from the euk folder copy all_SILVA128_99_otus.qza into the Bact folder
# and copy all_SILVA128_99_taxonomy.qza

#start 5:46pm, end 6:39pm
qiime feature-classifier extract-reads \
--i-sequences all_SILVA128_99_otus.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--o-reads ref-seqs_all_99_SILVA128.qza

#export to see what it trimmed
qiime tools export \
ref-seqs_all_99_SILVA128.qza \
--output-dir exported-ref-seqs_all_99_SILVA128

cd exported-ref-seqs_all_99_SILVA128
awk '{print length}' dna-sequences.fasta | sort | uniq -c
#Most lengths are 253, there are some longer and shorter from about 249-259

##Train classifier, this may have taken a long time, I forgot to log the times
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_SILVA128.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy.qza \
--o-classifier all_EMB_SILVA128_classifier.qza



##### Extract region for ITS ####
#They suggest not trimming the ref sequences prior to training the sklearn classifier for ITS: https://forum.qiime2.org/t/working-with-its-data-on-a-supercomputer-server/1061/15

##Train classifier
#start 9:25ish am, end 11:28pm (laptop)
#start 5:21pm (office)
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads all_SILVA128_99_otus.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy.qza \
--o-classifier all_EMB_SILVA128_classifier.qza



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
#I downloaded the file 12_11 alpha release from the qiime page here: http://qiime.org/home_static/dataFiles.html I didn't see it on th unite website...
#tutorial https://github.com/gregcaporaso/2017.06.23-q2-fungal-tutorial

#Import OTUs

#Error

#First correct some lowercase bp letters:
  grep 'a' 99_otus.fasta 
#https://forum.qiime2.org/t/plugin-error-feature-classifier-classify-sklearn/1532/2
tr 'acgt' 'ACGT' < 99_otus.fasta > 99_otus_uppercase.fasta 

#Then change the '_' to '-' note: i am assuming this is what _ means, it is at the beginning of one read and end of another read in the file
tr '_' '-' < 99_otus_uppercase.fasta > 99_otus_uppercase2.fasta 

#delete the lines with: '\x0b' and the header before it
#grep -v "\x0b" 99_otus_uppercase2.fasta > 99_otus_uppercase3.fasta #Only deletes the sequenc line
#arg i coudl just use an older release that actually works

#start 8:49pm, end 8:49
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path its_12_11_otus/rep_set/99_otus_uppercase3.fasta \
--output-path UNITE_12_11_99_otus.qza

#Import taxonomy
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path its_12_11_otus/taxonomy/99_otu_taxonomy.txt \
--output-path UNITE_12_11_99_taxonomy.qza

#no need to extract reference region

#train classifier
#start 11:35pm
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads UNITE_12_11_99_otus.qza \
--i-reference-taxonomy UNITE_12_11_99_taxonomy.qza \
--o-classifier UNITE_12_11_99_classifier.qza

#this release was messing up, so I downloaded an already trained classifir from unite here (UNITE version 7.2): https://forum.qiime2.org/t/unite-ver-7-2-2017-12-01-classifiers-for-qiime2-ver-2017-12-available-here/3020
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

#barcode is on the forward read (see word doc showing how the R1 has the reverse primers in the read, while the R2 has the forward primers, something got switched up. Or possibly, the initial primers got removed, but the sequence is so short that it kept reading the other half of the primer. However the mapping file is correct: "reverse primer" matches to the R2). looking at the barcodes and the barcode file, it looks like it is the reverse complement
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
#220 and 210 respectively are when the median quality scores hit 25, however reading the documentation it looks like reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step. Also, I don't think it matters as much for euks b/c the read is so short, if the primer is successfully removed, there is no issue with quality. reads are ~130bp, primer is 24 or 16bp, so the primers are well within the good quality read.
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

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat filename|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 244bp, and blasts to a fungus

#trying to grep a prokaryote sequence, it doesn't find it
grep "TCGATAAAAATGATTGGCGTATCCAACCTGCAGAGTTTTATCGCTTCCATG" dna-sequences.fasta


#export table to get otu table
qiime tools export \
table.qza \
--output-dir exported-table

#convert biom to otu table text file!
biom convert -i feature-table.biom -o otu_table.txt --to-tsv

#delete the space and # from '#OTU ID' on line 2 
#sed '2s/.......//' otu_table.txt > otu_table2.txt
#sed '2s/.//' otu_table.txt > otu_table2.txt
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

#generat tree from masked alignment
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
#start 8:09 pm, 8:15pm
# -2 njobs means all but 1 CPU is used
qiime feature-classifier classify-sklearn \
--i-classifier 18S_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy.qza

#second time with all taxa in database
#start 10:48, end 10:55
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifierRL150.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy2.qza

#third time with all taxa in database
#start 2:55, end 2:59
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier_majority.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy3.qza

#fourth time with blast, there are a lot of defaults that I accepted here because I couldn't figure out what they were in the last 97% analysis. These defaults might actually be wrong b/c thre is one about percent identity whose default is 80%. There is a post saying blast is terrible and will calssify anything, no matter how bad the fit is: https://github.com/benjjneb/dada2/issues/323
#I don't think I'll use blast or uclust, I feel like they work with OTU clusterign methods, not with amplicon sequence variant analysis like DADA2.
#start 11:42 pm, end 11:54
qiime feature-classifier classify-consensus-blast \
--i-query rep-seqs.qza \
--i-reference-reads all_SILVA128_99_otus.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy_plain.qza \
--o-classification taxonomy4.qza

#Fifth time with the non truncated classifier
#start 3:39, end 3:47
##USE THIS
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy5.qza

#Sixth time with the 111 release
#start 3:25, end 3:36
##USE THIS
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA111_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy6.qza



#explanation of sklearn and training 
#https://forum.qiime2.org/t/classifier-training-questions/1162/3

qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

qiime tools export \
taxonomy2.qza \
--output-dir exported-taxonomy2

qiime tools export \
taxonomy3.qza \
--output-dir exported-taxonomy3

qiime tools export \
taxonomy4.qza \
--output-dir exported-taxonomy4

qiime tools export \
taxonomy5.qza \
--output-dir exported-taxonomy5

qiime tools export \
taxonomy6.qza \
--output-dir exported-taxonomy6

#navigate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy2.qza \
--o-visualization taxonomy2.qzv

qiime metadata tabulate \
--m-input-file taxonomy3.qza \
--o-visualization taxonomy3.qzv

qiime metadata tabulate \
--m-input-file taxonomy4.qza \
--o-visualization taxonomy4.qzv

qiime metadata tabulate \
--m-input-file taxonomy5.qza \
--o-visualization taxonomy5.qzv

qiime metadata tabulate \
--m-input-file taxonomy6.qza \
--o-visualization taxonomy6.qzv

qiime tools view taxonomy6.qzv


#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy2.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots2.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy3.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots3.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy4.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots4.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy5.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots5.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy6.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots6.qzv

qiime tools view taxa-bar-plots2.qzv
qiime tools view taxa-bar-plots5.qzv
#you can download a csv file from this visualization, might be useful

#notes: there are unassigned taxa in the file that was blasted against the all database, there are no unassigned taxa in the file that blasted agains the Euk only database. this might be because when you blast against euks it just leavs out anything that doesnot have a hit (b/c they are also likely bacteria so it doesn't make sense to call them unassigned ??)
#looking only at eukaroyte classified reads, many of the numbers in the groups (level2) are identical, a few are a little different, with the Euk only database always having more reads (if they are not the same). the euk dataset also has many many more reads for the Eukaryota;__ classification. That is where the main difference is. This is because it looks like any read that is unclassified or classified as bacteria in the dataset (from the all database) is calssified as a Eukaryota;__ in the euk only database. (i.e. the total number of reads is the same). Looking at some of the taxonomic groups at Level 2, some of the numbers of reads are higher in the all data base, some higher in the euk-only database. I can't explain that, except it is a different training so there is variability
#comparing concensus and majority at level 2 they are nearly identical, does not affect the number of reads in Eukaryota;__
#comparing truncated vs not, they look almost identical
#comparing release 128 to 111, 111 has fewer unassigned, about the same number of Eukaryota;__



##### Bact #####
##### Demultiplexing #####

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/Bact

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip reverse.fastq 
gzip forward.fastq

#start 
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences \
--output-path emp-paired-end-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv

#barcode is on the forward read, it is not the reverse complement
#start 10:25pm, end 12:20am
qiime demux emp-paired \
--m-barcodes-file 515BC_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux 

#summarize, start 10:50am, end 11:24am
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#export it so you can look inside at the files
#started 10:51, end 10:54
qiime tools export \
demux.qza \
--output-dir exported-demux

#The amplicon size should be 390 bp (according to the earth microbiome website). So theoretically there shouldn't be primers in the reads, but the first one I checked did have primers in it. So I need to take them out.
library(Biostrings)
reverseComplement(DNAString("GGACTACNVGGGTWTCTAAT"))

###### Trim primers/adapters ######
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 9:28pm, error 10:07
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f ATTAGAWACCCBNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose

#Error:
Traceback (most recent call last):
  File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/plugin/model/file_format.py", line 24, in validate
self._validate_(level)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/q2_types/per_sample_sequences/_format.py", line 159, in _validate_
self._check_n_records(record_count_map[level])
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/q2_types/per_sample_sequences/_format.py", line 130, in _check_n_records
% (i * 4 + 1))
qiime2.plugin.model.base.ValidationError: Missing sequence for record beginning on line 17

The above exception was the direct cause of the following exception:
  
  Traceback (most recent call last):
  File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/q2cli/commands.py", line 224, in __call__
results = action(**arguments)
File "<decorator-gen-358>", line 2, in trim_paired
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/sdk/action.py", line 228, in bound_callable
output_types, provenance)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/sdk/action.py", line 391, in _callable_executor_
spec.qiime_type, output_view, spec.view_type, prov)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/sdk/result.py", line 239, in _from_view
result = transformation(view)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/core/transform.py", line 57, in transformation
self.validate(view)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/core/transform.py", line 131, in validate
view.validate('min')
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/plugin/model/directory_format.py", line 171, in validate
getattr(self, field)._validate_members(collected_paths, level)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/plugin/model/directory_format.py", line 101, in _validate_members
self.format(path, mode='r').validate(level)
File "/Users/farrer/miniconda3/envs/qiime2-2017.12/lib/python3.5/site-packages/qiime2/plugin/model/file_format.py", line 29, in validate
) from e
qiime2.plugin.model.base.ValidationError: /var/folders/d0/9mbd394w8xq_ch0059bpzfs80000gn/T/q2-CasavaOneEightSingleLanePerSampleDirFmt-n8n35_9r/N.47.2015_48_L001_R1_001.fastq.gz is not a(n) FastqGzFormat file:
  
  Missing sequence for record beginning on line 17

Plugin error from cutadapt:
  
  /var/folders/d0/9mbd394w8xq_ch0059bpzfs80000gn/T/q2-CasavaOneEightSingleLanePerSampleDirFmt-n8n35_9r/N.47.2015_48_L001_R1_001.fastq.gz is not a(n) FastqGzFormat file:
  
  Missing sequence for record beginning on line 17

See above for debug info.

#I can't find any existing help on this problem, a kindof related problm suggested that the could be a bug when the filtering parameters are so stringent that the fasta file is empty. I can't see how this could be the case with trimming, so I could try going back to the mapping file and deleting that sample b/c it is an N sample that I'm not interested in anyway
#note that for adapter 1 it said it was preceded by a "G" extrememly frequently, however, I think it is fine.

#Exploring Ns
grep --color -n "^N"  S.0.2015_5_L001_R1_001copy.fastq # ^ means at the beginning of the line. no lines start with N
grep --color -n "^N"  S.0.2015_5_L001_R2_001copy.fastq #only like 5 lines start with N. so I won't worry about it.

#going back to the mapping file and deleting N.47.2015
##### Demultiplexing #####
#barcode is on the forward read, it is not the reverse complement
#start 3:04pm, end 4:59
qiime demux emp-paired \
--m-barcodes-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
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
#start 5:45pm, end 6:20pm, ran fine, no error
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f ATTAGAWACCCBNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
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
#R1 almost all got trimmed to 253 bp. R2 some got trimmed and some didn't probably b/c there are lots of errors at the end of the read and it did not recognize the primer.

#Getting a histogram of line lengths
awk '{print length}' S.0.2015_5_L001_R1_001copy.fastq | sort | uniq -c
#read 1 80036 at 253, 572 at 252, 84 at 251, 2 at 250

awk '{print length}' S.99.2015_101_L001_R1_001copy.fastq | sort | uniq -c
#read 1 96184 at 253, 1066 at 252, 224 at 251, 34 at 250

awk '{print length}' S.0.2015_5_L001_R2_001copy.fastq | sort | uniq -c
#read 2 56518 at 253, 300 at 252, 54 at 251, 28 at 242

awk '{print length}' S.99.2015_101_L001_R2_001copy.fastq | sort | uniq -c
#read 2 65670 at 253, 550 at 252, 142 at 251, 24 at 250

##### Denoising with DADA2 #####
#I need to figure out if i need to trim things
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 275bp and 225 respectively are when the median quality scores hit 25 - based on previous r script. using visualizaiton above (after primer trimming - the cutoffs for quality 25 are 252 and 233 ish), however reading the documentation the reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step.
#the reads I'm seeing are mostly 253bp, if the primers are 20 and 19 bp, then for R1, the primer should be within the high quality region. but for R2 the primer will be in the poor quality region so it may not have been removed in the above step. thus trimming at 253 should take the primer off even if there were lots of bp errors in it. However, I will truncate at 251 to delete the potentially bad quality bp without losing too many small fragments
#For R2, bad reads start at 233, so I should probably truncate there 
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#start 8:41pm, end 7:25am
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 251 \
--p-trunc-len-r 233

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it
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
#start 10:51am, end 11:29am. --p-n-threads -1 means use all cores
qiime alignment mafft \
--i-sequences rep-seqs.qza \
--p-n-threads -1 \
--o-alignment aligned-rep-seqs.qza

#Mask (or filter) the alignment to remove positions that are highly variable. These positions are generally considered to add noise to a resulting phylogenetic tree.
#start 12:00pm, end 12:30
qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

#generate tree from masked alignment
#start 12:31, end 1:03
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

#using silva
#my computer stalled after about 8 hrs with p-n-jobs -2 (which should mean all but one core), trying using fewer cores
#start 8:55pm, end 7:59am
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy.qza

qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

#using greengenes
#start 3:45pm, end 4:03pm

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
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime tools view taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy_gg.qza \
--o-visualization taxonomy_gg.qzv

qiime tools view taxonomy_gg.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--o-visualization taxa-bar-plots.qzv

qiime tools view taxa-bar-plots.qzv
#nice! there are very few unassigned and very few assigned to bacteria;__

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_gg.qza \
--m-metadata-file 515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.tsv \
--o-visualization taxa-bar-plots_gg.qzv

qiime tools view taxa-bar-plots_gg.qzv
#I like greengenes, there are fewer bacteria;__ and unassigned compared to silva




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
#the first tim I did this, I got an error:
Plugin error from demux:
  Mismatched sequence ids: M01918:229:000000000-ALYY5:1:2102:16350:8485 and ACAGAGGAGAG7:CAACTGEGGG<E*CBT8GGGGGGGGTAGGTCCG,:@EGGGDEB;TGCGFFGGGGGEGGGGGGGGG4:00-ALYY5:1:2102:1GGGGGGG18:229:0000C*CGCTGGCTGACG3EFTAAAGGGFGF?,CB+C8*GGGGG0CGGGGGGCGFGGGGGFGGGGGGGCACCTATCCTTGCGCAGGGGGCGCACCTGDGGGGGFGF?,GGGGGGGGGGGG7GF+CGTGG25:1:2TGTAGCGGTGGAATG8>5/CCFBGGFGGGGGGGGGGG2/CCFBG::::GGGGFGGGGGGGGGAAGTGGGGGGGGGGG)7GGGGGGGGGG2
#then I realized that I couldn't gunzip the sequences file (it was corrupted or something) so I re-created the gzip file from the raw .fastq file and then it ran fine
qiime demux emp-single \
--m-barcodes-file 515BC_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-column BarcodeSequence \
--i-seqs emp-single-end-forward-sequences.qza \
--o-per-sample-sequences demux 

#summarize, start 11:44pm, end 11:54pm
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

qiime tools view demux.qzv

#export it so you can look inside at the files
#started , end 
qiime tools export \
demux.qza \
--output-dir exported-demux

#The amplicon size should be 390 bp (according to the earth microbiome website). So theoretically there shouldn't be primers in the reads, but the first one I checked did have primers in it. So I need to take them out.
library(Biostrings)
reverseComplement(DNAString("GGACTACNVGGGTWTCTAAT"))

###### Trim primers/adapters ######
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 11:54pm, error at 12:10ish
qiime cutadapt trim-single \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter ATTAGAWACCCBNGTAGTCC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose

#Error: the same error happend that I got with the paired end data, so I will try with removing N.47.2015

#going back to the mapping file and deleting N.47.2015
##### Demultiplexing #####
#barcode is on the forward read, it is not the reverse complement
#start 12:19pm, end 1:24 
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
#R1 almost all got trimmed to 253 bp. R2 some got trimmed and some didn't probably b/c there are lots of errors at the end of the read and it did not recognize the primer.

#Getting a histogram of line lengths. at least for S.99.2015, the "histogram" below is exactly the same as for paired reads
awk '{print length}' S.0.2015_5_L001_R1_001copy.fastq | sort | uniq -c
#read 1 80036 at 253, 572 at 252, 84 at 251, 2 at 250

awk '{print length}' S.99.2015_101_L001_R1_001copy.fastq | sort | uniq -c
#read 1 96184 at 253, 1066 at 252, 224 at 251, 34 at 250, 3862 at 301


##### Denoising with DADA2 #####
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 275bp is when the median quality scores hit 25 - based on previous r script. using visualizaiton above (after primer trimming - the cutoffs for quality 25 is 254), however reading the documentation the reads that are shorter than the truncation numbers are discarded
#the reads I'm seeing are mostly 253bp, if the primers are 20 and 19 bp, then for R1, the primer should be within the high quality region, so should be taken off but there were some reads still at 301bp. thus trimming at 253 should take the primer off even if there were lots of bp errors in it. So before I truncated at 251 to delete the potentially bad quality bp without losing too many small fragments
#used to be trunc-len 251, 
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

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it
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




##### ITS #####
##### Demultiplexing #####

#Import the data into an artifact
#Navigate to: /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figures\&Stats/kingdata/QIIME2/ITS

#first rename your files to (R1) forward.fastq, (R2) reverse.fastq, (I) and barcodes.fastq
#then gzip them, takes a couple minutes
gzip barcodes.fastq
gzip reverse.fastq 
gzip forward.fastq

#import paired end reads
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences \
--output-path emp-paired-end-sequences.qza

#move mapping file into directory, and delete .txt and replace with .tsv
grep "GCATCGATGAAGAACGCAGC" Undetermined_S0_L001_R1_001.fastq
grep "GTGTAGATCTCGGTGGTCGCCGTATCATT" Undetermined_S0_L001_R2_001.fastq
grep "TTACTTCCTCTAAATGACCAAG" Undetermined_S0_L001_R2_001.fastq
grep "CGCAAATTCGAC" Undetermined_S0_L001_R1_001.fastq

#it is the reverse complement of the barcode in the index file
reverseComplement(DNAString("GTCGAATTTGCG"))
grep "CGCAAATTCGAC" Undetermined_S0_L001_I1_001.fastq

#barcode is on the forward read (see word doc showing how the R1 has the reverse primers in the read, while the R2 has the forward primers. The initial primers got removed, but the sequence is so short that it kept reading the other half of the primer. However the mapping file is correct: "reverse primer" matches to the R2). looking at the barcodes and the barcode file, it looks like it is the reverse complement
#start 9:06pm, end 12:00am
qiime demux emp-paired \
--m-barcodes-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-category BarcodeSequence \
--i-seqs emp-paired-end-sequences.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize, start 7:53am, end 8:10am
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
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 9:11am, end 9:30am
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f GCATCGATGAAGAACGCAGC \
--p-adapter-r TTACTTCCTCTAAATGACCAAG \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose
#warning that on or mor of your adaptors may be incomplete, R2 adaptor is preceded by T extremely often. I'm leaving it, since I'm pretty sure we have the adapter sequence correct.

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs
#R1 almost all got trimmed to 253 bp. R2 some got trimmed and some didn't probably b/c there are lots of errors at the end of the read and it did not recognize the primer.

#Getting a histogram of line lengths
awk '{print length}' S.0.2015_53_L001_R1_001_copy.fastq | sort | uniq -c
#read 1 lots of variability, but most around 230 (8000 reads) and 301 (10808 reads)

awk '{print length}' S.0.2015_53_L001_R2_001_copy.fastq | sort | uniq -c
#read 2 lots of variability, 8000 at 251, 1606 at 301

##### Denoising with DADA2 #####
#I need to figure out if i need to trim things. if I wanted to use fastx_trimmer, I would have to see if it would work. Whe I used it before, I did it on the initial file from the sequencing company, prior to demultiplexing (a fastq file). I mgiht be able to do it on multiple files at a later stage (after demultiplexing), but hten I'd have to figure out how to do it on many files (each sample individually) and then convert those file back to a qza file. Not impossible but maybe not worth it.
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 293bp and 230 respectively are when the median quality scores hit 25 - based on previous r script. using visualizaiton above, however the reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step.
#the reads I'm seeing are 230-250bp, if the primers are 20 and 22 bp, then for R1, the primer should be within the high quality region. but for R2 the primer will be in the poor quality region so it may not have been removed in the above step. however, I don't think there is anything I can do without losing a bunch of small reads. If I trim, I lose small reads that are deleted (with truncation) and large reads that no longer overlap prior to trimming, If I don't trim, I lose all reads where the primer was not removed
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#start 12:15pm, end it was probably 1.5 hrs into it that I noticed that there was an error
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0

#There was an error:
Plugin error from dada2:
  
  An error was encountered while running DADA2 in R (return code -11), please inspect stdout and stderr to learn more.

Debug info has been saved to /var/folders/d0/9mbd394w8xq_ch0059bpzfs80000gn/T/qiime2-q2cli-err-at5itkms.log

#I reinstalled biocLite and dada2 and tried it again
CDPATH= R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("dada2")'
#that didn't fix it, same error

#then I went to /Users/farrer/ and created a .Rprofile file with below in it. before there was no .Rprofile file at all, see https://forum.qiime2.org/t/dada2-exit-code-11-and-rprofile/1386/14

echo ".libPaths(.libPaths()[2])" > $HOME/.Rprofile

#then tried again, start 2:40pm, same error. 

#then I tried deleting sample N.103.2015 b/c also within the error was a note that there were no samples passing filter

#then I found that they released a patch: https://forum.qiime2.org/t/qiime-2-2017-12-release-is-now-live/2308/10?u=thermokarst
#I just needed to uninstall then reinstall qiime2.
#denoising again
#start: 3:24pm, end 5:10pm. NICE No errors!

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat dna-sequences.fasta|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 565bp, and blasts to nothing, must be an error of some kind

#lots are 301 which is suspicious, othrs are 230-250bp
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
#start 12:41pm, end 2:26pm
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy.qza

#start 8:50pm, end 9:07pm
qiime feature-classifier classify-sklearn \
--i-classifier unite-ver7-99-classifier-01.12.2017.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy_unite.qza


qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

qiime tools export \
taxonomy_unite.qza \
--output-dir exported-taxonomy_unite

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy_unite.qza \
--o-visualization taxonomy_unite.qzv

qiime tools view taxonomy.qzv
qiime tools view taxonomy_unite.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_unite.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots_unite.qzv

qiime tools view taxa-bar-plots.qzv
qiime tools view taxa-bar-plots_unite.qzv








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
#I need to figure out if i need to trim things. if I wanted to use fastx_trimmer, I would have to see if it would work. Whe I used it before, I did it on the initial file from the sequencing company, prior to demultiplexing (a fastq file). I mgiht be able to do it on multiple files at a later stage (after demultiplexing), but then I'd have to figure out how to do it on many files (each sample individually) and then convert those file back to a qza file. Not impossible but maybe not worth it.
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 290bp is when the median quality scores hit 25 - based on previous r script. using visualizaiton above, however the reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step.
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
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
#start 8:42pm, end 1:22 am
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier1.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy.qza

#start 11:47pm, end 12:13 pm 
qiime feature-classifier classify-sklearn \
--i-classifier unite-ver7-99-classifier-01.12.2017.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy_unite.qza

qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

qiime tools export \
taxonomy_unite.qza \
--output-dir exported-taxonomy_unite

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy_unite.qza \
--o-visualization taxonomy_unite.qzv

qiime tools view taxonomy.qzv
qiime tools view taxonomy_unite.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime tools view taxa-bar-plots.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_unite.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots_unite.qzv

qiime tools view taxa-bar-plots_unite.qzv




###### ITS - Processing only the reverse reads #####

#importing only the reverse read 
#the sequences need to be called "sequences.fastq.gz"
qiime tools import \
--type EMPSingleEndSequences \
--input-path emp-single-end-reverse-sequences \
--output-path emp-single-end-reverse-sequences.qza

###### Demultiplexing single reads ######
#start 1:24pm, end 2:17
qiime demux emp-single \
--m-barcodes-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-column BarcodeSequence \
--i-seqs emp-single-end-reverse-sequences.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize
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
--p-adapter TTACTTCCTCTAAATGACCAAG \
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
#lots of variability, but most around 234 (10,000 reads) and 301 (24834 reads)

##### Denoising with DADA2 #####
#I need to figure out if i need to trim things. if I wanted to use fastx_trimmer, I would have to see if it would work. Whe I used it before, I did it on the initial file from the sequencing company, prior to demultiplexing (a fastq file). I mgiht be able to do it on multiple files at a later stage (after demultiplexing), but then I'd have to figure out how to do it on many files (each sample individually) and then convert those file back to a qza file. Not impossible but maybe not worth it.
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 229bp is when the median quality scores hit 25 - based on previous r script. using visualizaiton above, however the reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step.
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#start 3:12pm, end 4:13 pm
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

#most are between 233-251 bp, but a lot at 301 too
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

#start 6:03pm, end 6:17 pm 
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







###### ITS - trimming reads first ####
#look at quality
qualityr <- read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/fastaqual/Undetermined_S0_L001_R2_001quality_by_cycle.txt", header=T, sep='\t')

plot(qualityr$mean ~ qualityr$column, type="l", xlab="Cycle number", ylab="mean quality score")
abline(h=25, col="red")

plot(qualityr$med ~ qualityr$column, type="l", xlab="Cycle number", ylab="mean quality score")
abline(h=25, col="red")

qualityr$column[which(qualityr$med<25)]
qualityr$column[which(qualityr$mean<25)]

qualityf <- read.table("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/fastaqual/Undetermined_S0_L001_R1_001quality_by_cycle.txt", header=T, sep='\t')

plot(qualityf$mean ~ qualityf$column, type="l", xlab="Cycle number", ylab="mean quality score")
abline(h=25, col="red")

plot(qualityf$med ~ qualityf$column, type="l", xlab="Cycle number", ylab="mean quality score")
abline(h=25, col="red")

qualityf$column[which(qualityf$med<25)]
qualityf$column[which(qualityf$mean<25)]

gunzip forward.fastq.gz
#start 6:52, end 7:15
/Users/farrer/Dropbox/EmilyComputerBackup/Desktop/bin/fastx_trimmer -Q33 -l 292 -i forward.fastq -o forwardt.fastq
gzip forwardt.fastq

gunzip reverse.fastq.gz
#start 7:23, end 7:46
/Users/farrer/Dropbox/EmilyComputerBackup/Desktop/bin/fastx_trimmer -Q33 -l 241 -i reverse.fastq -o reverset.fastq
gzip reverset.fastq

#Move zipped files to new folder and rename them

#import paired end reads
qiime tools import \
--type EMPPairedEndSequences \
--input-path emp-paired-end-sequences-trimmed \
--output-path emp-paired-end-sequences-trimmed.qza

#barcode is on the forward read (see word doc showing how the R1 has the reverse primers in the read, while the R2 has the forward primers. The initial primers got removed, but the sequence is so short that it kept reading the other half of the primer. However the mapping file is correct: "reverse primer" matches to the R2). looking at the barcodes and the barcode file, it looks like it is the reverse complement
#start 9:41pm, end 10:46pm
qiime demux emp-paired \
--m-barcodes-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--m-barcodes-column BarcodeSequence \
--i-seqs emp-paired-end-sequences-trimmed.qza \
--o-per-sample-sequences demux \
--p-rev-comp-mapping-barcodes

#summarize
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
#To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 5:49pm, end pm
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f GCATCGATGAAGAACGCAGC \
--p-adapter-r TTACTTCCTCTAAATGACCAAG \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose
#warning that on or mor of your adaptors may be incomplete, R2 adaptor is preceded by T extremely often. I'm leaving it, since I'm pretty sure we have the adapter sequence correct.

#error in N.129.2015_104_L001_R1_001.fastq.gz, M00517:307:000000000-AM49E:1:2119:18214:8101 1:N:0:0 saying the read and quality length don't match. that doesn't make sense, when I grep it the read and quality socres are both 292bp
grep -C 5 "M00517:307:000000000-AM49E:1:2119:18214:8101 1:N:0:0" N.129.2015_104_L001_R1_001_copy.fastq
if i want to I could delete that sample and rerun everything...
XXXXXXXXX

#summarize
qiime demux summarize \
--i-data trimmed-seqs.qza \
--o-visualization trimmed-seqs.qzv

qiime tools view trimmed-seqs.qzv

#export
qiime tools export \
trimmed-seqs.qza \
--output-dir exported-trimmed-seqs
#R1 almost all got trimmed to 253 bp. R2 some got trimmed and some didn't probably b/c there are lots of errors at the end of the read and it did not recognize the primer.

#Getting a histogram of line lengths
awk '{print length}' S.0.2015_53_L001_R1_001_copy.fastq | sort | uniq -c
#read 1 lots of variability, but most around 230 (8000 reads) and 301 (10808 reads)

awk '{print length}' S.0.2015_53_L001_R2_001_copy.fastq | sort | uniq -c
#read 2 lots of variability, 8000 at 251, 1606 at 301

##### Denoising with DADA2 #####
#I need to figure out if i need to trim things. if I wanted to use fastx_trimmer, I would have to see if it would work. Whe I used it before, I did it on the initial file from the sequencing company, prior to demultiplexing (a fastq file). I mgiht be able to do it on multiple files at a later stage (after demultiplexing), but hten I'd have to figure out how to do it on many files (each sample individually) and then convert those file back to a qza file. Not impossible but maybe not worth it.
#DADA2 webpage https://benjjneb.github.io/dada2/tutorial.html suggests not truncating if you're using ITS data b/c of the large variability in read lengths. Since DADA2 does take into account bp quality it should be fine. however if you can truncate it will increase sensitivity to rare variants.
#at 293bp and 230 respectively are when the median quality scores hit 25 - based on previous r script. using visualizaiton above, however the reads that are shorter than the truncation numbers are discarded, thus I should keep everything or do truncation prior to this step.
#the reads I'm seeing are 230-250bp, if the primers are 20 and 22 bp, then for R1, the primer should be within the high quality region. but for R2 the primer will be in the poor quality region so it may not have been removed in the above step. however, I don't think there is anything I can do without losing a bunch of small reads. If I trim, I lose small reads that are deleted (with truncation) and large reads that no longer overlap prior to trimming, If I don't trim, I lose all reads where the primer was not removed
#note for n-threads specifying 0 means use all cores
#note for trunc-len specifying 0 means don't truncate
#start 12:15pm, end it was probably 1.5 hrs into it that I noticed that there was an error
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed-seqs.qza \
--o-table table \
--o-representative-sequences rep-seqs \
--p-n-threads 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0

#There was an error:
Plugin error from dada2:
  
  An error was encountered while running DADA2 in R (return code -11), please inspect stdout and stderr to learn more.

Debug info has been saved to /var/folders/d0/9mbd394w8xq_ch0059bpzfs80000gn/T/qiime2-q2cli-err-at5itkms.log

#I reinstalled biocLite and dada2 and tried it again
CDPATH= R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("dada2")'
#that didn't fix it, same error

#then I went to /Users/farrer/ and created a .Rprofile file with below in it. before there was no .Rprofile file at all, see https://forum.qiime2.org/t/dada2-exit-code-11-and-rprofile/1386/14

echo ".libPaths(.libPaths()[2])" > $HOME/.Rprofile

#then tried again, start 2:40pm, same error. 

#then I tried deleting sample N.103.2015 b/c also within the error was a note that there were no samples passing filter

#then I found that they released a patch: https://forum.qiime2.org/t/qiime-2-2017-12-release-is-now-live/2308/10?u=thermokarst
#I just needed to uninstall then reinstall qiime2.
#denoising again
#start: 3:24pm, end 5:10pm. NICE No errors!

#export rep seqs just to take a look at. especially to look at how the long sequences got identified, I don't think I can look at this, since there isn't a file that has the ID DNA header, the sequence, and the feature ID in it
qiime tools export \
rep-seqs.qza \
--output-dir exported-rep-seqs

#to print the line with the longest number of characters:
cat dna-sequences.fasta|awk '{print length, $0}'|sort -nr|head -1
#the longest line in the rep set is 565bp, and blasts to nothing, must be an error of some kind

#lots are 301 which is suspicious, othrs are 230-250bp
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
#start 12:41pm, end 2:26pm
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy.qza

#start 8:50pm, end 9:07pm
qiime feature-classifier classify-sklearn \
--i-classifier unite-ver7-99-classifier-01.12.2017.qza \
--i-reads rep-seqs.qza \
--p-n-jobs 2 \
--o-classification taxonomy_unite.qza


qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

qiime tools export \
taxonomy_unite.qza \
--output-dir exported-taxonomy_unite

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy_unite.qza \
--o-visualization taxonomy_unite.qzv

qiime tools view taxonomy.qzv
qiime tools view taxonomy_unite.qzv

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_unite.qza \
--m-metadata-file ITS_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots_unite.qzv

qiime tools view taxa-bar-plots.qzv
qiime tools view taxa-bar-plots_unite.qzv

