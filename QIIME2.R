##QIIME2

## first I needed to uninstall and then reinstall miniconda3. it would not update by itself with the update code
https://conda.io/docs/user-guide/install/macos.html

## how to install qiime2
https://docs.qiime2.org/2017.12/install/native/


#Activate Qiime env, must do this in every new terminal tab that is opened
source activate qiime2-2017.12 #for laptop
source activate qiime2-2018.2 #for desktop

#Paired end tutorial
https://docs.qiime2.org/2017.12/tutorials/atacama-soils/
#General tutorial
https://docs.qiime2.org/2017.12/tutorials/moving-pictures/


#discussion about whether it is worth it to join paired end reads when the reads completely overlap
https://forum.qiime2.org/t/question-about-dada2-denoise-paired-analysis/464
#antoher discussion of relaxing the maxEE filtering parameter when using full reads with poor quality scores near the ends (in qiime dada2 denoise-paired). however the default does delete everyting after the first instance of a bp with quality score of 2 (this can be changed as well)
https://github.com/qiime2/q2-dada2/issues/48






##### Taxonomic databases #####

#####Silva#####
#notes about the silva release: https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_notes.txt
#workflow from this thread: https://forum.qiime2.org/t/18s-classifier-using-silva-database-and-emb-primers/361

#OTUs
#first time
SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta
#second time
SILVA_128_QIIME_release/rep_set/rep_set_all/99/99_otus.fasta

#Taxonomy
#first time
SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_all_levels.txt
#second time
SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_all_levels.txt

#Import OTUs, try it with the 18S only first
#start 6:23pm, end 6:24
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_128_QIIME_release/rep_set/rep_set_18S_only/99/99_otus_18S.fasta \
--output-path 18S_SILVA128_99_otus.qza

#second time, with all taxa, start 7:33pm, end 7:36
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_128_QIIME_release/rep_set/rep_set_all/99/99_otus.fasta \
--output-path all_SILVA128_99_otus.qza


#Import taxonomy, again only with 18S. I'm using consensus taxonomy, which means all of the sequences need to have the same exact species taxonomy to be identified as a species, otherwise it is identified as "ambiguous taxa". This is very conservative. the alternative is using majority_taxonomy which 90% of the sequences in the database have to match. One comment in a forum is that the majority should be conservative enough for most uses, but if you are looking into a particular taxon, you should probably take the sequence and blast it yourself separately to make sure you are referencing the correct taxonomy. B/c I might do this when looking at hubs/connectors, I'll leave it a the highly conservative consensus, since I don't want to reblast things and I don't want to mis-report.
#Need to duplicate and save the txt file as tsv
#start 6:50, end 6:50
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/18S_only/99/consensus_taxonomy_7_levels.tsv \
--output-path 18S_SILVA128_99_taxonomy.qza

#second time, all taxa. this is odd, there is also a file called just "taxonomy_all_levels.txt" not sure how that is different from consnsus_taxonomy_all_levels.txt
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_all_levels.tsv \
--output-path all_SILVA128_99_taxonomy.qza

#third time, all taxa majority or I could also try the vsearch classifier as that is comparable to uclust https://forum.qiime2.org/t/import-i-reference-taxonomy-taxonomy-tsv-to-qza/781/14. i used blast though. or I could do a different train classifier thing below
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--source-format HeaderlessTSVTaxonomyFormat \
--input-path SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_all_levels.tsv \
--output-path all_SILVA128_99_taxonomy.qza

#Extract EMBP variable region.
#I'm not entirely sure what length to use, probably should use the median length from my euk data below (actually I can't easily extract this, it doesn't show up in the visualization summary I would need to figure out how to code this in unix)? or 150, since that is exactly what the person writing this workflow used. I'm actually not sure I need a trun-len if we give it the primers below...
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

#export to see what it trimmed. there are definitely squences that are 101 bp, thus lower than 150, so it's not like it deleted short sequences.
qiime tools export \
ref-seqs_18S_99_SILVA128_RL150.qza \
--output-dir exported-ref-seqs_18S_99_SILVA128_RL150

#Train classifier
#start 8:02, end 8:03
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_18S_99_SILVA128_RL150.qza \
--i-reference-taxonomy 18S_SILVA128_99_taxonomy.qza \
--o-classifier 18S_EMB_SILVA128_classifier.qza

#start 10:25, end 10:26
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs_all_99_SILVA128_RL150.qza \
--i-reference-taxonomy all_SILVA128_99_taxonomy.qza \
--o-classifier all_EMB_SILVA128_classifier.qza

#talking about what th different outputs (truncatd taxonomy vs blank genus/species levels means: 
https://forum.qiime2.org/t/consensus-blast-taxonomy-strings/586/2




#####Euk data#####

#Demultiplexing

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
  
#Trim primers/adapters. To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
#start 9:44pm, end 10:02
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 4 \
--p-adapter-f GTAGGTGAACCTGCAGAAGGATCA \
--p-adapter-r GACGGGCGGTGTGTAC \
--o-trimmed-sequences trimmed-seqs.qza \
--verbose
#note that it said that for the reverse read, the primer was preceeded by C extremely often, which might mean it is part of the primer. it shouldn't be though, so I'm leaving it in.
#not I could have used --p-error-rate 0 \  meaning no errors in the adapter, the default used above is 10%, not sure what difference it would make

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

# note that for some of the forward reads start with N for a basepair, this is actually "true" when you look at the reverse and forward reads, the N represents a certain bp that was apparntly unknown in the sequencing. hopfully that won't present a problem. so I won't try to delete this N now. (also when it starts with a g rather than an n (gctac) the g blasts to something so the g is correct). Update: the DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html states that no N are allowable in DADA2, so I need to get rid of that first basepair.

grep --color -n "^N"  N.0.2015_49_L001_R1_001copy.fastq # ^ means at the beginning of the line. there actually aren't that many, only like 30 reads start with N
grep --color -n "^N"  N.0.2015_49_L001_R2_001copy.fastq #not any Ns at the beginning of reads like in R1

#Denoising with DADA2 (within qiime2)
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


#Create a phylogenetic tree

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

#Assign taxonomy
#start 8:09 pm, 8:15pm
# -2 njobs means all but 1 CPU is used
qiime feature-classifier classify-sklearn \
--i-classifier 18S_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy.qza

#second time with all taxa in database
#stat 10:48, end 10:55
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy2.qza

qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

#second time
qiime tools export \
taxonomy2.qza \
--output-dir exported-taxonomy2

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy2.qza \
--o-visualization taxonomy2.qzv

qiime tools view taxonomy.qzv
qiime tools view taxonomy2.qzv
#this is not visualizing for some reason

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime tools view taxa-bar-plots.qzv







##### Bact #####
#Demultiplexing

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

#Trim primers/adapters. To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
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

#Trim primers/adapters. To get these adapters I looked in the raw sequnce data and found and made sure they were there. for example, the adapter-f ended up being the reverse complement of the reverse primer
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

#Denoising with DADA2 (within qiime2)
#I ned to figure out if i need to trim things
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


#Create a phylogenetic tree

#do the alignment
#start 10:51am, end 11:29am. --p-n-threads -1 means use all cors
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


#Assign taxonomy
#start 8:09 pm, 8:15pm
# -2 njobs means all but 1 CPU is used
qiime feature-classifier classify-sklearn \
--i-classifier 18S_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy.qza

#second time with all taxa in database
#stat 10:48, end 10:55
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy2.qza

qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

#second time
qiime tools export \
taxonomy2.qza \
--output-dir exported-taxonomy2

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy2.qza \
--o-visualization taxonomy2.qzv

qiime tools view taxonomy.qzv
qiime tools view taxonomy2.qzv
#this is not visualizing for some reason

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime tools view taxa-bar-plots.qzv


#Assign taxonomy
#start 8:09 pm, 8:15pm
# -2 njobs means all but 1 CPU is used
qiime feature-classifier classify-sklearn \
--i-classifier 18S_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy.qza

#second time with all taxa in database
#stat 10:48, end 10:55
qiime feature-classifier classify-sklearn \
--i-classifier all_EMB_SILVA128_classifier.qza \
--i-reads rep-seqs.qza \
--p-n-jobs -2 \
--o-classification taxonomy2.qza

qiime tools export \
taxonomy.qza \
--output-dir exported-taxonomy

#second time
qiime tools export \
taxonomy2.qza \
--output-dir exported-taxonomy2

#navegate to new directory, take out all spaces so it can be read into R (even the space in "unculutred eukaryote")
sed 's/[ ]//' taxonomy.tsv > taxonomy2.tsv
#then I still need to open the file in excel and save as a .csv - not sure why the import of the txt file is screwing up

#start 8:15, end 8:15
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime metadata tabulate \
--m-input-file taxonomy2.qza \
--o-visualization taxonomy2.qzv

qiime tools view taxonomy.qzv
qiime tools view taxonomy2.qzv
#this is not visualizing for some reason

#visualize barplots
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file EukBr_Niwot_20072015_All_MapFilenewlomehi.tsv \
--o-visualization taxa-bar-plots.qzv

qiime tools view taxa-bar-plots.qzv





