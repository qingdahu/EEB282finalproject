# Final project for DevBio282

Qingda Hu 
December 12, 2018


## Outline

RNA-seq data was obtained from publicly available resources and analyzed using a published pipeline. The goal of this project is for me to gain experience working with RNA-seq data which I will have to do for parts of my PhD research project. Pipeline was taken from 'https://www.nature.com/articles/nprot.2012.016' and data was taken from 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=DRP003328'




```
#make some directories for the project
cd /pub/$USER
mkdir finalproject
cd finalproject
mkdir {data,temp,output}
cd /pub/$USER/finalproject/data/
mkdir paperdata
cd paperdata

# get fastq files
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32038/suppl/GSE32038%5Fsimulated%5Ffastq%5Ffiles%2Etar%2Egz
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Drosophila_melanogaster/Ensembl/BDGP5.25/Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
tar -xvzf GSE32038_simulated_fastq_files.tar.gz
tar -xvzf Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz
#for some strange reason even with the -z option, the output was still .fq.gz
gunzip *.gz # it unzips the tar.gz file as well which if I open with tar will still generate fq.gz files.

#get the igenome required for the analysis
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Drosophila_melanogaster/Ensembl/BDGP5/Drosophila_melanogaster_Ensembl_BDGP5.tar.gz
tar -xvzf Drosophila_melanogaster_Ensembl_BDGP5.tar.gz 

# going to a multicpu interactive session will probably help from here on out
# qrsh -q bio -pe openmp 32
# I found it better to not use suspend-able nodes

#HPC already has bowtie and samtools
module load samtools/0.1.18 #the paper uses 0.1.17 but that is not available
module load bowtie/0.12.7
module load tophat/1.4.0 # the paper uses 1.3.2 but that is not available
module load cufflinks/2.2.1 # the paper uses 1.3.0

GENOME="/pub/$USER/finalproject/data/Drosophila_melanogaster/Ensembl/BDGP5/Sequence/BowtieIndex/genome"
GENEGTF="/pub/$USER/finalproject/data/Drosophila_melanogaster/Ensembl/BDGP5/Annotation/Genes/genes.gtf"
DATA="/pub/$USER/finalproject/data/paperdata"
TEMP="/pub/$USER/finalproject/temp"

tophat -p 32 -G $GENEGTF -o $TEMP/C1_R1_thout $GENOME $DATA/GSM794483_C1_R1_1.fq $DATA/GSM794483_C1_R1_2.fq
tophat -p 32 -G $GENEGTF -o $TEMP/C1_R2_thout $GENOME $DATA/GSM794484_C1_R2_1.fq $DATA/GSM794484_C1_R2_2.fq
tophat -p 32 -G $GENEGTF -o $TEMP/C1_R3_thout $GENOME $DATA/GSM794485_C1_R3_1.fq $DATA/GSM794485_C1_R3_2.fq
#each run was taking about 1hr 28mins so I tried to run the next 3 together to see check the speed difference
tophat -p 10 -G $GENEGTF -o $TEMP/C2_R1_thout $GENOME $DATA/GSM794486_C2_R1_1.fq $DATA/GSM794486_C2_R1_2.fq &
tophat -p 10 -G $GENEGTF -o $TEMP/C2_R2_thout $GENOME $DATA/GSM794487_C2_R2_1.fq $DATA/GSM794487_C2_R2_2.fq &
tophat -p 10 -G $GENEGTF -o $TEMP/C2_R3_thout $GENOME $DATA/GSM794488_C2_R3_1.fq $DATA/GSM794488_C2_R3_2.fq &


# in fact these were faster 1hr 09 mins each so more than 3 times faster.


cufflinks -p 8 -o $TEMP/C1_R1_clout $TEMP/C1_R1_thout/accepted_hits.bam &
cufflinks -p 8 -o $TEMP/C1_R2_clout $TEMP/C1_R2_thout/accepted_hits.bam &
cufflinks -p 8 -o $TEMP/C1_R3_clout $TEMP/C1_R3_thout/accepted_hits.bam &
cufflinks -p 8 -o $TEMP/C2_R1_clout $TEMP/C2_R1_thout/accepted_hits.bam &
cufflinks -p 8 -o $TEMP/C2_R2_clout $TEMP/C2_R2_thout/accepted_hits.bam &
cufflinks -p 8 -o $TEMP/C2_R3_clout $TEMP/C2_R3_thout/accepted_hits.bam &


touch $TEMP/assemblies.txt 
```

/pub/$USER/finalproject/temp/C1_R1_clout/transcripts.gtf
/pub/$USER/finalproject/temp/C2_R2_clout/transcripts.gtf
/pub/$USER/finalproject/temp/C1_R2_clout/transcripts.gtf
/pub/$USER/finalproject/temp/C2_R1_clout/transcripts.gtf
/pub/$USER/finalproject/temp/C1_R3_clout/transcripts.gtf
/pub/$USER/finalproject/temp/C2_R3_clout/transcripts.gtf


```
cd $TEMP
cuffmerge -g $GENEGTF -s $GENOME.fa -p 8 $TEMP/assemblies.txt



# Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate:
cuffdiff -o $TEMP/diff_out -b $GENOME.fa -p 8 -L C1,C2 -u $TEMP/merged_asm/merged.gtf \
$TEMP/C1_R1_thout/accepted_hits.bam,$TEMP/C1_R2_thout/accepted_hits.bam,$TEMP/C1_R3_thout/accepted_hits.bam \
$TEMP/C2_R1_thout/accepted_hits.bam,$TEMP/C2_R3_thout/accepted_hits.bam,$TEMP/C2_R2_thout/accepted_hits.bam

cp $TEMP/diff_out /pub/qingdah/finalproject/output/diff_out1

cd /pub/qingdah/finalproject/output

module load R
R
```

```
library(cummeRbund)
cuff_data <- readCufflinks('/pub/qingdah/finalproject/temp/diff_out')
jpeg('density1.jpg')
csDensity(genes(cuff_data))
dev.off()
jpeg('scatter1.jpg')
 csScatter(genes(cuff_data), 'C1', 'C2')
dev.off()
jpeg('volcano1.jpg')
csVolcano(genes(cuff_data), 'C1', 'C2')
dev.off()
jpeg('expression1.jpg')
mygene<- getGene(cuff_data,'regucalcin')
expressionBarplot(mygene)
dev.off()
jpeg('isoform1.jpg')
expressionBarplot(isoforms(mygene))
dev.off()

```

I am not sure why but my density and volcano plot look very different from the results in the paper. See folder output. I used a different igenome and different versions of the software which might have caused this or there was a mistake in my analysis. 







## Tentative procedure for DRP003328
NOTE: Currently running CUFFLINKS

Here is how I am processing data from 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=DRP003328'. I don't have enough time to do all the runs as I did not realize how big data sets can get. I will add the results when I get them but it will be after the deadline of this project. However, I think it would still be valuable for my own learning to continue this. 

Prepare assembly2.txt in advance

```
/pub/qingdah/finalproject/temp/DRR055272_thout/transcripts.gtf
/pub/qingdah/finalproject/temp/DRR055273_thout/transcripts.gtf
/pub/qingdah/finalproject/temp/DRR055274_thout/transcripts.gtf
/pub/qingdah/finalproject/temp/DRR055274_thout/transcripts.gtf
```


```
# obtain data from DRP003328
module load SRAToolKit/2.3.2-5
# for some reason fastq-dump doens't work with 'DRP003328' but for this project, I think it is ok to just work with 2 conditions each with 2 replicates. 2 and 3 are male adults of Drosophila pseudoobscura while 4 and 5 are female
#I did this part in tmux to avoid getting disconnected from SSH
cd /pub/$USER/finalproject/data/
fastq-dump -F --gzip --split-files DRR055272 2>/dev/null &
fastq-dump -F --gzip --split-files DRR055273 2>/dev/null &
fastq-dump -F --gzip --split-files DRR055274 2>/dev/null &	
fastq-dump -F --gzip --split-files DRR055275 2>/dev/null &


#get the igenome required for the analysis
#wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Drosophila_melanogaster/Ensembl/BDGP5/Drosophila_melanogaster_Ensembl_BDGP5.tar.gz
#tar -xvzf Drosophila_melanogaster_Ensembl_BDGP5.tar.gz 


# going to a multicpu interactive session will probably help from here on out
# qrsh -q bio -pe openmp 32
# I found it better to not use suspend-able nodes


#HPC already has bowtie and samtools
module load samtools/0.1.18 #the paper uses 0.1.17 but that is not available
module load bowtie/0.12.7
module load tophat/1.4.0 # the paper uses 1.3.2 but that is not available
module load cufflinks/2.2.1 # the paper uses 1.3.0

GENOME="/pub/$USER/finalproject/data/Drosophila_melanogaster/Ensembl/BDGP5/Sequence/BowtieIndex/genome"
GENEGTF="/pub/$USER/finalproject/data/Drosophila_melanogaster/Ensembl/BDGP5/Annotation/Genes/genes.gtf"
DATA="/pub/$USER/finalproject/data/compressed"
TEMP="/pub/$USER/finalproject/temp"

# Map the reads for each sample to the reference genome:
tophat -p 8 -G $GENEGTF -o $TEMP/DRR055272_thout $GENOME $DATA/DRR055272_1.fastq $DATA/DRR055272_2.fastq
tophat -p 8 -G $GENEGTF -o $TEMP/DRR055273_thout $GENOME $DATA/DRR055273_1.fastq $DATA/DRR055273_2.fastq
tophat -p 8 -G $GENEGTF -o $TEMP/DRR055274_thout $GENOME $DATA/DRR055274_1.fastq $DATA/DRR055274_2.fastq
tophat -p 8 -G $GENEGTF -o $TEMP/DRR055275_thout $GENOME $DATA/DRR055275_1.fastq $DATA/DRR055275_2.fastq
# Assemble transcripts for each sample:
cufflinks -p 8 -o $TEMP/DRR055272_chout $TEMP/DRR055272_thout/accepted_hits.bam
cufflinks -p 8 -o $TEMP/DRR055273_chout $TEMP/DRR055273_thout/accepted_hits.bam
cufflinks -p 8 -o $TEMP/DRR055274_chout $TEMP/DRR055274_thout/accepted_hits.bam
cufflinks -p 8 -o $TEMP/DRR055275_chout $TEMP/DRR055275_thout/accepted_hits.bam

touch TEMP/assemblies.txt 
vim TEMP/assemblies.txt #add following lines into file


# Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation:
cuffmerge -g $GENEGTF -s $GENOME.fa -p 8 $TEMP/assemblies.txt
# can i get this to write to specific location?

# Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate:
cuffdiff -o $TEMP/diff_out -b $GENOME.fa -p 8 -L C1,C2 -u merged_asm/merged.gtf \
$TEMP/DRR055272_thout/accepted_hits.bam,$TEMP/DRR055273_thout/accepted_hits.bam \
$TEMP/DRR055274_thout/accepted_hits.bam,$TEMP/DRR055275_thout/accepted_hits.bam

cp $TEMP/diff_out /pub/qingdah/finalproject/output/diff_out2
cd /pub/qingdah/finalproject/output
module load R
R
```

I switch to R to do the plotting.

```
library(cummeRbund)
cuff_data <- readCufflinks('/pub/qingdah/finalproject/temp/diff_out')
jpeg('density2.jpg')
csDensity(genes(cuff_data))
dev.off()
jpeg('scatter2.jpg')
csScatter(genes(cuff_data), 'female', 'male')
dev.off()
jpeg('volcano2.jpg')
csVolcano(genes(cuff_data), 'female', 'male')
dev.off()
jpeg('expression2.jpg')
mygene<- getGene(cuff_data,'regucalcin')
expressionBarplot(mygene)
dev.off()
jpeg('isoform2.jpg')
expressionBarplot(isoforms(mygene))
dev.off()
```





## Thoughts about this project

This was a very easy to follow set of instructions with the knowledge gained from devbio282. There were some useful tricks that were mentioned in class such as pushing tasks to background and piping stderr to dev/null. Experience with HPC modules made this pipeline much more streamlined. However, I was not use to handling the size of a real data set.



