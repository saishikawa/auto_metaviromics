#!/bin/bash
#SBATCH -J META_PIPELINE
#SBATCH -o Logs/log.pipeline.txt
#SBATCH -e Logs/err.pipeline.txt
#SBATCH -c 4
#SBATCH --mem=30000
#SBATCH -p geva
#SBATCH --qos=geva

#set PATH
source /pasteur/homes/sishikaw/.bashrc

source /local/gensoft2/adm/etc/profile.d/modules.sh
module load Trimmomatic bowtie2/2.3.4.3 samtools/1.9 diamond/0.9.24 Python/3.7.2 htslib/1.9 tabix/0.2.6 primer3/1.1.4 ClustalW golden
module load megahit/1.1.2 EMBOSS R/3.6.0 vcftools/0.1.13 SPAdes/3.12.0

srun perl metagenomics_pipeline_v0.2.2.pl \
--CPU=4 \
--MEM=30 \
--READ-TYPE=1 \
--SAMPLE-PATH=/path/to/your/sample \
--MOSQ-NAME=Aedes_aegypti \
--DATE=01-Apr-2019 \
--LOCATION=Madagascar \
--FILTERING=yes \
--HOST-DB=/path/to/DB/of/Host/sequences \
--RNA-DB=/path/to/DB/of/rRNA/sequences \
--ASSEMBLING=yes \
--ASSEMBLER=metaSPAdes \
--ASSEMBLE-OUTPUT=/path/to/assembled/contigs/file \
--VIRAL-CONTIG-SEARCH=yes \
--VIRFINDER-MODEL=../../../pipeline//VF.trainModUser.Riboviria_AnFunes.rda \
--VIRAL-CONTIG-LIST=/path/to/viral/contig/list \
--READ-MAPPING=yes \
--ORF-PRED=yes \
--DIAMOND-DB=/path/to/DIAMOND/DB/file \
--LINEAGE=/path/to/taxonomy/dump/file 

# Options: 
#	--CPU=(integer)                 specify number of CPUs to be used
#	--MEM=(integer)                 specify total amount of MEMORY to be used (e.g. specify 16 if you want to use 16GB)
#	--READ-TYPE=(integer)           specify 0 if your dataset is the single-end reads and 1 for the pair-end reads
#	--SAMPLE-PATH=(path)            specify path to the sample file to be analyzed, e.g. /path/to/sample/XXX if the pair-end read sequence files are named as XXX_R1.fastq.gz (single-end files should be names as XXX.fastq.gz) 
#   --MOSQ-NAME=(string)            specify the full name of your mosquito sampled
#   --DATE=(DD-MON-YEAR)            specify the sampling date
#   --LOCATION=(string)             specify the sampling location
#	--FILTERING=(boolean)           specify yes or no if you want/don't want to filter host & rRNA reads from your raw data
#	--HOST-DB=(path)                specify path to BOWTIE2 database of the host mosquito reference sequences
#	--RNA-DB=(path)                 specify path to BOWTIE2 database of the mosquito & Human rRNA reference sequences
#	--ASSEMBLING=(boolean)          specify yes or no if you want/don't want to assemble your sample
#   --ASSEMBLER=(string)            specify name of the assembler to be used, megahit or metaSPAdes
#	--ASSEMBLE-OUTPUT=(path)        specify path to the output directory of previously performed assemble, if you specified --ASSEMBLING=no
#   --VIRAL-CONTIG-SEARCH=(boolen)  specify yes or no if you want/don't want to search possible viral contigs by VirFInder
#   --VIRFINDER-MODEL=(path)        specify path to the user defined VirFinder model, if you specified --VIRAL-CONTIG-SEARCH=yes
#   --VIRAL-CONTIG-LIST=(path)      specify path to the user defined viral contig list if you specified --VIRAL-CONTIG-SEARCH=no and want to subject your list to the read mapping
#	--READ-MAPPING=(boolean)        specify yes or no if you want/don't want to map reads on your contigs
#   --ORF-PRED=(boolen)             specify yes or no if you want/don't want to predict ORFs on your contigs and perform their functional annotation via DIAMOND
#	--DIAMOND-DB=(path)             specify path to DIAMOND formatted database which must be assigned with NCBI taxonomy, if you specified --ORF-PRED=yes
#	--LINEAGE=(path)                specify path to NCBI taxonomy dump file converted in lineage format, if you specified --ORF-PRED=yes
