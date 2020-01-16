#!/bin/bash
#SBATCH -J META_PIPELINE
#SBATCH -o Logs/log.pipeline.MosqMadaMascaren_2_S10.txt
#SBATCH -e Logs/err.pipeline.MosqMadaMascaren_2_S10.txt
#SBATCH -c 4
#SBATCH --mem=30000
#SBATCH -p geva
#SBATCH --qos=geva

source /pasteur/homes/sishikaw/.bashrc
source /local/gensoft2/adm/etc/profile.d/modules.sh
module load Trimmomatic bowtie2/2.3.4.3 samtools/1.9 bedtools/2.25.0 diamond/0.9.24 Python/3.7.2 fasta htslib/1.9 primer3/1.1.4 ClustalW golden tabix/0.2.6 perl/5.30.1 aragorn/1.2.38 blast+/2.9.0 infernal/1.1.2 barrnap/0.7 prodigal/2.6.3 ncbitools/20170106 minced/0.2.0 signalp/4.1
module load mafft/7.407 megahit/1.1.2 EMBOSS R/3.6.0 vcftools/0.1.13 SPAdes/3.12.0 hmmer/3.2.1 prokka/1.14.5

srun perl metagenomics_pipeline_v0.2.7.pl \
--CPU=2 \
--MEM=30 \
--READ-TYPE=1 \
--SAMPLE-PATH=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/SEQ_DATA_Anopheles_Culex/Madagascar_140419/MosqMadaMascaren_2_S10 \
--MOSQ-NAME=Anopheles_mascarensis \
--DATE=23-OCT-2019 \
--LOCATION=Madagascar \
--FILTERING=yes \
--HOST-DB=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/DB_RNA_Mosquito/Mosquito_Anopheles_scaf \
--RNA-DB=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/DB_RNA_Mosquito/Anopheles_Human_rRNA \
--ASSEMBLING=yes \
--ASSEMBLER=metaSPAdes \
--ASSEMBLE-OUTPUT=../metaSPAdes/spades_output_MosqMadaMascaren_2_S10/scaffolds.fasta \
--SUM-DIAMOND=yes \
--DIAMOND-DB=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/DB_DIAMOND/U-RVDBv16.0-prot_tax_091919 \
--VIRAL-CONTIG-SEARCH=no \
--VIRFINDER-MODEL=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/Metagenomics_Pipeline/CODEs/VirFinder/pipeline/VF.trainModUser.AnFunes_Riboviria_CDS.rda \
--VIRAL-CONTIG-LIST=../VirFinder/viral_contigs_list.txt \
--READ-MAPPING=no \
--HUMMER=no \
--HMMER-DB=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/DB_HMMER/U-RVDBv16.0-prot.fasta \
--LINEAGE=/pasteur/projets/policy01/GEVAmetagenomics/users/Sota/Metagenomics_Pipeline/Madagascar/metaSPAdes/pipeline/lineages-2019-09-18.csv

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
#	--SUM-DIAMOND=(boolean)         specify yes or no if you want/don't want to identify possible viral sequences from your assembled contigs using DIAMOND+BLASTX
#	--DIAMOND-DB=(path)             specify path to DIAMOND formatted database which must be assigned with NCBI taxonomy
#   --VIRAL-CONTIG-SEARCH=(boolen)  specify yes or no if you want/don't want to search possible viral contigs by VirFInder
#   --VIRFINDER-MODEL=(path)        specify path to the user defined VirFinder model, if you specified --VIRAL-CONTIG-SEARCH=yes
#   --VIRAL-CONTIG-LIST=(path)      specify path to the user defined viral contig list if you specified --VIRAL-CONTIG-SEARCH=no and want to subject your list to the read mapping
#	--READ-MAPPING=(boolean)        specify yes or no if you want/don't want to map reads on your contigs
#   --HUMMER=(boolen)               specify yes or no if you want/don't want to perform functional annotation of predicted ORFs by HUMMER
#	--HUMMER-DB=(path)              specify path to HUMMER formatted database if you specified --HUMMER=yes
#	--LINEAGE=(path)                specify path to NCBI taxonomy dump file converted in lineage format
