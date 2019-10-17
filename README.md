# auto_metaviromics
An automated pipeline working on metagenomics analyses of mosquito virome  
The current version: 0.2.2  
Contact: Sohta Ishikawa Ph.D (sota.ishikawa@pasteur.fr)  

***
## Flowchart of the pipeline  
![](./image/auto_metaviromics.svg)  

## Usage
**NOTE:the pipeline assumes that it works on a cluster computer managed with Slurm, and dependencies required by pipeline are loaded from the Environment Modules package.**  
**Otherwise you need to install the below dependencies by yourself and set them in your PATH**  
+ Dependencies
    + [SeqKit](https://github.com/shenwei356/seqkit)
    + [Trimmomatic](https://github.com/timflutre/trimmomatic)
    + [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    + [samtools](http://www.htslib.org/doc/samtools.html)
    + [Megahit](https://github.com/voutcn/megahit)
    + [metaSPAdes](http://cab.spbu.ru/software/spades/)
    + [R](https://www.r-project.org/)
    + [VirFinder](https://github.com/jessieren/VirFinder)
    + [GeneMark.hmm with Heuristic models](http://exon.gatech.edu/index.html)
        + Packaged in MetaGeneMark software
    + [EMBOSS](http://www.bioinformatics.nl/emboss-explorer)
    + [vcftools](https://vcftools.github.io/index.html)  
    
```
cd ./pipeline
find ./ -type f -name "*.tar.gz" -exec tar zxf {} \;
# follow the tutorial in [NCBItax2lin](https://github.com/zyxue/ncbitax2lin) to dump the latest NCBI taxonomy list to 'lineages-YEAR-MM-DD.csv'

# to run the pipeline on the Slurm management
sbatch run_pipeline.sh
```  
The job batch file require the below options;
+ Options: 
  +	**--CPU**=(integer)                 specify number of CPUs to be used
  +	**--MEM**=(integer)                 specify total amount of MEMORY to be used (e.g. specify 16 if you want to use 16GB)
  +	**--READ-TYPE**=(integer)           specify 0 if your dataset is the single-end reads and 1 for the pair-end reads
  +	**--SAMPLE-PATH**=(path)            specify path to the sample file to be analyzed, e.g. /path/to/sample/XXX if the pair-end read sequence files are named as XXX_R1.fastq.gz (single-end files should be names as XXX.fastq.gz) 
  + **--MOSQ-NAME**=(string)            specify the full name of your mosquito sampled
  + **--DATE**=(DD-MON-YEAR)            specify the sampling date
  + **--LOCATION**=(string)             specify the sampling location
  +	**--FILTERING**=(boolean)           specify yes or no if you want/don't want to filter host & rRNA reads from your raw data
  +	**--HOST-DB**=(path)                specify path to BOWTIE2 database of the host mosquito reference sequences
  +	**--RNA-DB**=(path)                 specify path to BOWTIE2 database of the mosquito & Human rRNA reference sequences
  +	**--ASSEMBLING**=(boolean)          specify yes or no if you want/don't want to assemble your sample
  + **--ASSEMBLER**=(string)            specify name of the assembler to be used, megahit or metaSPAdes
  +	**--ASSEMBLE-OUTPUT**=(path)        specify path to the output directory of previously performed assemble, if you specified --ASSEMBLING=no
  + **--VIRAL-CONTIG-SEARCH**=(boolen)  specify yes or no if you want/don't want to search possible viral contigs by VirFInder
  + **--VIRFINDER-MODEL**=(path)        specify path to the user defined VirFinder model, if you specified --VIRAL-CONTIG-SEARCH=yes
  + **--VIRAL-CONTIG-LIST**=(path)      specify path to the user defined viral contig list if you specified --VIRAL-CONTIG-SEARCH=no and want to subject your list to the read mapping
  +	**--READ-MAPPING**=(boolean)        specify yes or no if you want/don't want to map reads on your contigs
  + **--ORF-PRED**=(boolen)             specify yes or no if you want/don't want to predict ORFs on your contigs and perform their functional annotation via DIAMOND
  +	**--DIAMOND-DB**=(path)             specify path to DIAMOND formatted database which must be assigned with NCBI taxonomy, if you specified --ORF-PRED=yes
  +	**--LINEAGE**=(path)                specify path to NCBI taxonomy dump file converted in lineage format, if you specified --ORF-PRED=yes
+ The pipeline provide some checkpoints by turn on/off the --FILTERING, --ASSEMBLING, --VIRAL-CONTIG-SEARCH, --READ-MAPPING and --ORF-PRED options. You can stop and restart the pipeline at any of them.
+ References of the pipileine: 
  + [Zakrzewski et al. 2018](https://www.nature.com/articles/s41598-018-22945-y)
  + [Belda et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6702732)

## Procedures
### Screening low-quality, rRNA- and host-derived reads
+ Remove low-quality reads by **Trimmomatic**
+ Remove host & human rRNA reads by mapping the low-quality filtered reads on the reference using **BOWTIE2** with '--sensitive-local' option
    + The reference of 5.8S, 28S, 18S, 16S, 12S rRNA was built by collecting corresponding sequences of mosquito and human from GENBANK in a FASTA format. 
+ Remove host mosquito reads by mapping the low-quality and rRNA filtered reads on the reference
    + The reference was build from NGS scaffolds from the [VectorBase](https://www.vectorbase.org/organisms/aedes-aegypti) in a FASTA format.
+ Total number of filtered reads can be computed at each step of screening by [SeqKit](https://github.com/shenwei356/seqkit)  

### Assembling contigs and identifying viral contigs
+ Assemble contigs from the non-low quality and non-rRNA, non-mosquito reads set using **Megahit** or **metaSPAdes**
+ Subject assembled contigs to the viral sequence prediction by their k-mer frequency using **VirFinder**
    + VirFinder uses an user-trained model to score the possibility of each contig as a viral sequence
    + The pipeline embedded model was trained based on the CDS sequences of global taxonomy of RNA virus
        + In GENBANK, all available CDS data was retrieved from Riboviria, of which genome was completed.
        + Then the collected sequences were clustered based on their nucleotide identity (>95%) using CD-HIT
        + Clustered sequences set was subjected to the model training by VirFinder, following the [instruction](https://github.com/jessieren/VirFinder)
        + CDS sequences of host mosquito species were also collected and clustered in the same way and used in the model training.  

### Mapping reads on the selected contigs themselves
+ For each contig selected from BlastX hits
    +  Map non-low quarity and non-mosquito reads to the contig itself using **BOWTIE2** with relaxed parameters of the similarity and alignment length
    +  a preset parameters of '--very-sensitive-local' option is applied
    +  Create a consensus sequence from the mapped reads
+ Compute read count, coverage depth (read per position), and call variants on the consensus sequence.  

### ORF prediction and protein gene annotation
+ Contigs from the above mapping procedure is then subjected to the gene prediction by [GeneMark.hmm](http://exon.gatech.edu/index.html). 
+ Predicted genes (ORFs) are translated into protein sequences and then subjected to BlastP search on NCBI NR database using **DIAMOND**
+ ~~DIAMOND output file (.DAA) is loaded in **MEGAN6** to perform gene ontology annotation via InterPro2GO and visualize functional annotation of the viral contigs in interest~~
+ BlastP result is also summarized in a tab-delimited table.
+ Results of ORF prediction and functional annotation on the contigs of interest are summarized and visualized by [UGENE](http://ugene.net/)

TO DO
***

### Phylogenetic analysis of RdRp protein sequences
+ Extract translated protein sequences which are annotated with the RNA-dependent RNA polymerase
+ Align them with RdRp protein sequences of known ssRNA (positive and negative) and dsRNA viruses exhensively collected from the [ViPR](https://www.viprbrc.org/brc/home.spg?decorator=vipr) database
+ Sequences are aligned automatically using **MAFFT v7.407** and then manually checked to trim ambiguously aligned positions.
+ Resultant MSAs are subjected to the maximum-likelihood phylogenetic analysis IQ-TREE with the best model selection option (-m TEST).

### Genome comparison 
+ Nucleotide sequences of RNA viral contigs found in samples are translated and compared with each other using tblastX approach (Belda et al. 2019)
+ Comparison of genomic similarity between viral contigs is visualized using **[Easyfig v2.2.2](https://mjsull.github.io/Easyfig/)**

## What we obtain from this pipeline?
+ Detection of RNA viral full/partial genomes, with taxonomic annotations, ORF predictions, comparison of gene repertory/order/similarity, and variant callings
+ Sequencing quality check on a target viral genome, e.g. rRNA read counts to chek the success of the rRNA depletion 
+ Viral taxonomic diversity and abundance in mosquito sample(s) of interest.
+ Detection of possible 'new' viral sequences
