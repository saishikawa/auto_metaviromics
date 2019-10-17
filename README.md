# auto_metaviromics
An automated pipeline working on metagenomics analyses of mosquito virome
The current version: 0.2.2
Contact: Sohta Ishikawa Ph.D (sota.ishikawa@pasteur.fr)

***
## Flowchart of the pipeline  
![](./image/auto_metaviromics.svg)  

## Usage
**NOTE:the pipeline assumes that it works on a cluster computer managed with Slurm, and dependencies required by pipeline are loaded from the Environment Modules package.  Otherwise you need to install the below dependencies by yourself**
```
cd ./pipeline
# follow the tutorial in [NCBItax2lin](https://github.com/zyxue/ncbitax2lin) to dump the latest NCBI taxonomy list to 'lineages-XXXX-XX-XX.csv'

# to run the pipeline on TARS cluster
sbatch run_pipeline.sh
```  
The job batch file require the below options;
+ Options: 
    +	--CPU=(integer)             specify number of CPUs to be used
    +	--MEM=(integer)             specify total amount of MEMORY to be used (e.g. specify 16 if you want to use 16GB)
    +	--SINGLE-PAIR=(integer)     specify 0 if your dataset is the single-end reads and 1 for the pair-end reads
    +	--SAMPLE-NAME=(string)      specify the name of sample to be analyzed
    +	--FILTERING=(boolean)       specify yes or no if you want/don't want to filter host & rRNA reads from your raw data
    +	--HOST-DB=(path)            specify path to BOWTIE2 database of the host mosquito reference sequences
    +	--RNA-DB=(path)             specify path to BOWTIE2 database of the mosquito & Human rRNA reference sequences
    +	--ASSEMBLING=(boolean)      specify yes or no if you want/don't want to assemble your sample
    +	--ASSEMBLE-OUTPUT=(path)    specify path to the output directory of previously performed assemble, if you specified --ASSEMBLING=no
    +	--SUM-DIAMOND=(boolean)     specify yes or no if you want/don't want to run diamond analysis and its summarization on the host-removed reads set
    +	--DIAMOND-RVDB=(path)       specify path to DIAMOND formatted RVDB protein database which must be assigned with NCBI taxonomy
    +	--DIAMOND-NR=(path)         specify path to DIAMOND formatted NCBI nr protein database which must be assigned with NCBI taxonomy
    +	--LINEAGE=(path)            specify path to NCBI taxonomy dump file converted in lineage format
    +	--MAPPING=(boolean)         specify yes or no if you want/don't want to map reads on your contigs
    +	--CONTIGS=(path)            specify path to the contigs set to be subjected to the two-steps read mapping if --MAPPING=yes
+ The pipeline provide some checkpoints by turn on/off the --FILTERING, --ASSEMBLING, --SUM-DIAMOND, and --MAPPING options. You can stop and restart the pipeline at any points.
+ Reference of Pipileine: 
    + [Zakrzewski et al. 2018](https://www.nature.com/articles/s41598-018-22945-y)
    + [Belda et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/31429704)

## Procedures
### Screening RNAseq reads
+ Remove low-quality reads by **Trimmomatic v0.36**
+ Remove host & human rRNA reads by mapping the low-quality filtered reads on the reference using **BOWTIE2 v2.3.4.3** with '--sensitive-local' option
    + The reference of 5.8S, 28S, 18S, 16S, 12S rRNA was built by collecting corresponding sequences of mosquito and human in NCBI database. 
+ Remove host mosquito reads by  by mapping the low-quality and rRNA filtered reads on the reference
    + The reference was build from NGS scaffold sets from the [VectorBase](https://www.vectorbase.org/organisms/aedes-aegypti)
    + [Aedes aegypti cell line 2](https://www.vectorbase.org/organisms/aedes-aegypti) and 14 Anopheles species NGS scaffolds were used to screen the host-derived reads.
+ Total number of filtered reads can be computed at each step of screening by **[SeqKit](https://github.com/shenwei356/seqkit)**   

### Contig assemble, BlastX search adn taxonomic annotation
+ Assemble contigs from the non-low quality and non-mosquito reads set using **MEGAHIT v1.1.2**
+ Subject assembled contigs to the blastX similarity search on [RVDB](https://rvdb-prot.pasteur.fr/) protein database using **DIAMOND v0.9.24**
+ At present v16.0 of RVDB database assigned with the latest NCBI taxonomy (September 2019) is used.
+ Based on its blastx information (possiblly includes multple hits to different taxonomy), and using the LCA algorithm, each contig is annotated to a single taxonomic lineage unless it is 'unclassified'.

### Screening contigs
+ Detect RNA viral contigs from the DIAMOND result
    + Screen contigs by length > 1,000 bp, read coverage > 10, and blastX e-value > 1e-10.
    + Among remained contigs, based on their taxonomic annotations, pick possible RNA viral contigs.
    + Pick some interesting (possiblly viral) contigs to subject them to the 'two-steps read mapping'

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
