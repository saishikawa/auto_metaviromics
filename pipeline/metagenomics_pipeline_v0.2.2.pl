use strict;
use warnings;

my $CPU;
my $MEM;
my $read_type;
my $sample_name;
my $sample_name_short;
my $mosq_fullname;
my $date;
my $location;
my $filtering;
my $host_db;
my $rna_db;
my $assemble;
my $assembler;
my $assemble_output;
my $viral_contig_search;
my $virfinder_model;
my $vircontig_list;
my $mapping;
my $orf_prediction;
my $diamond_db;
my $lineage;

foreach my $argv (@ARGV){
    if ($argv =~ /CPU/) {
        $CPU = $argv;
        $CPU =~ s/--CPU=//;
    } elsif ($argv =~ /MEM/) {
        $MEM = $argv;
        $MEM =~ s/--MEM=//;
    } elsif ($argv =~ /READ-TYPE/) {
        $read_type = $argv;
        $read_type =~ s/--READ-TYPE=//;
    } elsif ($argv =~ /SAMPLE-PATH/) {
        $sample_name = $argv;
        $sample_name =~ s/--SAMPLE-PATH=//;
        $sample_name_short = $argv;
        $sample_name_short =~ s/[^\/]+\///g;
    } elsif ($argv =~ /MOSQ-NAME/) {
        $mosq_fullname = $argv;
        $mosq_fullname =~ s/--MOSQ-NAME=//;
    } elsif ($argv =~ /DATE/) {
        $date = $argv;
        $date =~ s/--DATE=//;
    } elsif ($argv =~ /LOCATION/) {
        $location = $argv;
        $location =~ s/--LOCATION=//;
    } elsif ($argv =~ /FILTERING/) {
        $filtering = $argv;
        $filtering =~ s/--FILTERING=//;
    } elsif ($argv =~ /HOST-DB/) {
        $host_db = $argv;
        $host_db =~ s/--HOST-DB=//;
    } elsif ($argv =~ /RNA-DB/) {
        $rna_db = $argv;
        $rna_db =~ s/--RNA-DB=//;
    } elsif ($argv =~ /ASSEMBLING/) {
        $assemble = $argv;
        $assemble =~ s/--ASSEMBLING=//;
    } elsif ($argv =~ /ASSEMBLER/) {
        $assembler = $argv;
        $assembler =~ s/--ASSEMBLER=//;
    } elsif ($argv =~ /ASSEMBLE-OUTPUT/) {
        $assemble_output = $argv;
        $assemble_output =~ s/--ASSEMBLE-OUTPUT=//;
    } elsif ($argv =~ /VIRAL-CONTIG-SEARCH/) {
        $viral_contig_search = $argv;
        $viral_contig_search =~ s/--VIRAL-CONTIG-SEARCH=//;
    } elsif ($argv =~ /VIRFINDER-MODEL/) {
        $virfinder_model = $argv;
        $virfinder_model =~ s/--VIRFINDER-MODEL=//;
    }elsif ($argv =~ /VIRAL-CONTIG-LIST/) {
        $vircontig_list = $argv;
        $vircontig_list =~ s/--VIRAL-CONTIG-LIST=//;
    } elsif ($argv =~ /READ-MAPPING/) {
        $mapping = $argv;
        $mapping =~ s/--READ-MAPPING=//;
    } elsif ($argv =~ /ORF-PRED/) {
        $orf_prediction = $argv;
        $orf_prediction =~ s/--ORF-PRED=//;
    } elsif ($argv =~ /DIAMOND-DB/) {
        $diamond_db = $argv;
        $diamond_db =~ s/--DIAMOND-DB=//;
    } elsif ($argv =~ /LINEAGE/) {
        $lineage = $argv;
        $lineage =~ s/--LINEAGE=//;
    }
}
main($CPU,$MEM,$read_type,$sample_name,$sample_name_short,$mosq_fullname,$date,$location,$filtering,$host_db,$rna_db,$assemble,$assembler,$assemble_output,$viral_contig_search,$virfinder_model,$vircontig_list,$mapping,$orf_prediction,$diamond_db,$lineage);


sub main {
    my ($CPU,$MEM,$read_type,$sample_name,$sample_name_short,$mosq_fullname,$date,$location,$filtering,$host_db,$rna_db,$assemble,$assembler,$assemble_output,$viral_contig_search,$virfinder_model,$vircontig_list,$mapping,$orf_prediction,$diamond_db,$lineage) = @_;

    print("\nStart with the below options\n");
    print("##################################################################################\n");
    print("$CPU CPUs and $MEM GB memory used\n");
    if ($read_type == 0) {
        print("READ TYPE: Single-End\n");
    } elsif($read_type == 1) {
        print("READ TYPE: Pair-End\n");
    }
    print("SAMPLE NAME : $sample_name_short\n");
    print("MOSQUIT NAME: $mosq_fullname\n");
    print("DO FILTERING\? : $filtering\n");
    print("HOST DB : $host_db\n");
    print("rRNA DB : $rna_db\n");
    print("DO ASSEMBLE \? : $assemble\n");
    print("ASSEMBLER: $assembler\n");
    print("ASSEMBLE OUTPUT : $assemble_output\n");
    print("DO VIRAL CONTIG SEARCH\? : $viral_contig_search\n");
    print("VIRFINDER MODEL : $virfinder_model\n");
    print("DO CONTIG MAPPING \? : $mapping\n");
    print("DO ORF PREDICTION AND FUNCTIONAL ANNOTATION\? : $orf_prediction\n");
    print("DIAMOND DB : $diamond_db\n");
    print("LINEAGE DB : $lineage\n");
    print("##################################################################################\n\n");

    mkdir("../analyzed_samples");
    chdir("../analyzed_samples");
    mkdir("$sample_name_short");
    chdir("$sample_name_short");
    my $vf_name;

    #Filtering and Assembling commands for single-end reads
    if ($read_type == 0) {
        if ($filtering =~ /yes/) {
            print("\n##################################################################################\n");
            print "Removing host reads and rRNA reads from raw data, then assemble remained reads\n";
            print("##################################################################################\n\n");

            #Computing statistics for raw data
            print "### Computing basic stats of the raw data\n";
            system("seqkit stats $sample_name.fastq.gz > ./Read_Stats.$sample_name_short.txt");

            mkdir("Remove_Host_RNA_Reads");
            chdir("Remove_Host_RNA_Reads");

            #Low-quality read filtering
            print "### Filtering low-quality reads\n";
            system("Trimmomatic SE -phred33 -threads $CPU $sample_name.fastq.gz ./$sample_name_short.fastq.qual.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10 > Trimmomatic.log");

            #rRNA reads removing via BOWTIE2
            print "### rRNA reads removing\n";
            system("bowtie2 --sensitive-local -p $CPU -x $rna_db -U ./$sample_name_short.fastq.qual.gz -S Reads.mapped_and_unmapped.SE.sam");
            system("samtools view -bS Reads.mapped_and_unmapped.SE.sam > Reads.mapped_and_unmapped.SE.bam");
            system("samtools view -b -f 4 Reads.mapped_and_unmapped.SE.bam > Reads.unmapped.SE.bam");
            system("samtools sort Reads.unmapped.SE.bam -O Reads.unmapped_sorted.SE.bam");
            system("bedtools bamtofastq -i Reads.unmapped_sorted.SE.bam -fq rRNA_removed.$sample_name_short.fastq");
            system("gzip rRNA_removed.$sample_name_short.fastq");
            system("rm *.sam *.bam");

            #Host reads removing via BOWTIE2
            print "### Host reads removing\n";
            system("bowtie2 --sensitive-local -p $CPU -x $host_db -U ./rRNA_removed.$sample_name_short.fastq.gz -S Reads.mapped_and_unmapped.SE.sam");
            system("samtools view -bS Reads.mapped_and_unmapped.SE.sam > Reads.mapped_and_unmapped.SE.bam");
            system("samtools view -b -f 4 Reads.mapped_and_unmapped.SE.bam > Reads.unmapped.SE.bam");
            system("samtools sort Reads.unmapped.SE.bam -O Reads.unmapped_sorted.SE.bam");
            system("bedtools bamtofastq -i Reads.unmapped_sorted.SE.bam -fq host_rRNA_removed.$sample_name_short.fastq");
            system("gzip host_rRNA_removed.$sample_name_short.fastq");
            system("rm *.sam *.bam");

            #Computing statistics for filtered reads
            print "### Computing read stas after filtering\n";
            system("seqkit stats rRNA_removed.$sample_name_short.fastq.gz host_rRNA_removed.$sample_name_short.fastq.gz >> ../Read_Stats.$sample_name_short.txt");

            chdir(".."); #CD=$sample_name
        }

        #Assembling
        if ($assemble =~ /yes/) {
            print("\n##################################################################################\n");
            print "Assembling contigs\n";
            print("##################################################################################\n\n");
            if($assembler =~ /megahit/){
                #Megahit assembling on the filtered reads
                mkdir("Megahit");
                chdir("Megahit");
                system("megahit --presets meta-sensitive -t $CPU -o megahit_output\_$sample_name_short -r ../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name_short.fastq.gz > Megahit.log");
                chdir("..");
                $assemble_output = "../Megahit/megahit_output\_$sample_name_short/final.contigs.fa";               
            }            
        } 
    } elsif ($read_type == 1) {
        if ($filtering =~ /yes/) {
            print("\n##################################################################################\n");
            print "Removing host reads from raw data\n";
            print("##################################################################################\n\n");

            #Computing statistics for raw data
            print "Computing basic stats of the raw data\n";
            system("seqkit stats $sample_name\_R1.fastq.gz $sample_name\_R2.fastq.gz > ./Read_Stats.$sample_name_short.txt");
            mkdir("Remove_Host_RNA_Reads");
            chdir("Remove_Host_RNA_Reads");
        
            #Low-quality read filtering
            print "Filtering low-quality reads\n";
            system("Trimmomatic PE -phred33 -threads $CPU $sample_name\_R1.fastq.gz $sample_name\_R2.fastq.gz ./$sample_name_short\_R1.fastq.qual.gz ./$sample_name_short\_R1.fastq.unpair.gz ./$sample_name_short\_R2.fastq.qual.gz ./$sample_name_short\_R2.fastq.unpair.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:10 > Trimmomatic.log");

            #rRNA reads removing via BOWTIE2
            print "rRNA reads removing\n";
            system("bowtie2 --sensitive-local -p $CPU -x $rna_db -1 ./$sample_name_short\_R1.fastq.qual.gz -2 ./$sample_name_short\_R2.fastq.qual.gz -S Reads.mapped_and_unmapped.PE.sam");
            system("samtools view -bS Reads.mapped_and_unmapped.PE.sam > Reads.mapped_and_unmapped.PE.bam");
            system("samtools view -b -f 12 -F 256 Reads.mapped_and_unmapped.PE.bam > Reads.bothEndsUnmapped.PE.bam");
            system("samtools sort Reads.bothEndsUnmapped.PE.bam -O Reads.bothEndsUnmapped_sorted.PE.bam");
            system("bedtools bamtofastq -i Reads.bothEndsUnmapped_sorted.PE.bam -fq rRNA_removed.$sample_name_short.R1.fastq -fq2 rRNA_removed.$sample_name_short.R2.fastq");
            system("gzip rRNA_removed.$sample_name_short.R1.fastq");
            system("gzip rRNA_removed.$sample_name_short.R2.fastq");
            system("rm *.sam *.bam");

            #Host reads removing via BOWTIE2
            print "Host reads removing\n";
            system("bowtie2 --sensitive-local -p $CPU -x $host_db -1 ./rRNA_removed.$sample_name_short.R1.fastq.gz -2 ./rRNA_removed.$sample_name_short.R2.fastq.gz -S Reads.mapped_and_unmapped.PE.sam");
            system("samtools view -bS Reads.mapped_and_unmapped.PE.sam > Reads.mapped_and_unmapped.PE.bam");
            system("samtools view -b -f 12 -F 256 Reads.mapped_and_unmapped.PE.bam > Reads.bothEndsUnmapped.PE.bam");
            system("samtools sort Reads.bothEndsUnmapped.PE.bam -O Reads.bothEndsUnmapped_sorted.PE.bam");
            system("bedtools bamtofastq -i Reads.bothEndsUnmapped_sorted.PE.bam -fq host_rRNA_removed.$sample_name_short.R1.fastq -fq2 host_rRNA_removed.$sample_name_short.R2.fastq");
            system("gzip host_rRNA_removed.$sample_name_short.R1.fastq");
            system("gzip host_rRNA_removed.$sample_name_short.R2.fastq");
            system("rm *.sam *.bam");

            #Computing statistics for filtered reads
            print "Computing read stas after filtering\n";
            system("seqkit stats rRNA_removed.$sample_name_short.R1.fastq.gz rRNA_removed.$sample_name_short.R2.fastq.gz host_rRNA_removed.$sample_name_short.R1.fastq.gz host_rRNA_removed.$sample_name_short.R2.fastq.gz >> ../Read_Stats.$sample_name_short.txt");
            chdir(".."); #CD=$sample_name
        }

        #Assembling
        if ($assemble =~ /yes/) {
            print("\n##################################################################################\n");
            print "Assembling contigs\n";
            print("##################################################################################\n\n");
            if($assembler =~ /megahit/){
                #Megahit assembling on the filtered reads
                mkdir("Megahit");
                chdir("Megahit");
                system("megahit --presets meta-sensitive -t $CPU -o megahit_output\_$sample_name_short -1 ../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name_short.R1.fastq.gz -2 ../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name_short.R2.fastq.gz > Megahit.log");
                chdir("..");
                $assemble_output = "../Megahit/megahit_output\_$sample_name_short/final.contigs.fa";
            } elsif($assembler =~ /metaSPAdes/){
                #metaSPAdes assembling on the filtered reads
                mkdir("metaSPAdes");
                chdir("metaSPAdes");
                system("metaspades.py -t $CPU -m $MEM --pe1-1 ../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name_short.R1.fastq.gz --pe1-2 ../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name_short.R2.fastq.gz -o spades_output\_$sample_name_short > metaSPAdes.log");
                chdir(".."); #CD=$sample_name
                $assemble_output = "../metaSPAdes/spades_output\_$sample_name_short/contigs.fasta";
            }
        }       
    }

    #DIAMOND analysis and summary
    if ($viral_contig_search =~ /yes/) {
        print("\n##################################################################################\n");
        print "Searching viral contigs by the k-mer frequency metric of VirFinder\n";
        print("##################################################################################\n\n");
        mkdir("VirFinder");
        chdir("VirFinder");
        system("grep \">\" $assemble_output > ./contigs.list.txt");
        system("sed -e \'s/\>//\' -i ./contigs.list.txt");

        #k-mer frequency based viral contig search via VirFinder
        open(NEWFILE, "> R_command_VirFinder.txt") or die "$!";
        print NEWFILE "library(VirFinder)\nmodFile <- \"$virfinder_model\"\nload(modFile)\ninFaFile <- \"$assemble_output\"\npredResultUser <- VF.pred.user(inFaFile, VF.trainModUser)\npredResultUser\$qvalue <- VF.qvalue(predResultUser\$pvalue)\nx <- predResultUser[order(predResultUser\$qvalue),]\nwrite.table(x, file=\"contigs.scored.csv\",quote=F,col.names=T)\n";
        close(NEWFILE);
        system("Rscript R_command_VirFinder.txt > VirFinder.log");
        $vf_name = "contigs.scored.csv";
        my $sum_rvdb = sum_VirFinder($vf_name,$lineage,$assembler);
        $vf_name = "../VirFinder/viral_contigs_list.txt";
        chdir("..");#CD=sample_name
    }
    if($viral_contig_search =~ /no/){
        $vf_name = $vircontig_list;
    }
    #Mapping on selected contigs
    if ($mapping =~ /yes/) {
        print("\n##################################################################################\n");
        print "Mapping reads on viral contigs\n";
        print("##################################################################################\n\n");
        mkdir("GENBANK");
        mkdir("Mapping_contigs");
        chdir("Mapping_contigs");
    
        #Complete contigs by the mapping reads on themselves
        my $chimera = chimera_contig($sample_name_short,$CPU,$vf_name,$assemble_output,$mosq_fullname,$date,$location,$orf_prediction,$diamond_db,$lineage,$assembler);
        chdir("..");#CD=#CD=sample_name
    }
    
    if($orf_prediction =~ /yes/){
        print("\n##################################################################################\n");
        print "Functional annotation of the predicted ORFs by BlastP (DIAMOND)\n";
        print("##################################################################################\n\n");
        system("diamond blastp -p $CPU -d $diamond_db -q ./PredictedProteins.$sample_name_short.fasta --top 10 --evalue 1.0e-10 -f 6 qseqid qlen stitle staxids slen length evalue pident qseq -o ./diamond.NR.PRO.$sample_name_short.tsv > DIAMOND_ORF.log");
        my $prot_name = "./diamond.NR.PRO.$sample_name_short.tsv";
        my $sum_contig = sum_diamond_prot($prot_name,$lineage,$sample_name_short);
        system("grep \">\" PredictedProteins.$sample_name_short.fasta > DIAMOND/ORFs_name.txt");
        chdir("DIAMOND");
        my $check_ORFs = check_ORFs($prot_name,$sample_name_short);
        chdir("..");
    }

    open(FILE, "./Read_Stats.$sample_name_short.txt") or die "$!";
    my @stats = <FILE>;
    close(FILE);
    my $line_num = 0;
    my $total_counts;
    my $rRNA_counts;
    my $host_counts;
    foreach my $stats (@stats){
        chomp $stats;
        my @list = split(/ +/, $stats);
        $list[3] =~ s/,//g;
        if($line_num == 1){
            $total_counts = $list[3];
        } elsif($line_num == 4) {
            $rRNA_counts = $list[3];
        } elsif($line_num == 6) {
            $host_counts = $list[3];
        }
        $line_num++;
    }
    my $tmp_rRNA_counts = $total_counts - $rRNA_counts;
    my $tmp_host_counts = $rRNA_counts - $host_counts;
    my $per_rRNA_counts =  ( ($total_counts - $rRNA_counts) / $total_counts) * 100;
    my $per_host_counts =  ( ($rRNA_counts - $host_counts) / $total_counts) * 100;
    my $stats_file = "../Read_Stats_Samples.tsv";
    if(-e $stats_file) {
        open(NEWFILE, ">> $stats_file") or die "$!";
        print NEWFILE "pair-end\t$sample_name_short\t$total_counts\t$tmp_rRNA_counts\t$tmp_host_counts\t$per_rRNA_counts\t$per_host_counts\n";
        close(FILE);
    } else {
        open(NEWFILE, ">> $stats_file") or die "$!";
        print NEWFILE "sequencing type\tsample ID\tTotal reads\trRNA reads\tHost reads\t% of rRNA reads\t% of host reads\n";
        print NEWFILE "pair-end\t$sample_name_short\t$total_counts\t$tmp_rRNA_counts\t$tmp_host_counts\t$per_rRNA_counts\t$per_host_counts\n";
        close(FILE);
    }

    chdir("..");#CD=/analyzed_samples/
}

# Summarizing DIAMOND results with taxonomy
sub sum_VirFinder {
    my ($filename,$lineage,$assembler) = @_;
    open(FILE, $filename) or die "$!";
    my @file = <FILE>;
    close(FILE);
    open(FILE, "contigs.list.txt") or die "$!";
    my @names = <FILE>;
    close(FILE);
    print("\n##################################################################################\n");
    print "Filtering VirFinder Hits: $filename...\n";
    print("##################################################################################\n\n");
    my @contig_ID;
    my @contig_cov;
    my @contig_length;
    my $count = 0;
    my @viral_contigs;
    my @contigs_list;

    foreach my $file (@file) {
        chomp $file;
        if($count == 0){
            push(@viral_contigs,"$file\n");
            $count ++;
            next;
        }
        if($assembler =~ /megahit/){
            my @list = split(/ /, $file);
            #contigID,coverage
            my $tmp_contigID = $list[1];
            my $tmp_contigcov = $list[3];
            $tmp_contigcov =~ s/multi=//;
            my $tmp_contiglen = $list[4];
            $tmp_contiglen =~ s/len=//;
            my $tmp_pval = $list[7];
            if($tmp_contiglen > 1000 and $tmp_contigcov > 10 and $tmp_pval < 0.01) {
                push(@viral_contigs,"$file\n");
                push(@contigs_list,"$tmp_contigID\n");
            }
            $count ++;
        } elsif($assembler =~ /metaSPAdes/){
            my @list = split(/ /, $file);
            #contigID,coverage
            my $tmp_contigID = $list[1];
            $tmp_contigID =~ s/\_length\_.+//;
            my $tmp_contigcov = $list[1];
            $tmp_contigcov =~ s/NODE[^c]+cov_//;
            my $tmp_contiglen = $list[2];
            my $tmp_pval = $list[3];
            if($tmp_contiglen > 1000 and $tmp_contigcov > 10 and $tmp_pval < 0.01) {
                push(@viral_contigs,"$file\n");
                push(@contigs_list,"$tmp_contigID\n");
            }
            $count ++;
        }
    }

    open(NEWFILE, "> selected_viral_contigs.csv") or die "$!";
    print NEWFILE @viral_contigs;
    close(NEWFILE);
    open(NEWFILE, "> viral_contigs_list.txt") or die "$!";
    print NEWFILE @contigs_list;
    close(NEWFILE);
    return 1;
}

#BlastP summary via DIAMOND
sub sum_diamond_prot {
    my ($filename,$lineage,$sample_name) = @_;
    open(FILE, $filename) or die "$!";
    my @file = <FILE>;
    close(FILE);
    print("\n##################################################################################\n");
    print "Filtering BlastP result: $filename...\n";
    print("##################################################################################\n\n");
    my @contig_ID;
    my @contig_length;
    my @subj_accession;
    my @subj_name;
    my @subj_length;
    my @align_length;
    my @evalue;
    my @pident;
    my @contig_orf;
    my @subj_lineage;
    my @URL_accession;
    my @URL_subj_lineage;
    my $count_kept_hit = 0;
#qseqid qlen qstart qend qframe stitle staxids slen length evalue pident qseq
    foreach my $file (@file) {
        chomp $file;
        my @list = split(/\t/, $file);
        #contigID,coverage
        my $tmp_contigID = $list[0];
        my $tmp_contiglen = $list[1];
        my $tmp_accession = $list[2];
        $tmp_accession =~ s/ .+//;
        my $tmp_name = $list[2];
        #$tmp_name =~ s/[ ]+ //;
        my $tmp_taxid = $list[3];
        my $tmp_subjlen = $list[4];
        my $tmp_alilen = $list[5];
        my $tmp_eval = $list[6];
        my $tmp_pident = $list[7];
        my $tmp_orf = $list[8];
        if($tmp_alilen > 100) {
            push(@contig_ID,$tmp_contigID);
            push(@contig_length,$tmp_contiglen);
            push(@subj_accession,$tmp_accession);
            push(@subj_name,$tmp_name);
            push(@subj_length,$tmp_subjlen);
            push(@align_length,$tmp_alilen);
            push(@evalue,$tmp_eval);
            push(@pident,$tmp_pident);
            push(@contig_orf,$tmp_orf);
            push(@URL_accession,"https\:\/\/www\.ncbi\.nlm\.nih\.gov\/protein\/$tmp_accession");
            push(@URL_subj_lineage,"https\:\/\/www\.ncbi\.nlm\.nih\.gov\/taxonomy\/\?term\=txid$tmp_taxid");
            #Extract lineage information
            open(FILE, $lineage) or die "$!";
            my @lineage_list = <FILE>;
            close(FILE);
            my $tmp_lineage = "empty\tempty";
            foreach my $lineage_list (@lineage_list){
                if($lineage_list =~ /^$tmp_taxid\,/){
                    chomp $lineage_list;
                    $lineage_list =~ s/^$tmp_taxid\,/$tmp_taxid\t/;
                    $lineage_list =~ s/\,{2,}/\|/g;
                    $tmp_lineage = $lineage_list;
                    last;
                }
            }
            push(@subj_lineage,$tmp_lineage);
            $count_kept_hit ++;
        }
    }
    mkdir("DIAMOND");
    system("mv $filename DIAMOND_ORF.log DIAMOND");
    open(NEWFILE, "> DIAMOND_BlastP_ORFs.$sample_name.tsv") or die "$!";
    print NEWFILE "Query_ID\tQuery_Length(aa)\tSubject_Name\tSubject_Len(aa)\tAlignment_Len(aa)\tE-value\tP-Ident\tSubject_taxonID\tSubject_Lineage\tURL_taxonomy\tURL_protein\tORF seq\n";
    for (my $i=0; $i<$count_kept_hit; $i++){
        print NEWFILE "$contig_ID[$i]\t$contig_length[$i]\t$subj_name[$i]\t$subj_length[$i]\t$align_length[$i]\t$evalue[$i]\t$pident[$i]\t$subj_lineage[$i]\t$URL_subj_lineage[$i]\t$URL_accession[$i]\t$contig_orf[$i]\n";
    }
    close(NEWFILE);
    return 1;
}

#Check viral ORFs with no siilarity in the BlastP search
sub check_ORFs {
    my ($filename,$sample_name) = @_;
    open(FILE, $filename) or die "$!";
    my @file = <FILE>;
    close(FILE);
    open(FILE, "ORFs_name.txt") or die "$!";
    my @ORFs = <FILE>;
    close(FILE);
    print("\n##################################################################################\n");
    print "Checking possible viral ORFs with no siilarity in the BlastP search: $filename...\n";
    print("##################################################################################\n\n");

    my @no_sim_ORFs;
    my $have_sim = 0;
    foreach my $ORFs (@ORFs) {
        $have_sim = 0;
        chomp $ORFs;
        $ORFs =~ s/\>//;
        foreach my $file (@file){
            chomp $file;
            if ($file =~ /$ORFs/){
                $have_sim = 1;
                last;
            }
        }
        if($have_sim == 0){
            push(@no_sim_ORFs,"$ORFs\n");
        }
    }
    open(NEWFILE, ">> ../ORFs_no_Similarity_BlastP.$sample_name.txt");
    print NEWFILE @no_sim_ORFs;
    close(NEWFILE);
    return 1;
}

#Two-step read mapping on selected contigs
sub chimera_contig {
    my ($sample_name,$CPU,$vf_name,$assemble_output,$mosq_fullname,$date,$location,$orf_prediction,$diamond_db,$lineage,$assembler) = @_;
    print("\n##################################################################################\n");
    print "Completeing viral contigs by mapping reads on themselves\n";
    print("##################################################################################\n\n");
    open(FILE,$vf_name) or die "$!";
    my @node_list = <FILE>;
    close(FILE);
    open(FILE,$assemble_output) or die "$!";
    my @seq = <FILE>;
    close(FILE);
    my $count_list = 0;

    foreach my $node_list (@node_list) {
        my @new_seq;
        chomp $node_list;
        my $count = 0;
        foreach my $seq (@seq){
            if($seq =~ />/){
                if($assembler =~ /megahit/){
                    if($seq =~ /$node_list flag/){
                        #my $tmp = $seq;
                        #chomp $tmp;
                        #$tmp =~ s/^[^\=]+\=[^\=]+\=[^\=]+\=//;
                        $seq =~ s/ /_/g;
                        push(@new_seq,$seq);
                        $count = 1;
                    } else {
                        $count = 0;
                    }
                } elsif($assembler =~ /metaSPAdes/){
                    if($seq =~ /$node_list\_length/){
                        #my $tmp = $seq;
                        #chomp $tmp;
                        #$tmp =~ s/^[^\=]+\=[^\=]+\=[^\=]+\=//;
                        push(@new_seq,$seq);
                        $count = 1;
                    } else {
                        $count = 0;
                    }
                }
            } else {
                if($count == 1){
                    push(@new_seq,$seq);
                }
            }
        }
        if(-d $node_list){
            next;
        } else {
            mkdir("$node_list");
            chdir("$node_list");
            open(NEWFILE,"> original_contig.fasta") or die "$!";
            print NEWFILE @new_seq;
            close(NEWFILE);

            system("cat original_contig.fasta >> ../Assembled_contigs.$sample_name.fasta");
            #First mapping
            system("bowtie2-build original_contig.fasta original_contig > log.index");

            if($read_type == 0){
                system("bowtie2 --very-sensitive-local -p $CPU -x original_contig -U ../../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name.fastq.gz -S second.reads.sam");
            } elsif($read_type == 1){
                system("bowtie2 --very-sensitive-local -p $CPU -x original_contig -1 ../../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name.R1.fastq.gz -2 ../../Remove_Host_RNA_Reads/host_rRNA_removed.$sample_name.R2.fastq.gz -S second.reads.sam");
            }

            system("samtools view -bS second.reads.sam > second.reads.bam");
            system("samtools sort second.reads.bam -O second.reads.sorted.bam");
            system("samtools index second.reads.sorted.bam");
            #Read_count
            system("samtools idxstats second.reads.sorted.bam > bam.stats.$sample_name.$node_list.txt");
            #Read_coverage
            system('samtools view -H second.reads.sorted.bam | perl -ne \'if ($_ =~ m/^\@SQ/){ print $_ }\' | perl -ne \'if ($_ =~ m/SN:(.+)\s+LN:(\d+)/) { print $1, "\t", $2, "\n"}\' > lengths.ref.txt');
            system("bedtools genomecov -ibam second.reads.sorted.bam -g lengths.ref.txt -d > mapped.reads.$sample_name.$node_list.cov");
            #Variant calling
            #system("bcftools mpileup -Ou -f original_contig.fasta second.reads.sorted.bam | bcftools call -mv -Oz -o calls.vcf");
            system("bcftools mpileup -Ou -f original_contig.fasta second.reads.sorted.bam | bcftools call -mv -Ov -o calls.vcf");
            system("bgzip calls.vcf");
            system("tabix -p vcf calls.vcf.gz");
            #system("bcftools index calls.vcf.gz");
            #Consensus
            #system("cat original_contig.fasta | bcftools consensus calls.vcf.gz > consensus.$sample_name.$node_list.fa");
            system("cat original_contig.fasta | vcf-consensus calls.vcf.gz > consensus.$sample_name.$node_list.fa");
            system("mv original_contig.fasta Original_contig.$sample_name.$node_list.fa");

            system("rm second.reads.bam *.sam");
            system("mv second.reads.sorted.bam Second.step.$sample_name.$node_list.sorted.bam");

            #Extract consensus
            system("sed -e \'s/>/>$sample_name./\' ./consensus.$sample_name.$node_list.fa > ./$node_list.consensus.second.$sample_name.fasta");
            open(FILE,"$node_list.consensus.second.$sample_name.fasta") or die "$!";
            my @cons = <FILE>;
            close(FILE);
            my $tmpname = $node_list;
            if($assembler =~ /megahit/){
                $tmpname =~ s/k\d+\_//;
            } elsif($assembler =~ /metaSPAdes/){
                $tmpname =~ s/NODE\_//;
            }
            $cons[0] = ">$sample_name\_contig\_$tmpname\n";
            open(NEWFILE,"> $node_list.consensus.second.$sample_name.fasta") or die "$!";
            print NEWFILE @cons;
            close(NEWFILE);
            system("cat ./$node_list.consensus.second.$sample_name.fasta >> ../../Completed_contigs.$sample_name.fasta");

            #Make summary of abundance (read counts) and coverage
            open(FILE, "bam.stats.$sample_name.$node_list.txt") or die "$!";
            my @info = <FILE>;
            close(FILE);
            my $count_info = 0;
            my $abundance = 0;
            my $total = 0;
            foreach my $info (@info) {
                chomp $info;
                if($count_info==0){
                    my @list = split(/\t/, $info);
                    $abundance = $list[2];
                } elsif ($count_info==1){
                    my @list = split(/\t/, $info);
                    $total = $list[3];
                }
                $count_info++;
            }
            open(FILE, "mapped.reads.$sample_name.$node_list.cov") or die "$!";
            @info = <FILE>;
            close(FILE);
            my $coverage = 0;
            my $positions = 0;
            foreach my $info (@info) {
                chomp $info;
                my @list = split(/\t/, $info);
                $coverage = $coverage + $list[2];
                $positions++;
            }
            my $avg_coverage = $coverage / $positions;

            #make GENBANK file
            if($orf_prediction =~ /yes/){
                system("gmhmmp -m source ../../../../pipeline/MetaGeneMark_v1.mod -o ProteinPredictionGFF.$sample_name.$node_list.txt -f 3 -a -A PredictedProteins.$sample_name.$node_list.fasta ./$node_list.consensus.second.$sample_name.fasta");
                open(FILE, "PredictedProteins.$sample_name.$node_list.fasta") or die "$!";
                my @proteins = <FILE>;
                close(FILE);
                my $linenum = 1;
                foreach my $proteins (@proteins){
                    if($proteins =~ />/){
                        $proteins =~ s/\>gene[^\>]+\>/\>ORF\_$linenum\|/;
                        $linenum++;
                    }
                    $proteins =~ s/^\n//;
                }
                open(NEWFILE, "> PredictedProteins.$sample_name.$node_list.fasta") or die "$!";
                print NEWFILE @proteins;
                close(NEWFILE);
                system("cat PredictedProteins.$sample_name.$node_list.fasta >> ../../PredictedProteins.$sample_name.fasta");
                system("seqret -sequence $node_list.consensus.second.$sample_name.fasta -feature -fformat gff -fopenfile ProteinPredictionGFF.$sample_name.$node_list.txt -osformat genbank -auto -outseq $sample_name.$node_list.gb");
                open(FILE, "$sample_name.$node_list.gb") or die "$!";
                my @genbank = <FILE>;
                close(FILE);
                my $i=1;
                my @new_genbank;
                foreach my $genbank (@genbank){
                    if($genbank =~ /LOCUS/){
                        $genbank =~ s/UNC [^\-]+\-[^\-]+\-\d\d\d\d/VRL $date/;
                    }
                    if($genbank =~ /FEATURES/){
                        $genbank =~ s/FEATURES\s+Location\/Qualifiers/FEATURES             Location\/Qualifiers\n     source          1..$positions\n                     \/organism=\"$location mosquito virus $sample_name.$tmpname\"\n                     \/strain=\"$sample_name\"\n                     \/host=\"$mosq_fullname\"\n                     \/country=\"$location\"/;
                    }
                    if($genbank =~ /mRNA/ or $genbank =~ /CDS/ or $genbank =~ /\/note/){
                        next;
                    }
                    if($genbank =~ /gene/){
                        $genbank =~ s/gene/CDS /;
                        my $tmp = $genbank;
                        my $tmp2 = "                     \/gene=\"ORF$i\"\n";
                        $genbank = $tmp.$tmp2;
                        $i++;
                    }
                    push(@new_genbank,$genbank);
                }
                open(NEWFILE, "> $sample_name.$node_list.mod.gb") or die "$!";
                print NEWFILE @new_genbank;
                close(NEWFILE);
                system("mv $sample_name.$node_list.mod.gb ../../GENBANK/$sample_name.$node_list.gb");

                #reverse_complement
                system("revseq ./$node_list.consensus.second.$sample_name.fasta ./$node_list.consensus.second.$sample_name.revcomp.fasta");
                system("gmhmmp -m source ../../../../pipeline/MetaGeneMark_v1.mod -o ProteinPredictionGFF.$sample_name.$node_list.revcomp.txt -f 3 -a -A PredictedProteins.$sample_name.$node_list.revcomp.fasta ./$node_list.consensus.second.$sample_name.revcomp.fasta");
                system("seqret -sequence $node_list.consensus.second.$sample_name.revcomp.fasta -feature -fformat gff -fopenfile ProteinPredictionGFF.$sample_name.$node_list.revcomp.txt -osformat genbank -auto -outseq $sample_name.$node_list.revcomp.gb");
                open(FILE, "$sample_name.$node_list.revcomp.gb") or die "$!";
                @genbank = <FILE>;
                close(FILE);
                $i=1;
                @new_genbank = ();
                foreach my $genbank (@genbank){
                    if($genbank =~ /LOCUS/){
                        $genbank =~ s/UNC [^\-]+\-[^\-]+\-\d\d\d\d/VRL $date/;
                    }
                    if($genbank =~ /FEATURES/){
                        $genbank =~ s/FEATURES\s+Location\/Qualifiers/FEATURES             Location\/Qualifiers\n     source          1..$positions\n                     \/organism=\"$location mosquito virus $sample_name.$tmpname\"\n                     \/strain=\"$sample_name\"\n                     \/host=\"$mosq_fullname\"\n                     \/country=\"$location\"/;
                    }
                    if($genbank =~ /mRNA/ or $genbank =~ /CDS/ or $genbank =~ /\/note/){
                        next;
                    }
                    if($genbank =~ /gene/){
                        $genbank =~ s/gene/CDS /;
                        my $tmp = $genbank;
                        my $tmp2 = "                     \/gene=\"ORF$i\"\n";
                        $genbank = $tmp.$tmp2;
                        $i++;
                    }
                    push(@new_genbank,$genbank);
                }
                open(NEWFILE, "> $sample_name.$node_list.mod.revcomp.gb") or die "$!";
                print NEWFILE @new_genbank;
                close(NEWFILE);
                system("mv $sample_name.$node_list.mod.gb ../../GENBANK/$sample_name.$node_list.revcomp.gb");
            }

            chdir("..");#CD=Mapping_contigs
            open(NEWFILE, ">> ..\/Statistics.chimera.contigs.$sample_name.tsv") or die "$!";
            if($count_list==0){
                print NEWFILE "\nSample ID\: $sample_name\nContig ID\tContig length\tmapped reads count\tmapped reads \/ total reads (%)\tRPKM\taveraged mapping coverage\n";
            }
            my $p_abundance = ($abundance / ($abundance+$total))*100;
            my $RPKM = $abundance / 1000000 / $positions / 1000;
            print NEWFILE "$sample_name\_contig\_$tmpname\t$positions\t$abundance\t$p_abundance\t$RPKM\t$avg_coverage\n";
            close(NEWFILE);
            $count_list++;
        }
    }
    return 1;
}