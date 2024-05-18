Detailed methods and intermediate files for study "Chromosome-level genome assembly of an auxotrophic strain of the pathogenic yeast *Candida parapsilosis*"

## Files

Folder `assembly`

* `SR23sample-flye.fa.gz` original Flye assembly
* `SR23-rDNAp.fa` polished copy of rDNA repeat
* `SR23_mtgenome.fa` mitochondrial genome (from RefSeq)
* `manual1-regions.list` manually created file listing which regions to connect to final chromosomes
* `manual1-fasta.list` list of filenames used in `manual1-regions.list`, both these files are inputs for script `assembly.pl`
* `canParSR23v1.fa.gz` final assembly after polishing

Folder `annotation`

* `transcripts.psl.gz` alignments of assembled transcripts to the genome
* `au-other-tr-train.tar.gz` Augustus parameters for *C. parapsilosis* used for prediction
* `au-1-tr.gp.gz` Augustus predictions in genepred format
* `canParCDC317-mapped.gp.gz` *C. parapsilosis* strain CDC317 genes mapped to SR23 assembly
* `canParB1.fa.gz` assembly of *C. parapsilosis* strain CLIB214
* `canParB1.gtf.gz` annotation of *C. parapsilosis* strain CLIB214
* `canParB1-mapped.gp.gz` *C. parapsilosis* strain CLIB214 genes mapped to SR23 assembly
* `to-remove-manual.gp` Augustus predictions removed from the final gene annotation
* `to-add-manual.gp` Genes based on mappings from other strains added to the final gene annotation
* `genes1-all.gp.gz` Final gene annotation

## Custom scripts

In addition to standard bioinformatics software tools, we used several custom scripts available at https://github.com/fmfi-compbio/assembly-scripts

* `extract_fasta.pl`
* `assembly.pl`
* `filter-gtf`, uses config file `about/add_cds.about`

## Assembly
```bash
# basecalling nanopore data (HAC mode)
guppy_basecaller --num_callers 4 --cpu_threads_per_caller 4 --device auto -i reads -s guppy -r --flowcell FLO-MIN106 --kit SQK-LSK109

# checking nanopore read quality with NanoPlot
NanoPlot --fastq SR23-N.fastq.gz -o nanoplot-SR23-N
NanoPlot --summary guppy/sequencing_summary.txt -o nanoplot-summary

# Illumina read quality by fastqc
fastqc SR23-I_1.fastq.gz
fastqc SR23-I_2.fastq.gz

# obtaining 40% of nanopore reads longer than 5kbp
# by random smapling using a Perl one0liner
zcat SR23-N.fastq.gz | perl -ne 'BEGIN { srand(42); } $s.=$_; if($.%4==0) { if(length($_)>5*1000 && rand(1)<0.4) { print $s; } $s=""; }' | gzip -c > SR23sample-N.fastq.gz

# running flye assembler on nanopore data
flye -t 4 --nano-raw SR23sample-N.fastq.gz --out-dir SR23sample-flye.fa.tmp 2>> SR23sample-flye.fa.log
mv SR23sample-flye.fa.tmp/assembly.fasta SR23sample-flye.fa

# aligning nanopore reads to assembly
# in both paf and bam formats
minimap2 -c -x map-ont --secondary=no -t 4 SR23-flye.fa SR23-N.fastq.gz > SR23-flye-MAP-SR23-N.paf
minimap2 -a -x map-ont --secondary=no -t 4 SR23-flye.fa SR23-N.fastq.gz | samtools view -S -b - | samtools sort - -o SR23-flye-MAP-SR23-N.bam
samtools index SR23-flye-MAP-SR23-N.bam

# extract one full rDNA repeat discovered based on
# similarity with other known strains
# using a custom script
extract_fasta.pl SR23sample-flye.fa contig_6 13465 21051 + rDNA-SR23 > SR23-rDNA.fa

# aligning illumina reads to rDNA repeat
bwa index SR23-rDNA.fa -p SR23-rDNA-I.bam.tmp
bwa mem -t 4 SR23-rDNA-I.bam.tmp illumina_1.fastq.gz illumina_2.fastq.gz | samtools view -b - | samtools sort - -o SR23-rDNA-I.bam
samtools index SR23-rDNA-I.bam

# polishing rDNA repeat by pilon
java -Xmx20G -jar /opt/broad/pilon-1.21.jar --genome SR23-rDNA.fa --outdir SR23-rDNA-pilon --changes --tracks --frags SR23-rDNA-I.bam > SR23-rDNAp.fa.log
perl -lne 's/^(>.*)_pilon\s*$/$1/; print' SR23-rDNA-pilon/pilon.fasta > SR23-rDNAp.fa

# putting together the assembly from pieces using a custom script
manual_assembly.pl -j manual1-joins.list manual1-fasta.list manual1-regions.list > manual1.fa

# polishing assembly by pilon, 3 iterations
# starting with manual.fa, result manual1ppp.fa

# iteration 1, align reads
bwa index manual1.fa -p manual1-I.bam.tmp
bwa mem -t 4 manual1-I.bam.tmp illumina_1.fastq.gz illumina_2.fastq.gz | samtools view -b - | samtools sort - -o manual1-I.bam
samtools index manual1-I.bam
# iteration 1, pilon
java -Xmx20G -jar /opt/broad/pilon-1.21.jar --genome manual1.fa --outdir manual1-pilon --changes --tracks --frags manual1-I.bam > manual1p.fa.log
perl -lne 's/^(>.*)_pilon\s*$/$1/; print'  manual1-pilon/pilon.fasta > manual1p.fa
# iteration 2, align reads
bwa index manual1p.fa -p manual1p-I.bam.tmp
bwa mem -t 4 manual1p-I.bam.tmp illumina_1.fastq.gz illumina_2.fastq.gz | samtools view -b - | samtools sort - -o manual1p-I.bam
samtools index manual1p-I.bam
# iteration 2, pilon
java -Xmx20G -jar /opt/broad/pilon-1.21.jar --genome manual1p.fa --outdir manual1p-pilon --changes --tracks --frags manual1p-I.bam > manual1pp.fa.log
perl -lne 's/^(>.*)_pilon\s*$/$1/; print' manual1p-pilon/pilon.fasta > manual1pp.fa
# iteration 3, align reads
bwa index manual1pp.fa -p manual1pp-I.bam.tmp
bwa mem -t 4 manual1pp-I.bam.tmp illumina_1.fastq.gz illumina_2.fastq.gz | samtools view -b - | samtools sort - -o manual1pp-I.bam
samtools index manual1pp-I.bam
# iteration 3, pilon
java -Xmx20G -jar /opt/broad/pilon-1.21.jar --genome manual1pp.fa --outdir manual1pp-pilon --changes --tracks --frags manual1pp-I.bam > manual1ppp.fa.log
perl -lne 's/^(>.*)_pilon\s*$/$1/; print' manual1pp-pilon/pilon.fasta > manual1ppp.fa
    
# file manual1ppp.fa was renamed to canParSR23v1.fa
# and used as the final assembly (submitted to ENA)
```


## Protein-coding gene annotation

```bash
# trimming RNA-seq reads (all three conditions combined) 
trimmomatic PE SR23_combined_1.fastq SR23_combined_2.fastq SR23_combined_trimmed_1.fastq SR23_combined_trimmed_u.fastq SR23_combined_trimmed_2.fastq SR23_combined_trimmed_u2.fastq ILLUMINACLIP:/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 2> SR23_combined_trimmed_1.fastq.log
gzip SR23_combined_trimmed_1.fastq SR23_combined_trimmed_u.fastq SR23_combined_trimmed_2.fastq SR23_combined_trimmed_u2.fastq

# assembling transcripts separately on paired and single reads
# Single reads come from short transcripts
# then both sets are given unique suffix and concatenated
Trinity --max_memory 40G --CPU 8 --seqType fq --left SR23_combined_trimmed_1.fastq.gz --right SR23_combined_trimmed_2.fastq.gz --full_cleanup --output SR23_combined_tmp_trinityp
mv SR23_combined_tmp_trinityp.Trinity.fasta SR23_combined_tr_paired.fa
Trinity --max_memory 40G --CPU 8 --seqType fq --single SR23_combined_trimmed_u.fastq.gz --full_cleanup --output SR23_combined_tmp_trinitys
mv SR23_combined_tmp_trinitys.Trinity.fasta SR23_combined_tr_single.fa
perl -ne 's/>TRINITY/>TRINITYP/; print' SR23_combined_tr_paired.fa > SR23_combined_tr.fa
perl -ne 's/>TRINITY/>TRINITYS/; print' SR23_combined_tr_single.fa >> SR23_combined_tr.fa

# mapping reads to genome, keeping only strong matches
blat -maxIntron=5000 genome.fa SR23_combined_tr.fa SR23_combined_tr.tmp_psl
pslCDnaFilter -minId=0.95 -minCover=0.75 -globalNearBest=0 SR23_combined_tr.tmp_psl transcripts.psl 2> SR23_combined_tr.psl.log

# creating hints for Augustus from aligned transcripts
sort -k14,14 -k16,16g transcripts.psl > transcripts.hints.gff.tmp.psl
/usr/local/share/augustus-3.2.3/augustus-3.2.3/scripts/blat2hints.pl --in=transcripts.hints.gff.tmp.psl --out=transcripts.hints.gff

# running augustus
augustus --alternatives-from-evidence=false --alternatives-from-sampling=false --uniqueGeneId=true --AUGUSTUS_CONFIG_PATH=/path/to/au-other-tr-train --species=canParB1 --hintsfile=transcripts.hints.gff --extrinsicCfgFile=/usr/local/share/augustus-3.2.3/augustus-3.2.3/config/extrinsic/extrinsic.ME.cfg genome.fa > au-1-tr.orig.gtf

# remove unneeded parts of the gtf file and convert to USC genePred format
perl -lane 'print unless /^#/ || $F[1] ne "AUGUSTUS" || $F[2] eq "gene" || $F[2] eq "transcript"' au-1-tr.orig.gtf > au-1-tr.gp.tmp.gtf
gtfToGenePred -genePredExt au-1-tr.gp.tmp.gtf au-1-tr.gp

# obtain genome and annotation of C.parapsilosis strain CDC317
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/765/GCA_000182765.2_ASM18276v2/GCA_000182765.2_ASM18276v2_genomic.fna.gz -O canParCDC317.fa.gz
gunzip canParCDC317.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/765/GCA_000182765.2_ASM18276v2/GCA_000182765.2_ASM18276v2_genomic.gff.gz -O canParCDC317.gff.gz
gunzip canParCDC317.gff.gz

# omit pseudogenes and map genes to our assembly using liftoff tool v1.6.3
grep -v -E '(gene_biotype=pseudogene|pseudo=true)' canParCDC317.gff > tmp.gff
gff3ToGenePred tmp.gff tmp.gp
perl -lane 'print if $F[5]<$F[6]' tmp.gp > tmp2.gp
genePredToGtf file tmp2.gp canParCDC317-filtered.gtf
# require 80% coverage and identity, look for additional gene copies
liftoff genome.fa canParCDC317.fa -g canParCDC317-filtered.gtf -o canParCDC317-mapped.gtf -u canParCDC317-unmapped.txt -dir canParCDC317-mapped-tmp -a=0.8 -s=0.8 -flank 0.1 -copies -sc 0.9 -infer_genes
# add CDS records for the longest ORF in each mapped gene
# and convert to genepred
perl -lane 's/(transcript_id ")rna-/$1/; print if $F[2] eq "exon"' canParCDC317-mapped.gtf > tmp-exons.gtf
filter-gtf -p "INPUT tmp-exons.gtf OUTPUT tmp-CDS.gtf" about/add_cds.about genome.fa
gtfToGenePred -genePredExt tmp-CDS.gtf canParCDC317-mapped.gp


# similar process for strain CLIB214, here denoted as canParB1
liftoff genome.fa canParB1.fa -g canParB1.gtf -o canParB1-mapped.gtf -u canParB1-unmapped.txt -dir canParB1-mapped-tmp -a=0.8 -s=0.8 -flank 0.1 -copies -sc 0.9 -infer_genes
# add CDS records for the longest ORF in each mapped gene
# and convert to genepred
perl -lane 'print if $F[2] eq "exon"' canParB1-mapped.gtf > tmp-exons.gtf
filter-gtf -p "INPUT tmp-exons.gtf OUTPUT tmp-CDS.gtf" about/add_cds.about genome.fa
gtfToGenePred -genePredExt tmp-CDS.gtf canParB1-mapped.gp
```

### Manual adjustments

Selected genes were visually inspected in a local copy of UCSC genome browser, using files `transcripts.psl`, `au-1-tr.gp`, `canParCDC317-mapped.gp`, `canParB1-mapped.gp`. Where appropriate, Augustus prediction was removed or replaced with a protein mapped from another *C. parapsilosis* strain (typically the reference genome CDC317 but in three cases from strain CLIB214). Some CDC317 genes not present in Augustus predictions were added (typically short genes). See the relevant files above.
