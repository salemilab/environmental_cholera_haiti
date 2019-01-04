# SNP Calling Pipeline on whole genome Vibrio cholerae sequences
# Pipeline created by Taylor K. Paisie
# can be run with indivdual samples or parallel on multiple samples
# This pipeline is meant to be run in parallel with multiple samples
########## Programs needed to run pipeline:
########## fastqc, trimmomatic, bowtie2, samtools, picard, and gatk(gatkv3.8.0 was used in this pipeline)

export _JAVA_OPTIONS="-Xmx10g"

REF=/ufrc/data/reference/bowtie2/v_cholerae_o1_2010el_1786

#trimmomatic and fastqc to analyze quality control of the fastq files before and after trimming
ls -1 | parallel 'fastqc -f fastq {}_1.fastq.gz {}_2.fastq.gz' ::: $(ls *_*.fastq.gz | rev | cut -c 12- | rev | uniq)
parallel 'trimmomatic PE {}_1.fastq.gz {}_2.fastq.gz {}_pair_1.fastq.gz U_{}_1.fastq.gz {}_pair_2.fastq.gz U_{}_2.fastq.gz ILLUMINACLIP:/apps/trimmomatic/0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20' ::: $(ls *_*.fastq.gz | rev | cut -c 12- | rev | uniq)
ls -1 | parallel 'fastqc -f fastq {}_pair_1.fastq.gz {}_pair_2.fastq.gz' ::: $(ls *_pair_*.fastq.gz | rev | cut -c 17- | rev | uniq)

# using bowtie2 for reference based mapping to a reference genome
# after sam file is created the sam file is converted to a bam file
parallel 'bowtie2 -p 8 -x /ufrc/data/reference/bowtie2/v_cholerae_o1_2010el_1786 -1 {}_pair_1.fastq.gz -2 {}_pair_2.fastq.gz -I 0 -X 1200 --un-conc-gz {}_unmapped | samtools view -bS - | samtools sort -O bam -T tmp{}.tmp > {}_sorted.bam' ::: $(ls HC248*_pair_*.fastq.gz | rev | cut -c 17- | rev | uniq)

# Add or replace read groups picard command
# Replace read groups in a BAM file. 
# This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
parallel 'picard AddOrReplaceReadGroups I={}_sorted.bam O={}_re.bam RGID=ID_{} RGLB=LB_{} RGPL=ILLUMINA RGPU=ILLUMINA RGSM=SM_{}' ::: $(ls *_sorted.bam | rev | cut -c 12- | rev | uniq)

# # mark duplicates - Identifies duplicate reads
# This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. 
parallel 'picard MarkDuplicates INPUT={}_re.bam OUTPUT={}_nodups.bam METRICS_FILE=marked_dup_metrics_{}.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true READ_NAME_REGEX="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" OPTICAL_DUPLICATE_PIXEL_DISTANCE=100' ::: $(ls *_re.bam | rev | cut -c 8- | rev | uniq)

# index the output bam file from MarkDuplicates
# Must be indexed for the following GATK commands
parallel 'samtools index {}' ::: *_nodups.bam

# Realigner target creator - Define intervals to target for local realignment
# Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches
parallel 'GenomeAnalysisTK -T RealignerTargetCreator -R /ufrc/salemi/tpaisie/cholera/ref_seq/vibrChol1_N16961.fa -I {}_nodups.bam -o forIndelRealigner_{}.intervals' ::: $(ls *_nodups.bam | rev | cut -c 12- | rev | uniq)

# Indel realigner - Perform local realignment of reads around indels
parallel 'GenomeAnalysisTK -T IndelRealigner -R /ufrc/salemi/tpaisie/cholera/ref_seq/vibrChol1.fa -I {}_nodups.bam -targetIntervals forIndelRealigner_{}.intervals -o {}_realn.bam' ::: $(ls *_nodups.bam | rev | cut -c 12- | rev | uniq)


# Fix mate information - Verify mate-pair information between mates and fix if needed
# This tool ensures that all mate-pair information is in sync between each read and its mate pair
parallel 'picard FixMateInformation I={}_realn.bam O={}_fix.bam SORT_ORDER=coordinate' ::: $(ls *_realn.bam | rev | cut -c 11- | rev | uniq)


# index final bam file for summary statistics
parallel 'samtools index {}' ::: *.bam


# summary statistics on bam files
# Produces a summary of alignment metrics from a SAM or BAM file
parallel 'picard CollectAlignmentSummaryMetrics R=/ufrc/salemi/tpaisie/vvulnificus/ref_seqs/M06-24.fasta I={}.bam O=AlignmentSummaryMetrics_{}.txt ASSUME_SORTED=true MAX_INSERT_SIZE=100000' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# Generate index statistics from a BAM fileThis tool calculates statistics from a BAM index (.bai) file
# The statistics collected include counts of aligned and unaligned reads as well as all records with no start coordinate
parallel 'picard BamIndexStats I={}.bam' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
parallel 'GenomeAnalysisTK -T DepthOfCoverage -R /ufrc/salemi/tpaisie/cholera/ref_seq/vibrChol1.fa -o {} --outputFormat table -I {}.bam' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# Next step is to run freebayes on all the bam files

# reference sequences used to call variants on the bam files
REF=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa


FILE1=$1   # List of bam files to call variants from
FILE2=$2   # output name of the VCF file

freebayes -L ${FILE1} -f $REF -T 0.001 -p 1 -i -X -n 0 -E 3 --min-repeat-size 5 -m 1 -q 20 -R 0 -Y 0 -e 1000 -F 0.5 -C 2 -3 0 -G 1 -! 0 > ${FILE2}
