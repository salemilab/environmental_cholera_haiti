#**********************************************#

#Author: Taylor K Paisie
#Must have bam files made
#Uses bam files that were made with the bam_pipeline.sh script
#This script will give you a SNP alignment with "?" as gaps

########## Programs needed to run pipeline:
########## samtools and gatk


#**********************************************#

module load intel/2016.0.109 
module load openmpi/1.10.2
module load gcc/5.2.0
module load parallel
module load gatk

# reference sequences used on the samples
REF=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa

export _JAVA_OPTIONS="-Xms1g -Xmx10g"

# VCF file that was produced by freebayes from the bam_pipeline.sh script
VCFFILE=$1

# step 2 - extract SNPs and Indels from each vcf file
# extracts SNPS
gatk SelectVariants -R $REF -V ${VCFFILE} --select-type-to-include SNP -O snps_${VCFFILE}
gatk SelectVariants -R $REF -V ${VCFFILE} --select-type-to-include INDEL -O indels_${VCFFILE}

# step 3 - filter the SNP and Indel files
# filters the SNP vcf
gatk VariantFiltration -R $REF -V snps_${VCFFILE} --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O filtered_snps_${VCFFILE}
gatk VariantFiltration -R $REF -V indels_${VCFFILE} --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O filtered_indels_${VCFFILE}


# step 4 - BQSR #1
parallel 'gatk BaseRecalibrator -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}.bam --known-sites filtered_snps_javiana_uncal_041418.vcf --known-sites filtered_indels_javiana_uncal_041418.vcf -O {}.table' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# step 6 - Analyze the BQSR reports from the base recalibration steps (steps 4 & 5)
parallel 'gatk ApplyBQSR -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}.bam --bqsr-recal-file {}.table -O recal_{}.bam' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)


# step 8 - rename bam files (if you want) and call variants with freebayes
# renaming recalibrated bam files
for i in $(ls recal_*.bam | rev | cut -c 5- | rev | uniq)
	do
		mv ${i}.bam ${i}.bam
		mv ${i}.bai ${i}.bai
	done


LIST=$2 # List of bam files to call variants from
VCF=$3 # output name of the VCF file

#calling variants on recalibrated bam files (will have only one vcf file as the output)
freebayes -L ${LIST} -v ${VCF} -f $REF 0.001 -p 1 -i -X -n 0 -E 3 --min-repeat-size 5 -m 1 -q 20 -R 0 -Y 0 -e 1000 -F 0.5 -C 2 -3 0 -G 1 -! 0

# step 9 - filter vcf file from freebayes for SNPs only
vcffilter -f "TYPE = snp" ${VCF} > snps_${VCF}

# step 10 - compressing and indexing the snp only vcf file & variant normalization of the snp only vcf file & decompressing the vcf file

#compress vcf
bgzip snps_${VCF}

#index vcf
tabix -p vcf snps_${VCF}.gz

# normalizing the variant vcf file
bcftools norm -f $REF -o norm_snps_${VCF}.gz snps_${VCF}.gz

#decompresses the output variant normalized file
bgzip -d norm_snps_${VCF}.gz


# step 11 - filter normalized vcf file and create SNP alignment fasta file
# Need the have the script vcf_fa_extractor.py in order to create the SNP alignment
# vcflib_pipeline_HPC2.sh uses the vcf_fa_extractor.py automatically to make the SNP alignment output
bash vcflib_pipeline_HPC2.sh norm_snps_${VCF}
