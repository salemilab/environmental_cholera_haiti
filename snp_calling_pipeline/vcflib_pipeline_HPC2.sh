#!/bin/bash
#SBATCH --job-name=vcv_cholera
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --output=vcf_fa_extraction_%j.log
#SBATCH --ntasks=2
#SBATCH --mem=1gb
#SBATCH --time=1:00:00
date;hostname;pwd

module load vcflib

#Defaults
BASES="3"
QUAL="20"
DP="10"
AODP=".75"
GT="0"
INFILE=$1

echo ${INFILE}

# Example submission line
# qsub -v INFILE="Study1.vcf",BASES="3",QUAL="20",DP="10",AODP=".75",GT="0"

function check_var () {
if [ "x" == "x$VAR" ]; then
    echo "Please provide the $NAME variable on the qsub command line"
    exit 1
fi
}

function verify_data_file () {
if [[ ! -s ${TEST_FILE}  ]]; then
    echo -e "\n${TEST_FILE} is either empty or missing. Exiting.\n"
    exit 1
fi
}

function verify_exit () {
if [[ $exit_code == 0 ]]; then
    echo -e "\nCommand: ${command} exited normally"
else
    echo -e "\nCommand: ${command} exited abnormally, exiting"
    exit $exit_code
fi
}

#Check variables:
VAR=$INFILE; NAME="INFILE"; check_var
VAR=$BASES; NAME="BASES"; check_var
VAR=$QUAL; NAME="QUAL"; check_var
VAR=$DP; NAME="DP"; check_var
VAR=$AODP; NAME="AODP"; check_var
VAR=$GT; NAME="GT"; check_var

echo "Starting from the input file: $INFILE"
TEST_FILE=${INFILE}
verify_data_file

#DIST
#vcfdistance < Study1.vcf > Study1_d.vcf
distfile="$(basename ${INFILE} .vcf)_d.vcf"
command="vcfdistance"
echo "Running command: vcfdistance < $INFILE > $distfile"
bash -c "vcfdistance < ${INFILE} > ${distfile}"
exit_code=$?
verify_exit
TEST_FILE=${distfile}
verify_data_file
echo -e "\n\tvcfdistance stage completed successfully\n"

#BASE
#vcffilter -f "BasesToClosestVariant > 3" Study1_d.vcf > Study1_d3.vcf
basefile="$(basename ${distfile} .vcf)_d3.vcf"
command="vcffilter_base"
args="BasesToClosestVariant > ${BASES}"
echo -e "\nRunning command: vcffilter -f \"BasesToClosestVariant > ${BASES}\" ${distfile} > ${basefile}"
bash -c "vcffilter -f \"$args\" ${distfile} > ${basefile}"
exit_code=$?
verify_exit
TEST_FILE=${basefile}
verify_data_file
echo -e "\n\tvcffilter BasesToClosestVariant stage completed successfully\n"

#QUAL
#vcffilter -f "QUAL > 20" Study1_d3.vcf > Study1_d3_q20.vcf
qualfile="$(basename ${basefile} .vcf)_q20.vcf"
command="vcffilter_qual"
args="QUAL > ${QUAL}"
echo -e "\nRunning command: vcffilter -f \"QUAL > ${QUAL}\" ${basefile} > ${qualfile}"
bash -c "vcffilter -f \"$args\" ${basefile} > ${qualfile}"
exit_code=$?
verify_exit
TEST_FILE=${qualfile}
verify_data_file
echo -e "\n\tvcffilter QUAL stage completed successfully\n"

#DP
#vcffilter -g " DP > 10 " Study1_d3_q20.vcf > Study1_d3_q20_DP5.vcf
dpfile="$(basename ${qualfile} .vcf)_DP10.vcf"
command="vcffilter_dp"
args=" DP > ${DP} "
echo "Running command: vcffilter -g \"DP > ${DP}\" ${qualfile} > ${dpfile}"
bash -c "vcffilter -g \"$args\" ${qualfile} > ${dpfile}"
exit_code=$?
verify_exit
TEST_FILE=${dpfile}
verify_data_file
echo -e "\n\tvcffilter DP stage completed successfully\n"


#AODP
#vcffilter -g "( ( AO / DP ) > .75 )" Study1_d3_q20_DP5.vcf > Study1_d3_q20_DP5_AOF.vcf
aodpfile="$(basename ${dpfile} .vcf)_AOF.vcf"
command="vcffilter_aodp"
args="( ( AO / DP ) > ${AODP} )"
echo "Running command: vcffilter -g \"( ( AO / DP ) > ${AODP} )\" ${dpfile} > ${aodpfile}"
bash -c "vcffilter -g \"$args\" ${dpfile} > ${aodpfile}"
exit_code=$?
verify_exit
TEST_FILE=${aodpfile}
verify_data_file
echo -e "\n\tvcffilter AO/DP stage completed successfully\n"

#GT
#vcffilter -g " GT = 0 " Study1_d3_q20_DP5.vcf > Study1_d3_q20_DP5_GT.vcf
gtfile="$(basename ${dpfile} .vcf)_GT.vcf"
command="vcffilter_gt"
args=" GT = ${GT} "
echo "Running command: vcffilter -g \" GT = ${GT} \" ${dpfile} > ${gtfile}"
bash -c "vcffilter -g \"$args\" ${dpfile} > ${gtfile}"
exit_code=$?
verify_exit
TEST_FILE=${gtfile}
verify_data_file
echo -e "\n\tvcffilter GT stage completed successfully\n"


#COMBINE
#vcfcombine Study1_d3_q20_DP5_AOF.vcf Study1_d3_q20_DP5_GT.vcf > Study1_Filtered.vcf
combinefile="$(basename ${INFILE} .vcf)_Filtered.vcf"
command="vcfcombine"
echo "Running command: vcfcombine ${aodpfile} ${gtfile} > ${combinefile}"
bash -c "vcfcombine ${aodpfile} ${gtfile} > ${combinefile}"
exit_code=$?
verify_exit
TEST_FILE=${combinefile}
verify_data_file
echo -e "\n\tvcfcombine stage completed successfully\n"


#EXTRACT
if [[ -f vcf_fa_extractor.py ]]; then
    exec="./vcf_fa_extractor.py"
elif [[ -d vcf_fa_extractor ]]; then
    if [[ -f vcf_fa_extractor/vcf_fa_extractor.py ]]; then
        exec="./vcf_fa_extractor/vcf_fa_extractor.py"
    else
        rm -rf vcf_fa_extractor
        module load git
        git clone https://github.com/moskalenko/vcf_fa_extractor
        echo -e "\n\tDownloading vcf_fa_extractor\n"
    fi
else
    module load git
    git clone https://github.com/moskalenko/vcf_fa_extractor
    echo -e "\n\tDownloading vcf_fa_extractor\n"
fi
fafile="$(basename ${combinefile} .vcf).fa"
command="vcf_fa_extractor.py"
${exec} -i ${combinefile}
exit_code=$?
verify_exit
TEST_FILE=${fafile}
verify_data_file
echo -e "\n\tvcf_fa_extractor stage completed successfully\n"
echo -e "\n\tAll stages completed successfully. Exiting.\n"
