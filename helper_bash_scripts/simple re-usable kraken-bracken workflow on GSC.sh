#Getting onto the GSC
#Submit jobs of 
#Choosing the remote host is numbers
#N105 - centos6
# or N106 - centos7 - set everything up here
#From numbers transition into high performance computer
ssh gphost08
#Resume the screen - work from within - if it says there are no screens, then create one
screen -S kraken
###########This is very important
mamba activate kraken2 
screen -r kraken 
#Work from within
#/projects/molonc/scratch/krekhno
#First remove all the encrypted files since they just take up space
rm egads/EGAD*/*.c4gh

#This is the script to sort, extract, kraken and bracken the reads
#No need to execute any of these commands, skip straight to executing the full script - the line that starts with bash
nano helper_bash_scripts/reusable_sort_class_counts.sh
#!/bin/bash
set -e 
set -u
set -o pipefail
for i in egads/EGAD*/*.bam*
do
    echo "file is $i"
    #Extract the file name and dataset name for easier naming
    file_name="$(basename $i | cut -d. -f1)"
    dataset_name="$(basename "$(dirname $i)")"
    dataset_file="${dataset_name}_${file_name}"
    #Extract these reads to a fastq file for input to kraken:
    echo 'extracting and sorting unmapped reads into a new bam file'
    samtools view -b -f 4 $i > unmapped/${file_name}_unmapped_sort.bam
    echo 'sorting unmapped reads'
    samtools sort -n -o unmapped/${file_name}_unmapped_true_sort.bam unmapped/${file_name}_unmapped_sort.bam
    echo 'extracting sorted unmapped reads into fastq files'
    bedtools bamtofastq -i unmapped/${file_name}_unmapped_true_sort.bam -fq unmapped/${file_name}_unmapped_1.fq -fq2 unmapped/${file_name}_unmapped_2.fq
    echo 'Classifying the reads with kraken2'
    kraken2 --db ./kraken2_library --report-zero-counts \
        --output kraken2_outputs/${dataset_file}_kraken_out.txt --report kraken2_outputs/${dataset_file}_kraken_report.txt \
        --report-minimizer-data --minimum-hit-groups 3 --threads 8 --paired unmapped/${file_name}_unmapped_1.fq unmapped/${file_name}_unmapped_2.fq 
    bracken -d ./kraken2_library -i kraken2_outputs/${dataset_file}_kraken_report.txt -r 100 \
        -l G -t 10 -o bracken_outputs/${dataset_file}_sample.bracken -w bracken_reports/${dataset_file}_sample.breport
done
#To execute the script
bash helper_bash_scripts/reusable_sort_class_counts.sh
#Check results by looking into the bracken_reports folder - if there are more folders there, then you are good.

