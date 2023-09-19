#Getting onto the GSC
#Submit jobs of 
#Choosing the remote host is numbers
#N105 - centos6
# or N106 - centos7 - set everything up here
#Install mamba
#Path to conda executable - /gsc/software/linux-x86_64-centos7/Anaconda3-4.4.0/bin
#Add samtools to path - if samtools is not working in 
#Added this command to ~/.bashrc to be executed every time I start - also add bedtools there
export PATH="$PATH:/gsc/software/linux-x86_64-centos7/samtools-1.9/bin/"
export PATH="$PATH:/gsc/software/linux-x86_64-centos7/bedtools-2.27.1/bin/"
#Work from within
#/projects/molonc/scratch/krekhno
#From home directory make a link
ln -s /projects/molonc/scratch/krekhno
#Use mamba-forge
#Download install script 
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
#Install to home directory - /projects/molonc/scratch/krekhno/mambaforge
#Install kraken2 using mamba
mamba create -n kraken2 -c conda-forge -c bioconda -c defaults kraken2 bracken
#From here work with mamba
####################################
#Have used Aspera to download files, and decrypted them - check using Aspera to download EGA files.sh in helper_bash_scripts
#Now check read lenght to build the database
samtools view egads/EGAD00001009184/EGAR00003355683_EGAS00001006343_SHAH_H000566_T01_01_WG01_A48967_4_lanes_dupsFlagged.bam.1660354221626 | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
#Confirm that the reads are paired
samtools view -f 0x1 egads/EGAD00001009184/EGAR00003355683_EGAS00001006343_SHAH_H000566_T01_01_WG01_A48967_4_lanes_dupsFlagged.bam.1660354221626 | head -n 1

mamba activate kraken2 
#Roughly follow this link
#https://hackmd.io/@AstrobioMike/kraken2-bracken-standard-build
#Download the kraken/bracken database and expand it here
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220607.tar.gz -O kraken2-std.tar.gz 
tar -xzf kraken2-std.tar.gz

#Figure out read length - most common is 125bp
#Next we need to extract unmapped reads from the bam file, extract fastq files from the unmapped bam file, and then do kraken2 alignment, and Bracken abundance estimation
#ssh into gphost08, run these commands on a screen
screen -S kraken
mamba activate kraken2 
#Extract unmapped reads
samtools view -b -f 4 egads/EGAD00001009184/EGAR00003355683_EGAS00001006343_SHAH_H000566_T01_01_WG01_A48967_4_lanes_dupsFlagged.bam.1660354221626 > test_unmapped.bam
#Sort the reads
samtools sort -n -o test_unmapped_sort.bam test_unmapped.bam
#Extract these reads to a fastq file for input to kraken:
#Extract the unmapped sorted reads into fastq files
bedtools bamtofastq -i test_unmapped_sort.bam -fq test_unmapped_sort_1.fq -fq2 test_unmapped_sort_2.fq
#Use Kraken2 to classify these reads
kraken2 --db ./kraken2_library --report-zero-counts \
 --output kraken2_outputs/test_kraken_out.txt --report kraken2_outputs/test_kraken_report.txt \
 --report-minimizer-data --minimum-hit-groups 3 --threads 8 --paired test_unmapped_sort_1.fq test_unmapped_sort_2.fq 
 #Use Bracken to estimate bacterial abundance at the genus level setting read length at 125bp, and setting read threshold of 10 reads for re-restimation
 #Don't have the 125bp kmers in default, so use 100 instead
bracken -d ./kraken2_library -i kraken2_outputs/test_kraken_report.txt -r 100 \
 -l G -t 10 -o bracken_outputs/test_sample.bracken -w bracken_reports/test_sample.breport
##################################
#All this worked on gphost08, so just run a script  for the files we have
#All seem to have worked so let's run a script
nano helper_bash_scripts/kraken_class_count.sh
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
    echo 'extracting and sorting unmapped reads into a new bam file'
    samtools view -b -f 4 $i | samtools sort -o unmapped/${file_name}_unmapped_sort.bam
    #For whatever reason this does not actually sort the reads, so have to re-do sorting and getting fastq files and the rest.
done
#Exit nano
#Write a new script to sort, extract, kraken and bracken the reads
nano helper_bash_scripts/kraken_sort_class_count.sh
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
#Something is really odd with extracting the fastq files - very few reads extracted - try manually extracting the reads - this should be the same file as the test files
bedtools bamtofastq -i unmapped/EGAR00003355683_EGAS00001006343_SHAH_H000566_T01_01_WG01_A48967_4_lanes_dupsFlagged_unmapped_sort.bam -fq test_683_1.fq -fq2 test_683_2.fq
#Same result as the script - try sorting the bam file again
samtools sort -n -o true_683_sort.bam unmapped/EGAR00003355683_EGAS00001006343_SHAH_H000566_T01_01_WG01_A48967_4_lanes_dupsFlagged_unmapped_sort.bam
bedtools bamtofastq -i true_683_sort.bam -fq test_683_1.fq -fq2 test_683


#Try copying files from the server using xfer and rsync
#In mobaxterm - local terminal move to a a reasonable directory and try copying files
#Either option depending on the computer
rsync -a zkrekhno@xfer.bcgsc.ca:/home/zkrekhno/krekhno/bracken_outputs ./
rsync -a zkrekhno@xfer.bcgsc.ca:/home/zkrekhno/krekhno/bracken_reports ./

#Download a pyega3 log
rsync -a zkrekhno@xfer.bcgsc.ca:/home/zkrekhno/krekhno/pyega3_output.log ./






