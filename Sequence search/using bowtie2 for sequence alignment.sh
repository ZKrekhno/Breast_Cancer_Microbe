#Getting onto the GSC
#Submit jobs of 
#Choosing the remote host is numbers
#N105 - centos6
# or N106 - centos7 - set everything up here
#From numbers transition into high performance computer
ssh gphost08
#This is where the bowtie2 software is located 
#/gsc/software/linux-x86_64-centos7/bowtie2-2.3.4.1/bin/
screen -S test_bowtie
#Add it to PATH for easy calling
export PATH="/gsc/software/linux-x86_64-centos7/bowtie2-2.3.4.1/bin":$PATH;
#Make directory pathogenic_genes - inside the directory create individual fasta files for each sequence
bowtie2-build Sample_1_NoTaxonomy_Rast.CDS.1956.fasta,Sample_1_NoTaxonomy_Rast.CDS.746.fasta test_index
#Try making an index from a conjoined file
cat Sample_1_NoTaxonomy_Rast.CDS.1956.fasta Sample_1_NoTaxonomy_Rast.CDS.746.fasta > test_big_file.fasta
bowtie2-build test_big_file.fasta test_big_index
#Test the big index with synthetic reads (really basic just some hits and misses done manually)
bowtie2 -x test_big_index -U alignment_test.fasta -f --local -S test_big_bowtie.sam
#Create a sub-directory alignment_counting_tests and move everything there to keep it clean
mv * alignment_counting_tests/
#Try immediately converting the alignment output into a bam file and. Use no-unal flag to only keep aligned reads
bowtie2 -x Curated -U alignment_counting_tests/alignment_test.fasta -f --local --no-unal | samtools view -bS - > test_big_bowtie.bam
#We can try counting this with featurecounts


#Install featurecounts (subread)
mamba create -c bioconda -n featurecounts subread 
mamba activate featurecounts
#Create a simple annotation file and use it for this analysis
#Run a test count with 5 threads
featureCounts -a test_annot.saf -F SAF -T 5 -g GeneID --verbose -o test_big_counts.txt test_big_bowtie.bam 
#All this works fine - so can move to preparing the index and saf file for all the genes of interest.
#The file I got from Zack needed some massaging which was done in R to make a big fasta file and a proper saf file
#Make the index using the curated_CDS - this is done in the pathogenic genes folder
bowtie2-build Curated_cds.fasta Curated
#Test the featureCounts with complete annotation
featureCounts -a Annot_cds.saf -F SAF -T 5 -g GeneID --verbose -o Sample_1_big_counts.txt alignment_counting_tests/test_big_bowtie.bam 
#All this works, so we can transition to using actual files for mapping and counting the reads
#Set up a smart script to run the alignment and feature counting that uses names of files and checks where these analyses have already been done
nano helper_bash_scripts/map_and_count.sh
#!/bin/bash
set -e 
set -u
set -o pipefail
for i in unmapped/*unmapped_true_sort.bam
do
    echo "file is $i"
    #Extract the file name and dataset name for easier naming
    file_name="$(basename $i | cut -d. -f1)"
    file_name="${file_name//_true_sort/}"
    echo "short name is $file_name"
        # Check if the file file_name_mapped_cds.bam already exists - and proceed to map against cds of itnerest if it doesn't
    if [ -e "pathogenic_genes/${file_name}_mapped_cds.bam" ]; then
        echo "File ${file_name}_mapped_cds.bam already exists in the directory."
    else #We can align the reads to the cds of interes
        echo "Aligning reads to cds for file: $file_name"
        bowtie2 -x pathogenic_genes/Curated  -1 unmapped/${file_name}_1.fq \
         -2 unmapped/${file_name}_2.fq -q --local --no-unal | samtools view -bS - > pathogenic_genes/${file_name}_mapped_cds.bam
    fi
    #Now we can count the number of reads
    if [ -e "pathogenic_genes/count_results/${file_name}_counts.txt" ]; then
        echo "File ${file_name}_counts.txt already exists in the directory."
    else #We can count
        echo "Counting aligned reads to cds for file: $file_name"
        #Use -p flag to indicate that reads are paired
        featureCounts -p -a pathogenic_genes/Annot_cds.saf -F SAF -T 5 -g GeneID --verbose -o pathogenic_genes/count_results/${file_name}_counts.txt pathogenic_genes/${file_name}_mapped_cds.bam
    fi
done
#That is done - only 21 files had 1 or more reads map to our genes of interest
#Do this in pathogenic_genes/count_results
mamba deactivate
paste *counts.txt > comb_succ_counts.txt


rsync -a zkrekhno@xfer.bcgsc.ca:/home/zkrekhno/krekhno/pathogenic_genes/count_results/comb_succ_counts.txt ./



