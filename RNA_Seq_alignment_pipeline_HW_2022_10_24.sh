#Before everyting, make sure you have all the fastq.gz files in one dir. Make sure your annotation file is .gtf file (there are many tutorial online showing how to transform .gff3 to .gtf). Make sure your genome file is a .fasta file.

#make a workdir under the dir of fastq.gz files
mkdir workdir


#copy all the files the the workdir (.fastq.gz and genome file and annotation file to workdir). 
cp *fastq.gz workdir
cd workdir
cp 'your genome file dir (No quote when using the code)' .
cp 'your annotation file dir (No quote when using the code)' .


#Under workdir, make the genome dir
mkdir Genome


#Index the reference genome in STAR (for a 40-core machine). Edit dir and file name accordingly
export PATH=/programs/STAR:$PATH
STAR --runMode genomeGenerate --runThreadN 32 --genomeDir Genome --genomeFastaFiles genome.fasta \
--sjdbGTFfile annotation.gtf --sjdbOverhang 50


#Setting share memory in STAR
STAR --genomeLoad LoadAndExit --genomeDir Genome


#Trimming
#Trimming adaptors for 3' RNA seq. Each trimming uses up to 20 cores.
mkdir samples_for_trimming
ls *.fastq.gz | while read file; do echo /programs/bbmap-38.45/bbduk.sh in=${file} out=${file%%.*}_trimmed.fastq.gz ftl=10 ftr=75; done > barcode_remove.sh
parallel -j 2  < barcode_remove.sh >& log &
##for 85bp 3'RNA seq this trimming process trimes the frist 10bp and the last 10bp. The sequence length is reduced from 85bp to 65 bp.
##Use 'top' to monitor progress


#QC
export PATH=/programs/FastQC-0.11.8:$PATH
ls *_trimmed.fastq.gz | while read file; do echo fastqc ${file} ; done > fastqc.sh
parallel -j 36  < fastqc.sh >& log &
mkdir fastqc_result
mv *.html fastqc_result/
mv *.zip fastqc_result/
zip -r fastqc_result.zip fastqc_result
##download the fastqc_result.zip through filezilla to check the quality of your sequencing result


#Write a script to run STAR with all samples
ls *_trimmed.fastq.gz | while read file; do echo STAR --quantMode GeneCounts --genomeDir Genome --readFilesCommand zcat --runThreadN 4 --readFilesIn ${file} --outFileNamePrefix ${file%%.*} --outFilterMultimapNmax 2  --outFilterMismatchNmax 10 --outSAMtype BAM SortedByCoordinate --genomeLoad LoadAndKeep --limitBAMsortRAM 4000000000; done > STAR.sh


#STAR script name: STAR.sh
#For a 40-core machine, run 9 jobs simultaneously
parallel -j 9  < STAR.sh >& log &
##all process record are stored in log. use 'less log' to check if there is any problem


#Create a gene count matrix:
cat 'the first ReadsPerGene.out.tab in the file list (No quote when using the code)' | cut -f1 > 100.txt
ls *ReadsPerGene.out.tab | while read file; do cat ${file} | cut -f3 > ${file%%.*}_genecount.txt ; done
paste *.txt | tail -n +5 > Gene_count_matrix.txt


#Remove the genome database from the shared memory
STAR --genomeLoad Remove --genomeDir Genome


#file names as column
ls *_trimmed.fastq.gz | while read file; do echo ${file%%.*}; done > column_names.txt


#You can work on gene expression analysis using Gene_count_matrix.txt and column_names.txt