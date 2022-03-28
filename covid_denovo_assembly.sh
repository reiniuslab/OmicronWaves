#!/bin/bash

trim_data(){
	# $1 in1
	# $2 in2
	# $3 folder clean
	# $4 basename

	mkdir ${3}/fastp_reports/${4}
	
	fastp --in1 $1 --in2 $2 --out1 ${3}/${1##*/} --out2 ${3}/${2##*/} --json ${3}/fastp_reports/${4}/fastp.json --html ${3}/fastp_reports/${4}/fastp.html -R ${4}
}
export -f trim_data

fix_fasta_header(){
	awk '/>/{sub(">","&"FILENAME" ");sub(/denovo\//,x)}1' $1
}
export -f fix_fasta_header

align_covid(){
	# $1 index
	# $2 Read1
	# $3 Read2
	# $4 output
	minimap2 -ax sr $1 $2 $3 | samtools sort -O bam -o ${4}
	samtools index ${4}
}
export -f align_covid

##

## extract fastq

mkdir fastq

ulimit -n 4000 # increase open file limit
bcl2fastq --no-lane-splitting -R 220204_NB552561_0016_AHG7YJBGXJ --sample-sheet SampleSheet_NB552561_0016_AHG7YJBGXJ.csv -o fastq -r 64 -p 64 -w 64
bcl2fastq --no-lane-splitting -R 220204_NB552564_0016_AHLMC7BGXJ --sample-sheet SampleSheet_NB552564_0016_AHLMC7BGXJ.csv -o fastq -r 64 -p 64 -w 64
bcl2fastq --no-lane-splitting -R 220305_NB552561_0017_AHLMCFBGXJ --sample-sheet SampleSheet_NB552561_0017_AHLMCFBGXJ.csv -o fastq -r 64 -p 64 -w 64

## quality trim data

mkdir clean clean/fastp_reports

parallel -j 32 --link --rpl '{%} s:.*/::; s:_S\d{1,3}.*::;' \
trim_data {1} {2} clean {1%} \
::: fastq/COVIDSeq/*_R1_001.fastq.gz ::: fastq/COVIDSeq/*_R2_001.fastq.gz

multiqc -d -n qc_report clean/fastp_reports

## de novo assembly

mkdir denovo

parallel -j 12 --link --rpl '{%} s:.*/::; s:_S\d{1,3}.*::;' \
./MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
-t 16 \
--min-contig-len 5000 \
-1 {1} \
-2 {2} \
-o denovo/{1%} \
--out-prefix {1%} \
::: clean/*_R1_001.fastq.gz ::: clean/*_R2_001.fastq.gz > denovo.log

# scaffold contigs together
parallel -j 32 ./RagTag-2.1.0/ragtag.py scaffold -o {//} Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz {} ::: denovo/*/*.contigs.fa

# rename fasta headers and merge
parallel --plus fix_fasta_header {} ">" {..}.header.fa ::: denovo/*/ragtag.scaffold.fasta
cat denovo/*/*.header.fa > denovo/merged.scaffold.fa

# Assembly N50 stats
grep -o "N50.*bp" denovo/*/*.log > denovo/n50.txt

## Align raw data to SARS-CoV-2 reference genome and generate alignment statistics

mkdir bam

parallel -j 64 --link --rpl '{%} s:.*/::; s:_S\d{1,3}.*::;' \
align_covid Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz {1} {2} bam/{1%}.bam \
::: clean/*_R1_001.fastq.gz ::: clean/*_R2_001.fastq.gz

mkdir alignment_stats

parallel -j 64 samtools flagstat {} ">" alignment_stats/flagstat_{/.}.txt ::: bam/*.bam
parallel -j 64 samtools coverage {} ">" alignment_stats/coverage_{/.}.txt ::: bam/*.bam
