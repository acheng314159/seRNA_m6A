##1.download SRA
bash loop_prefetch.sh SRPxxxx 10 100 >loop_prefetch.log   ### SRPxxxx.txt in this directory

##2.SRA to fastq
ls *sra |while read id ;do /Sum/lirui/soft/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --split-3 $id;done

##3.alignment
bowtie2 -p  5 -x /Sum/lirui/database/bowtie2_index/hg19bowtie2index -U /Sum/lirui/project/H3K27ac/fastq/SRRxxx.fastq.gz | samtools sort  -O bam  -@ 5 -o - >/Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.bam

samtools markdup /Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.bam /Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.rmdup.bam

samtools index /Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.rmdup.bam > /Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.rmdup.bam.bai

##4.callpeak
macs2 callpeak -t /Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.rmdup.bam -c /Sum/lirui/project/H3K27ac/bowtie2/input.rmdup.bam  -f BAM -g hs -n SRRxxx -B -p 1e-4

cat SRRxxx_peaks.narrowPeak |awk -F "\t" '{print $1"\t"$4"\t"".""\t"$2"\t"$3"\t"".""\t"".""\t"".""\t"$4}' > /Sum/lirui/project/H3K27ac/gff/SRRxxx.gff

##5.call super-enhancers
python ROSE_main.py -g hg19 -i /Sum/lirui/project/H3K27ac/gff/SRRxxx.gff -r /Sum/lirui/project/H3K27ac/bowtie2/SRRxxx.rmdup.bam -c /Sum/lirui/project/H3K27ac/bowtie2/input.rmdup.bam -o /Sum/lirui/project/H3K27ac/SRRxxx_ROSE/ -s 12500 -t 2500
