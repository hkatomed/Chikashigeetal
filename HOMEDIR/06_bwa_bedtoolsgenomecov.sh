
# Specify the name of home directory before run.
HOMEDIR=XXX


cd $HOMEDIR
cd reference_transcriptome
bwa index sp_spliced_transcripts.fasta

cd $HOMEDIR
cd TXT
SEQIDS=$(<seqIDs3.txt)

cd $HOMEDIR
mkdir bwa_aln
mkdir bwa_samse
for SEQID in $SEQIDS
do
bwa aln -t 32 reference_transcriptome/sp_spliced_transcripts.fasta \
fastq3/cutadapt_trim25_${SEQID}.fastq > bwa_aln/${SEQID}.sai
bwa samse reference_transcriptome/sp_spliced_transcripts.fasta \
bwa_aln/${SEQID}.sai fastq3/cutadapt_trim25_${SEQID}.fastq > bwa_samse/${SEQID}.sam
done

cd $HOMEDIR
mkdir flagstat_2
mkdir idxstats_2
mkdir mappedReads_2
mkdir bedGraph_5end_2
mkdir RPKM_2
cd bwa_samse
for SEQID in $SEQIDS
do
sambamba_v0.6.6 view -h -f=bam -t=8 -S ${SEQID}.sam -o ${SEQID}.bam
sambamba_v0.6.6 sort -t=8 ${SEQID}.bam -o ${SEQID}.sorted.bam
sambamba_v0.6.6 view -F "unmapped" -h -f=sam -t=8 ${SEQID}.sorted.bam -o ${SEQID}.unmapped.sam
sambamba_v0.6.6 flagstat -t=8 ${SEQID}.sorted.bam > \
../flagstat_2/${SEQID}_flagstat.txt
samtools idxstats ${SEQID}.sorted.bam > \
../idxstats_2/${SEQID}_idxstats.txt
awk 'NR==1 {print $1}' ../flagstat_2/${SEQID}_flagstat.txt > \
../totalReads_2/${SEQID}_totalReads.txt
awk 'NR==5 {print $1}' ../flagstat_2/${SEQID}_flagstat.txt > \
../mappedReads_2/${SEQID}_mappedReads.txt
bedtools genomecov -ibam ${SEQID}.sorted.bam \
-d -5 -strand + -g ../genome_info/sp_spliced_transcriptome_info.txt \
> ../bedGraph_5end_2/${SEQID}_Sense.bedGraph
done

cd $HOMEDIR
for SEQID in $SEQIDS
do
Rscript Rscripts/getShiftedRPKM.R ${SEQID}
done
