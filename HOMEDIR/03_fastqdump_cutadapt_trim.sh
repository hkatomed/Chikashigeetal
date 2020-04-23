
# Specify the name of home directory before run.
HOMEDIR=XXX

cd $HOMEDIR
cd TXT
SEQIDS=$(<seqIDs3.txt)

cd $HOMEDIR
cd SRA
for SEQID in $SEQIDS
do
fastq-dump ${SEQID} > fastqdumplog_${SEQID}.txt
done
cat fastqdumplog_* > fastqdumplog.txt
rm fastqdumplog_*
mv *.fastq ../fastq1

ADAPT=NNNNTGGAATTCTCGGGTGCCAAGG
cd $HOMEDIR
mkdir fastq3
cd fastq1
for SEQID in $SEQIDS
do
cat ${SEQID}.fastq | paste - - - - | LC_ALL=C sort -t$'\t' -k2,2 -u | \
tr "\t" "\n" > unique_${SEQID}.fastq
cutadapt -a $ADAPT -j 8 -m 29 -M 35 --discard-untrimmed \
-o cutadapt_${SEQID}.fastq unique_${SEQID}.fastq \
> cutadaptlog_${SEQID}.txt
fastx_trimmer -f 5 -l 29 -i cutadapt_${SEQID}.fastq -o cutadapt_trim25_${SEQID}.fastq
done
cat cutadaptlog_* > cutadaptlog.txt
rm cutadaptlog_*
mv cutadapt_trim25_* ../fastq3
