# Chikashigeetal
in-house scripts for Chikashige et al.

This is a set of in-house scripts for Chikashige et al.

Before running scripts, fetch required files indicated in the files 01_NGS_reads_preparation.txt from SRA and 02_pombase_data_preparation.txt from Pombase, and set HOMEDIR in the scripts 03-11.

To perform analyses, run the scripts below.

cd $HOMEDIR
bash 03_fastqdump_cutadapt_trim.sh
Rscript 04_uORFsearch.R
Rscript 05_spliced_transcripts_preparation.R
bash 06_bwa_bedtoolsgenomecov.sh
Rscript 07_RPKM_calculation_for_each_sample.R
Rscript 08_weighted_average_calculation.R
Rscript 09_data_summarization.R
Rscript 10_figure_preparation1.R
Rscript 11_figure_preparation2.R


Contact:
Hiroaki Kato
hkato@med.shimane-u.ac.jp
