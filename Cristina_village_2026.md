List files
```
===== andatsakala village =====
PhiX_S1_L001_R1_001.fastq.gz 
PhiX_S1_L001_R2_001.fastq.gz

===== mandena village =====
pool1 pool2
---- Subdirectory: mandena village/pool1/ ----PhiX_S1_L001_R1_001.fastq.gz 
PhiX_S1_L001_R2_001.fastq.gz
---- Subdirectory: mandena village/pool2/ ----PhiX_S1_L001_R1_001.fastq.gz 
PhiX_S1_L001_R2_001.fastq.gz

===== sarahandrano village =====
PhiX_L001-ds.a2d8ffbdc7f04d428212dae5984f48c5_sara.zip
```

```
#check file corruption: 
#!/bin/bash  
  
dirs=("andatsakala_village" "mandena_village" "sarahandrano_village")  
  
for d in "${dirs[@]}"; do  
echo "===== $d ====="  
  
find "$d" -type f | while read f; do  
echo "---- $f ----"  
  
if [[ "$f" == *.gz ]]; then  
zcat "$f" | head -n 20  
elif [[ "$f" == *.zip ]]; then  
unzip -p "$f" | head -n 20  
else  
head -n 20 "$f"  
fi  
  
echo  
done  
done


```


===== andatsakala village =====
~={orange}PhiX_S1_L001_R1_001.fastq.gz = looks ok possibly drop in quality at 3' end for the reads=~
~={orange}PhiX_S1_L001_R2_001.fastq.gz = looks like bad sequencing maybe? a lot of Ns.=~

===== mandena village =====
pool1 pool2
~={green}---- Subdirectory: mandena village/pool1/ ----=~
~={green}PhiX_S1_L001_R1_001.fastq.gz = good=~
~={green}PhiX_S1_L001_R2_001.fastq.gz = good=~

~={orange}---- Subdirectory: mandena village/pool2/ ----=~
~={orange}PhiX_S1_L001_R1_001.fastq.gz = looks ok possibly drop in quality at 3' end for the reads=~
~={orange}PhiX_S1_L001_R2_001.fastq.gz = looks like bad sequencing maybe?=~

```
zcat "mandena village/pool2/PhiX_S1_L001_R2_001.fastq.gz" | head -n 400

zcat "mandena village/pool2/PhiX_S1_L001_R2_001.fastq.gz" | awk 'NR%4==2' | head -100000 | grep -c "^N*$"
#224, so 224 reads out of 100k are all Ns. 

not all reads are all Ns, so lets try to demux and see what happens.
```
 

===== sarahandrano village =====
PhiX_L001-ds.a2d8ffbdc7f04d428212dae5984f48c5_sara.zip

```
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE  
unzip PhiX_L001-ds.a2d8ffbdc7f04d428212dae5984f48c5_sara.zip


```

~={green}---- sarahandrano village/PhiX_S1_L001_R1_001.fastq.gz ---- looks good=~
~={green}---- sarahandrano village/PhiX_S1_L001_R2_001.fastq.gz ---- looks good=~


## import, testing on sarahandrano first...
- cannot import like usual via EMP methods b/c we dont have barcodes.fastq, try this method

```
#create raw_reads dir and move the reads into it
mkdir raw_reads

#now we have: 
#sarahandrano_village/raw_reads 
#PhiX_S1_L001_R1_001.fastq.gz  
#PhiX_S1_L001_R2_001.fastq.gz

#rename the files to: 
forward.fastq.gz and reverse.fastq.gz

#!/bin/bash
#SBATCH --job-name=import_sarahandrano_village
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/sarahandrano_village

qiime tools import \
--type MultiplexedPairedEndBarcodeInSequence \
--input-path raw_reads \
--input-format MultiplexedPairedEndBarcodeInSequenceDirFmt \
--output-path sarahandrano_reads.qza

```

sbatch import_2.sh
Submitted batch job 24619415
Works!

# Extract barcodes in order to demultiplex and create a new file that lists forward and reverse barcodes ("tags") for each sample
```
#extract barcodes/tags by removing primer sequences from the barcodes fasta

awk '
/^>/ {tag=substr($0,2)}
/^\^/ {
  seq=$0
  gsub("\\^","",seq)
  sub("GTGYCAGCMGCCGCGGTAA.*","",seq)
  print tag"\t"seq
}' SmMa16S_forward_barcodes.fasta > forward_barcodes.tsv


awk '  
/^>/ {tag=substr($0,2)}  
/^\^/ {  
seq=$0  
gsub("\\^","",seq)  
sub("GGACTACNVGGGTWTCTAAT.*","",seq)  
print tag"\t"seq  
}' SmMa16S_reverse_barcodes.fasta > reverse_barcodes.tsv

awk '  
BEGIN{  
FS=OFS="\t"  
while((getline<"forward_barcodes.tsv")>0) f[$1]=$2  
while((getline<"reverse_barcodes.tsv")>0) r[$1]=$2  
}  
NR==1{print "sample-id","barcode-forward","barcode-reverse"; next}  
{  
print $1, f[$3], r[$4]  
}  
' metadata_microgale_microbiome_VS.txt > barcode_metadata.txt
```

# Then in excel i subset the barcode_metadata.txt file into each run: 
- also had to append numbers to duplicate sample names (e.g. ECs)
- Files: 
	- barcode_metadata_mandena1.txt
	- barcode_metadata_and.txt
	- barcode_metadata_mandena2.txt
	- barcode_metadata_sara.txt

## Demultiplex using the barcode_metadata files 

# sarahandrano_village
```
#!/bin/bash
#SBATCH --job-name=demux_sara
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/sarahandrano_village

qiime cutadapt demux-paired \
--i-seqs sarahandrano_village/sarahandrano_reads.qza \
--m-forward-barcodes-file barcode_metadata_sara.txt \
--m-forward-barcodes-column barcode-forward \
--m-reverse-barcodes-file barcode_metadata_sara.txt \
--m-reverse-barcodes-column barcode-reverse \
--p-anchor-forward-barcode \
--p-anchor-reverse-barcode \
--p-error-rate 0 \
--p-cores 4 \
--o-per-sample-sequences sarahandrano_village/demux/demux_sara.
qza \
--o-untrimmed-sequences sarahandrano_village/demux/unassigned_reads.qza

qiime demux summarize \
--i-data sarahandrano_village/demux/demux_sara.qza \
--o-visualization sarahandrano_village/demux/demux_sara.qzv
```
Submitted batch job 24620809


# Import the remaining runs

```
#!/bin/bash
#SBATCH --job-name=import_other_runs
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=amilan
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1

qiime tools import \
--type MultiplexedPairedEndBarcodeInSequence \
--input-path raw_reads \
--input-format MultiplexedPairedEndBarcodeInSequenceDirFmt \
--output-path mandena_pool1_reads.qza

#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2

qiime tools import \
--type MultiplexedPairedEndBarcodeInSequence \
--input-path raw_reads \
--input-format MultiplexedPairedEndBarcodeInSequenceDirFmt \
--output-path mandena_pool2_reads.qza

#andatsakala Village
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/andatsakala_village

qiime tools import \
--type MultiplexedPairedEndBarcodeInSequence \
--input-path raw_reads \
--input-format MultiplexedPairedEndBarcodeInSequenceDirFmt \
--output-path andatsakala_reads.qza

```
sbatch import.sh
Submitted batch job 24627480

# Demux remaining runs

```
#!/bin/bash
#SBATCH --job-name=demux_other_runs
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1

qiime cutadapt demux-paired \
--i-seqs mandena_pool1_reads.qza \
--m-forward-barcodes-file barcode_metadata_mandena1.txt \
--m-forward-barcodes-column barcode-forward \
--m-reverse-barcodes-file barcode_metadata_mandena1.txt \
--m-reverse-barcodes-column barcode-reverse \
--p-anchor-forward-barcode \
--p-anchor-reverse-barcode \
--p-error-rate 0 \
--p-cores 4 \
--o-per-sample-sequences demux/demux_mandena1.qza \
--o-untrimmed-sequences demux/unassigned_reads.qza

qiime demux summarize \
--i-data demux/demux_mandena1.qza \
--o-visualization demux/demux_mandena1.qzv


#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2

qiime cutadapt demux-paired \
--i-seqs mandena_pool2_reads.qza \
--m-forward-barcodes-file barcode_metadata_mandena2.txt \
--m-forward-barcodes-column barcode-forward \
--m-reverse-barcodes-file barcode_metadata_mandena2.txt \
--m-reverse-barcodes-column barcode-reverse \
--p-anchor-forward-barcode \
--p-anchor-reverse-barcode \
--p-error-rate 0 \
--p-cores 4 \
--o-per-sample-sequences demux/demux_mandena2.qza \
--o-untrimmed-sequences demux/unassigned_reads.qza

qiime demux summarize \
--i-data demux/demux_mandena2.qza \
--o-visualization demux/demux_mandena2.qzv

#andatsakala Village
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/andatsakala_village

qiime cutadapt demux-paired \
--i-seqs andatsakala_reads.qza \
--m-forward-barcodes-file barcode_metadata_and.txt \
--m-forward-barcodes-column barcode-forward \
--m-reverse-barcodes-file barcode_metadata_and.txt \
--m-reverse-barcodes-column barcode-reverse \
--p-anchor-forward-barcode \
--p-anchor-reverse-barcode \
--p-error-rate 0 \
--p-cores 4 \
--o-per-sample-sequences demux/demux_and.qza \
--o-untrimmed-sequences demux/unassigned_reads.qza

qiime demux summarize \
--i-data demux/demux_and.qza \
--o-visualization demux/demux_and.qzv
```

sbatch demux.sh
Submitted batch job 24627762

## Notes
Transferring the barcodes files, slurm, demux folders for all 3 villages NOT the reads because they're huge and you wont need them 
transferring the raw reads only for the sarah village since those were not originally unzipped before


## my suggestion for denoising 
Here is what i would do for denoising
make sure to edit file paths, email, etc
```
# denoise 

#!/bin/bash
#SBATCH --job-name=denoise_mand
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=EMAIL@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1/dada2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../demux/demux_mandena1.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 243 \
--p-trunc-len-r 242 \
--p-n-threads 6 \
--o-representative-sequences mandena1_seqs_dada2.qza \
--o-denoising-stats mandena1_dada2_stats.qza \
--o-table mandena1_table_dada2.qza

#Visualize the denoising results:
qiime metadata tabulate \
--m-input-file mandena1_dada2_stats.qza \
--o-visualization mandena1_dada2_stats.qzv

qiime feature-table summarize \
--i-table mandena1_table_dada2.qza \
--m-sample-metadata-file ../metadata/metadata.txt \
--o-visualization mandena1_table_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data mandena1_seqs_dada2.qza \
--o-visualization mandena1_seqs_dada2.qzv

#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../demux/demux_mandena2.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 243 \
--p-trunc-len-r 243 \
--p-n-threads 6 \
--o-representative-sequences mandena2_seqs_dada2.qza \
--o-denoising-stats mandena2_dada2_stats.qza \
--o-table mandena2_table_dada2.qza

#Visualize the denoising results:
qiime metadata tabulate \
--m-input-file mandena2_dada2_stats.qza \
--o-visualization mandena2_dada2_stats.qzv

qiime feature-table summarize \
--i-table mandena2_table_dada2.qza \
--m-sample-metadata-file ../metadata/metadata.txt \
--o-visualization mandena2_table_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data mandena2_seqs_dada2.qza \
--o-visualization mandena2_seqs_dada2.qzv
```



## I accidentally didn't include the pos/neg/ec controls in mandena pool1, plus the pos/neg/ec were duplicate names between pool 1 and 2 (because the same samples were used across both pools), so that needs to be fixed. i updated the samples names in the metadata and barcodes_metadata files for each pool. re-demultiplexing now: 

```
#!/bin/bash
#SBATCH --job-name=demux_mand_redo
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1

qiime cutadapt demux-paired \
--i-seqs mandena_pool1_reads.qza \
--m-forward-barcodes-file barcode_metadata_mandena1.txt \
--m-forward-barcodes-column barcode-forward \
--m-reverse-barcodes-file barcode_metadata_mandena1.txt \
--m-reverse-barcodes-column barcode-reverse \
--p-anchor-forward-barcode \
--p-anchor-reverse-barcode \
--p-error-rate 0 \
--p-cores 4 \
--o-per-sample-sequences demux/demux_mandena1.qza \
--o-untrimmed-sequences demux/unassigned_reads.qza

qiime demux summarize \
--i-data demux/demux_mandena1.qza \
--o-visualization demux/demux_mandena1.qzv


#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2

qiime cutadapt demux-paired \
--i-seqs mandena_pool2_reads.qza \
--m-forward-barcodes-file barcode_metadata_mandena2.txt \
--m-forward-barcodes-column barcode-forward \
--m-reverse-barcodes-file barcode_metadata_mandena2.txt \
--m-reverse-barcodes-column barcode-reverse \
--p-anchor-forward-barcode \
--p-anchor-reverse-barcode \
--p-error-rate 0 \
--p-cores 4 \
--o-per-sample-sequences demux/demux_mandena2.qza \
--o-untrimmed-sequences demux/unassigned_reads.qza

qiime demux summarize \
--i-data demux/demux_mandena2.qza \
--o-visualization demux/demux_mandena2.qzv
```
Submitted batch job 24756081


## Trim out primers (515f and 806r) from the reads

```
#!/bin/bash
#SBATCH --job-name=trim_mand
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1

qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux/demux_mandena1.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--p-adapter-r TTACCGCGGCMGCTGYCAC \
--p-adapter-f ATTAGAWACCCVNGTAGTCC \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--o-trimmed-sequences trimmed_reads/mandena1_trimmed_reads.qza \
--verbose


#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2

qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux/demux_mandena2.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--p-adapter-r TTACCGCGGCMGCTGYCAC \
--p-adapter-f ATTAGAWACCCVNGTAGTCC \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--o-trimmed-sequences trimmed_reads/mandena2_trimmed_reads.qza \
--verbose

#visualize trimmed reads 
qiime demux summarize \
--i-data trimmed_reads/mandena2_trimmed_reads.qza \
--o-visualization trimmed_reads/mandena2_trimmed_reads.qzv

cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1

#visualize trimmed reads 
qiime demux summarize \
--i-data trimmed_reads/mandena1_trimmed_reads.qza \
--o-visualization trimmed_reads/mandena1_trimmed_reads.qzv
```
Submitted batch job 24785474

# Results
the output for this looks good!! e.g. that the primers were trimmed from 356 samples which is what we expect. 

see the cutadapt_output.txt for that info
data_shareable/Cristina_village_2026/mandena_village/cutadapt_output.txt


# Denoise the trimmed sequences using dada2
```
# denoise 

#!/bin/bash
#SBATCH --job-name=denoise_mand
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=amilan
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1/dada2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../trimmed_reads/mandena1_trimmed_reads.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 223 \
--p-trunc-len-r 222 \
--p-n-threads 6 \
--o-representative-sequences mandena1_seqs_dada2.qza \
--o-denoising-stats mandena1_dada2_stats.qza \
--o-table mandena1_table_dada2.qza

#Visualize the denoising results:
qiime metadata tabulate \
--m-input-file mandena1_dada2_stats.qza \
--o-visualization mandena1_dada2_stats.qzv

qiime feature-table summarize \
--i-table mandena1_table_dada2.qza \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt \
--o-visualization mandena1_table_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data mandena1_seqs_dada2.qza \
--o-visualization mandena1_seqs_dada2.qzv

#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2/dada2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../trimmed_reads/mandena2_trimmed_reads.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 225 \
--p-trunc-len-r 221 \
--p-n-threads 6 \
--o-representative-sequences mandena2_seqs_dada2.qza \
--o-denoising-stats mandena2_dada2_stats.qza \
--o-table mandena2_table_dada2.qza

#Visualize the denoising results:
qiime metadata tabulate \
--m-input-file mandena2_dada2_stats.qza \
--o-visualization mandena2_dada2_stats.qzv

qiime feature-table summarize \
--i-table mandena2_table_dada2.qza \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt \
--o-visualization mandena2_table_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data mandena2_seqs_dada2.qza \
--o-visualization mandena2_seqs_dada2.qzv
```

Submitted batch job 24794780

these denoising stats are terrible, losing most of the reads at the quality filtering step, probably because of the quality distribution of the end of the reads
-renamed these dada2 folders to dada2_bad, but did not transfer to sharepoint

# mand pool 2 is the limiting run, because its quality is poorer than mand1, so we have to use truncate params that match pool2

--p-trunc-len-f 188 \
--p-trunc-len-r 181 \

## Try these truncation parameters: 

```
#!/bin/bash
#SBATCH --job-name=denoise_mand1
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=amilan
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 1
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1/dada2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../trimmed_reads/mandena1_trimmed_reads.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 188 \
--p-trunc-len-r 181 \
--p-n-threads 6 \
--o-representative-sequences mandena1_seqs_dada2.qza \
--o-denoising-stats mandena1_dada2_stats.qza \
--o-table mandena1_table_dada2.qza

#Visualize the denoising results:
qiime metadata tabulate \
--m-input-file mandena1_dada2_stats.qza \
--o-visualization mandena1_dada2_stats.qzv

qiime feature-table summarize \
--i-table mandena1_table_dada2.qza \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt \
--o-visualization mandena1_table_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data mandena1_seqs_dada2.qza \
--o-visualization mandena1_seqs_dada2.qzv
```
sbatch denoise_mand1.sh
Submitted batch job 24817713

# These look WAY better. repeat for pool 2 

```
#!/bin/bash
#SBATCH --job-name=denoise_mand2
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=amilan
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

#Mandena Village Pool 2
cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2/dada2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../trimmed_reads/mandena2_trimmed_reads.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 188 \
--p-trunc-len-r 181 \
--p-n-threads 6 \
--o-representative-sequences mandena2_seqs_dada2.qza \
--o-denoising-stats mandena2_dada2_stats.qza \
--o-table mandena2_table_dada2.qza

#Visualize the denoising results:
qiime metadata tabulate \
--m-input-file mandena2_dada2_stats.qza \
--o-visualization mandena2_dada2_stats.qzv

qiime feature-table summarize \
--i-table mandena2_table_dada2.qza \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt \
--o-visualization mandena2_table_dada2.qzv

qiime feature-table tabulate-seqs \
--i-data mandena2_seqs_dada2.qza \
--o-visualization mandena2_seqs_dada2.qzv
```
sbatch denoise_mand2.sh
Submitted batch job 24811441

## Denoising results

- There are contamination long amplicons in the representative sequences file. You will need to use the amplicon removal/filtering code from homework 2 to remove anything larger than 300bp)
    
- Read retention after denoising is great! Looks like the average for pool1 real samples (not controls) is 74% and for pool2 real samples (not controls) is 79%. Amazing!

```
# trim out long amplicons from each pool

#pool1

cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool1/dada2

qiime feature-table filter-seqs \
--i-data mandena1_seqs_dada2.qza \
--m-metadata-file mandena1_seqs_dada2.qza \
--p-where 'length(sequence) < 256' \
--o-filtered-data mandena1_seqs_dada2_filtered256.qza

qiime feature-table tabulate-seqs \
--i-data mandena1_seqs_dada2_filtered256.qza \
--o-visualization mandena1_seqs_dada2_filtered256.qzv

qiime feature-table filter-features \
--i-table mandena1_table_dada2.qza \
--m-metadata-file mandena1_seqs_dada2_filtered256.qza \
--o-filtered-table mandena1_table_dada2_filtered256.qza
  
qiime feature-table summarize \
--i-table mandena1_table_dada2_filtered256.qza \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt \
--o-visualization mandena1_table_dada2_filtered256.qzv
```


```
# repeat the filtering for pool2


cd /scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/pool2/dada2

qiime feature-table filter-seqs \
--i-data mandena2_seqs_dada2.qza \
--m-metadata-file mandena2_seqs_dada2.qza \
--p-where 'length(sequence) < 256' \
--o-filtered-data mandena2_seqs_dada2_filtered256.qza

qiime feature-table tabulate-seqs \
--i-data mandena2_seqs_dada2_filtered256.qza \
--o-visualization mandena2_seqs_dada2_filtered256.qzv

qiime feature-table filter-features \
--i-table mandena2_table_dada2.qza \
--m-metadata-file mandena2_seqs_dada2_filtered256.qza \
--o-filtered-table mandena2_table_dada2_filtered256.qza
  
qiime feature-table summarize \
--i-table mandena2_table_dada2_filtered256.qza \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt \
--o-visualization mandena2_table_dada2_filtered256.qzv
```

```
# Merge the seqs file and the table, into 1 seqs file and 1 table. 

# make a new directory for the merged data
mkdir merged_data
/scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/merged_data/

#then make a new dada2 directory for the merged files
mkdir dada2_merged
/scratch/alpine/lindsval@colostate.edu/Cristina_village_2026/mandena_village/merged_data/dada2_merged

#then run this once you are in the new merged_data/dada2_merged directory

qiime feature-table merge \
--i-tables ../../pool1/dada2/mandena1_table_dada2_filtered256.qza \
--i-tables ../../pool2/dada2/mandena2_table_dada2_filtered256.qza \
--o-merged-table merged_table.qza  
  
qiime feature-table summarize \
--i-table merged_table.qza \
--o-visualization merged_table.qzv \
--m-sample-metadata-file ../../metadata/metadata_microgale_microbiome_VS_new.txt

qiime feature-table merge-seqs \
--i-data ../../pool1/dada2/mandena1_seqs_dada2_filtered256.qza \
--i-data ../../pool2/dada2/mandena2_seqs_dada2_filtered256.qza \
--o-merged-data merged_seqs.qza  
  
qiime feature-table tabulate-seqs \
--i-data merged_seqs.qza \
--o-visualization merged_seqs.qzv
```

## now use these files to assign taxonomy, run the tree and do core metrics. this all should run the same way as the class friday tutorials. the rest of the commands should be run in the merged_data folders

so make sure to create the taxonomy, tree folders in the mandena_village/merged_data folder.

