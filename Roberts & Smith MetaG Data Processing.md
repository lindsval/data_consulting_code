
# Directory Structure 
- per sample directories with: raw_reads, processed_reads, fastqc, assembly subdirectories
- CURC software https://curc.readthedocs.io/en/latest/software/curc_provided_software.html

```

# unzip the files

#!/bin/bash
#SBATCH --job-name=gunzip_rob_metag
#SBATCH --partition=amilan  
#SBATCH --qos=normal  
#SBATCH --time=04:00:00  
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=gunzip_%j.out  
#SBATCH --error=gunzip_%j.err  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=lindsval@colostate.edu  
  
# Go to your working directory  
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG  
  
# Load pigz  
module purge  
module load pigz  
  
# Decompress all .fastq.gz files in parallel  
pigz -d -p $SLURM_CPUS_PER_TASK *.fastq.gz

```
sbatch 01_unzip.sh 
Submitted batch job 24697291

## Rename the files to their sample names
- see /Users/valerielindstrom/Documents/PostDoc/data_consulting/roberts_metagenomics/metadata/20250926_CCE_Metagenomics_Metadata_VS.xlsx
- i updated the file names to be more informative, they follow this format: [treatment]_[soilType]_[timepoint]_[plot]
- example:  Control_BulkSoil_Post_10
- previous samples names were not informative enough (e.g. Control_BulkSoil_Post_10 was just A10B, and the raw reads file name was A10B_cleaned_S1_R1_001.fastq.gz
	- for assembly and co-assembly having more informative names will help)

for renaming, i just used the ```mv``` function, see "rename_ files" tab in the metadata file. and then i copy and pasted that in the terminal. That code only produce one error, since there was a duplicate sample name for A38R, however there was a r1 and r2 for A38R and A38B
A38B_cleaned_S25_R1_001.fastq.gz and A38B_cleaned_S25_R2_001.fastq.gz
A38R_cleaned_S26_R1_001.fastq.gz and A38R_cleaned_S26_R2_001.fastq.gz

so in the metadata i just changed the first instance of A38R to A38B, since the clear patter was that the B samples always come before the R samples, and that the raw data files indicated that a38b was sample 25 and a38r was sample 26. Good to go!


## create directories for each sample
take that new sample list from the metadata file (all the samples listed in the sample_name column) and paste into a new .txt file, call it sample_list.txt

/Users/valerielindstrom/Documents/PostDoc/data_consulting/roberts_metagenomics/metadata/sample_list.txt

```
while read sample; do  
mkdir -p "$sample"  
done < sample_list.txt

dos2unix sample_list.txt

#create subdirectories
for d in */; do  
mkdir -p "${d}raw_reads"  
done

for d in */; do  
mkdir -p "${d}processed_reads"  
done

for d in */; do  
mkdir -p "${d}assembly"  
done

#move files into their sample/raw_reads folder
for f in *_R1.fastq *_R2.fastq; do  
sample=${f%_R*.fastq}  
mv "$f" "$sample/raw_reads/"  
done

#make fastq directories
for d in */; do  
mkdir -p "${d}fastqc"  
done

for d in */; do  
mkdir -p "${d}fastqc/raw" "${d}fastqc/trimmed"  
done
```
## Perform fastqc on the raw files (FastQC v 0.11.9)

```
#!/bin/bash  
#SBATCH --job-name=fastqc  
#SBATCH --partition=amilan  
#SBATCH --qos=normal  
#SBATCH --time=08:00:00  
#SBATCH --cpus-per-task=16  
#SBATCH --mem=64G  
#SBATCH --output=fastqc_%j.out  
#SBATCH --error=fastqc_%j.err  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=lindsval@colostate.edu  
  
module load fastqc  
  
# Go to your working directory  
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG

while read sample; do
    RAW_DIR="${sample}/raw_reads"
    OUT_DIR="${sample}/fastqc/raw"

    for fq in "$RAW_DIR"/*_R[12].fastq*; do
        [ -e "$fq" ] || continue
        fastqc -t 16 "$fq" -o "$OUT_DIR"
    done
done < sample_list.txt
```
sbatch 02_fastqc_raw.sh
Submitted batch job 24784019

make sure this worked for all samples (so it should print 88 samples) 

```
wc -l sample_list.txt
#88 samples total
```

```
count=0  
while read sample; do  
compgen -G "${sample}/fastq/raw/*R1_fastqc.html" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#88

count=0  
while read sample; do  
compgen -G "${sample}/fastq/raw/*R2_fastqc.html" > /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count

```
Great! fastQC ran on all of our samples. 

## Now check the read quality using MultiQC (v.1.19)
https://github.com/MultiQC/MultiQC

```
#install Multiqc

acompile --ntasks=4 
module load anaconda

#make sure youre chanels are set correctly: 
conda config --add channels bioconda
conda config --add channels conda-forge

# alternative to running that code, you can edit the .condarc file or the environment.yaml to make sure they are in that order

conda create -n multiqc
conda activate multiqc
conda install multiqc
multiqc --version #this should print v1.19
```

```
#use multiqc to generate one report for all fastqc_data.txt files
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG

multiqc . -o multiqc_report
```
Output: 
/// MultiQC 🔍 | v1.19| multiqc | MultiQC Version v1.33 now available!| 
multiqc | Search path : /gpfs/alpine1/scratch/lindsval@colostate.edu/roberts_soils_metaG| 

searching |  100% 726/726 
fastqc | Found 176 reports
multiqc | Report : multiqc_report/multiqc_report.html
multiqc | Data : multiqc_report/multiqc_data
multiqc | Flat-image plots used. Disable with '--interactive'. See docs.
multiqc | MultiQC complete

take a look at the output .html file:///Users/valerielindstrom/Downloads/multiqc_report.html
YouTube video on the outputs: https://www.youtube.com/watch?v=qPbIlO_KWN0 


## Fastqc results from the raw files: 
- read quality is fantastic, and the negative controls contain little reads and the positive controls look good too. 

## Install Sickle on Alpine

```
# https://github.com/najoshi/sickle
# version 1.33

acompile --ntasks=4 
module load anaconda
conda create -n sickle-trim
conda activate sickle-trim
conda install sickle-trim
```

# Run sickle on the raw reads, we will use the sickle default of 20 quality score for trimming

```
#!/bin/bash

while read -r element
do
  (
    cd "/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/$element/raw_reads" || exit 1

    sickle pe \
      -f "${element}_R1.fastq" \
      -r "${element}_R2.fastq" \
      -t sanger \
      -o "${element}_R1_trimmed.fastq" \
      -p "${element}_R2_trimmed.fastq" \
      -s "${element}_singles.fastq"

    rm -f "${element}_singles.fastq"

    mv *trimmed.fastq ../processed_reads/
  ) &

  # 4 parallel jobs
  if [[ $(jobs -r -p | wc -l) -ge 4 ]]; then
    wait
  fi

done < "$1"

wait

#save
03_trim_sickle_loop.sh

#!/bin/bash
#SBATCH --job-name=sickle 
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm_output/sickle_%j.out
#SBATCH --error=slurm_output/sickle_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu

module load anaconda
conda activate sickle-trim
bash 03_trim_sickle_loop.sh ../sample_list.txt

#save
03_trim_sickle.sh
```
sbatch 03_trim_sickle.sh
Submitted batch job 24812106, finished in a few hours



## Run bbduk (from the BBTools Package, bbtools/v39.01) to remove adapters
https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

```
#loop

#!/bin/bash

while read -r element
do
  (
    cd "/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/$element/processed_reads" || exit 1

	bbduk.sh \
    threads=4 \
    overwrite=t \
    in1="${element}_R1_trimmed.fastq" \
    in2="${element}_R2_trimmed.fastq" \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe \
    tbo \
    ref=/curc/sw/install/bio/bbtools/bbmap/resources/adapters.fa \
    out1="${element}_R1_bbduktrimmed.fastq \
    out2="${element}_R2_bbduktrimmed.fastq

  ) &

  # 4 parallel jobs
  if [[ $(jobs -r -p | wc -l) -ge 4 ]]; then
    wait
  fi

done < "$1"

wait

```
04_bbduk_loop.sh

```
#!/bin/bash
#SBATCH --job-name=bbduk
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm_output/bbduk_%j.out
#SBATCH --error=slurm_output/bbduk_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu

module load bbtools
bash 04_bbduk_loop.sh ../sample_list.txt
```

sbatch 04_bbduk.sh
Submitted batch job 24977428

## Notes on bbduk parameters: 
- **ktrim=r means k-mer based trimming**: In ktrim=r mode, once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; this is the normal mode for adapter trimming. this goes hand-in-hand with using a reference (this line ref=/opt/bbtools/bbmap/resources/adapters.fa)

- bbduk & memory: Most operations such as adapter-trimming and quality-trimming need only a tiny amount of memory.

- hdist is hamming distance, 1 is good, this allows one mismatch.

- flags “tbo”, which specifies to also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences), and “tpe”, which specifies to trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).


```
#count the bbduk trimmed files: 
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG

count=0  
while read sample; do  
compgen -G "${sample}/processed_reads/*R1_bbduktrimmed.fastq" > /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count
#88

count=0  
while read sample; do  
compgen -G "${sample}/processed_reads/*R2_bbduktrimmed.fastq" >  /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count
#88 yay!
```


## redo fastqc and make sure we pass checks and adapters were removed. 

```
#!/bin/bash  
#SBATCH --job-name=fastqc_trimmed  
#SBATCH --partition=amilan  
#SBATCH --qos=normal  
#SBATCH --time=08:00:00  
#SBATCH --cpus-per-task=16  
#SBATCH --mem=64G  
#SBATCH --output=fastqc_%j.out  
#SBATCH --error=fastqc_%j.err  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=lindsval@colostate.edu  
  
module load fastqc  
  
# Go to your working directory  
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG

while read sample; do
    TRIMMED_DIR="${sample}/processed_reads"
    OUT_DIR="${sample}/fastqc/trimmed"

    for fq in "$TRIMMED_DIR"/*_R[12]_bbduktrimmed.fastq; do
        [ -e "$fq" ] || continue
        fastqc -t 16 "$fq" -o "$OUT_DIR"
    done
done < sample_list.txt
```
Submitted batch job 24986389

```
count=0  
while read sample; do  
compgen -G "${sample}/fastqc/trimmed/*R1_bbduktrimmed_fastqc.html" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#88


count=0  
while read sample; do  
compgen -G "${sample}/fastqc/trimmed/*R2_bbduktrimmed_fastqc.html" > /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count
#88
```

## multiqc

```
#use multiqc to generate one report for all fastqc_data.txt files
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG
module load anaconda
conda activate multiqc


multiqc \
/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/*/fastqc/trimmed/*bbduktrimmed*_fastqc.zip \
--outdir /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_trimmed \
--filename trimmed_multiqc_report.html

```

 100% 176/176 | fastqc | Found 176 reports| multiqc | Report : ../../../../../scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_trimmed/trimmed_multiqc_report.html| 
 multiqc | 
 Data : ../../../../../scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_trimmed/trimmed_multiqc_report_data| 
 multiqc | Flat-image plots used. Disable with '--interactive'. See docs.| 
 multiqc | MultiQC complete

### Questions/Notes:
1. How many reads did we lose?
	1. looks like for the samples, we lost around 800k reads per sample. but sample reads are still around 20mil 
2. were adapters actually removed? 
	1. Yes!
3. Does the number of R1 and R2 reads match? 
	1. yes they do match, some of the controls are less like around 132 and 145bp
4. Note the quality metrics- have they improved? 
	1. quality looks good
	2. Most samples have greater than 10 million reads per read/sample
5. the per base sequence count is still a little wonky for the first few base pairs; sounds like its a normal result of library prep methods which use tagmentation. see [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html) and kit with tagmentations protocol [here](). 

file:///Users/valerielindstrom/Downloads/trimmed_multiqc_report%20(1).html 

# zip or delete raw reads
At this point, we will only proceed with the trimmed reads. As such, let's either zip the raw reads to save space or delete them from the working directory.

```
#!/bin/bash
#SBATCH --job-name=zip_reads
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=23:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_output/zip_%j.out
#SBATCH --error=slurm_output/zip_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu

module load pigz

while read SAMPLE; do
find /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/${SAMPLE}/raw_reads \ 
-type f ! -name "*.gz" -exec pigz {} +
done < /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt
```
sbatch 06_zip_raw_reads.sh
Submitted batch job 25178913

```
#check they were all zipped

count=0  
while read sample; do  
compgen -G "${sample}/raw_reads/*R1*.gz" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#88


count=0  
while read sample; do  
compgen -G "${sample}/raw_reads/*R2*.gz" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#88

```


## install MEGAHIT 

```
acompile --ntasks=4 --time=03:00:00
module load anaconda
conda create -n megahit
conda activate megahit
conda install -c bioconda megahit
megahit -v
#this should print: MEGAHIT v1.2.9
```

## Assemble trimmed reads using MEGAHIT (v1.2.9)

```
#!/bin/bash

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# number of samples to run at once
MAX_JOBS=5

while read SAMPLE; do
  (
    R1="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R1_bbduktrimmed.fastq"
    R2="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R2_bbduktrimmed.fastq"
    OUTDIR="${BASE_DIR}/${SAMPLE}/assembly/megahit_out"

    # check files exist
    if [[ -f "$R1" && -f "$R2" ]]; then
      echo "Running MEGAHIT for $SAMPLE"

      megahit \
        -1 "$R1" \
        -2 "$R2" \
        --k-min 31 --k-max 121 --k-step 10 \
        -m 0.4 \
        -t 10 \
        -o "$OUTDIR"

    else
      echo "Missing reads for $SAMPLE" >&2
    fi
  ) &

  # limit number of concurrent jobs
  if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
    wait -n
  fi

done < "$SAMPLE_LIST"

wait
```
07_megahit_individual_assembly_loop.sh

```

#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --nodes=1
#SBATCH --cpus-per-task=55
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/megahit_%j.out
#SBATCH --error=slurm_output/megahit_%j.err


module load anaconda
conda activate megahit

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/slurm

bash 07_megahit_individual_assembly_loop.sh 
```
07_megahit_individual_assembly.sh
Submitted batch job 25224839


```
#test assembly stats on one metaG
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/

/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts/contig_stats.pl -i /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/megahit_out_testing_on_Drought_Rhizo_Post_11/final.contigs.fa -o Drought_Rhizo_Post_11_final.contigs_STATS
```