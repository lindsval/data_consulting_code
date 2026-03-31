Part 1: Data upload through trimming reads, and removing adapters; reruning final QC check

Training Doc: https://docs.google.com/presentation/d/1hXJD7P0WjOE-iaKu9HhdtIVDVLCY8hIjD71KNBJ04WY/edit?slide=id.g3cd05613a6c_0_62#slide=id.g3cd05613a6c_0_62

----------------------------------------------------------------

## Steps
 Step 1: upload your file(s) to alpine
 
 Step 2: create a working directory and unzip the files 
 
 Step 3: Rename the fastq files
 
 Step 4: create directories for each sample you chose. Within the sample directories, create the follow subdirectories:  raw_reads, processed_reads, assembly, and fastqc. within the fastqc directory, create a 'raw' and 'trimmed' directory
 
Step 5: move fastq files into their sample/raw_reads folder

Step 6: Perform fastqc on the raw files

Step 7: Now check the read quality using MultiQC (v.1.19)

Step 8: Install Sickle on Alpine

Step 9: Run sickle on the raw reads

Step 10: Run bbduk (from the BBTools Package, bbtools/v39.01) to remove adapters

Step 11: redo fastqc and make sure we pass checks and adapters were removed

----------------------------------------------------------------


```
# Step 2: create a working directory, and go into it and unzip the files 

cd /scratch/alpine/$USER$@colostate.edu/DIR_NAME

### load nodes
ainteractive --ntasks=4 
module purge

### unzip the file(s) using pigz
module load pigz   
pigz -d -p *.fastq.gz

```

## Rename the files to their sample names
- see /Users/valerielindstrom/Documents/PostDoc/data_consulting/roberts_metagenomics/metadata/20250926_CCE_Metagenomics_Metadata_VS.xlsx
- i updated the file names to be more informative, they follow this format: [treatment]_[soilType]_[timepoint]_[plot]
- example:  Control_BulkSoil_Post_10
	- Control_BulkSoil_Post_10 was just A10B, and the raw reads file name was A10B_cleaned_S1_R1_001.fastq.gz
	- for assembly and co-assembly having more informative names will help

for renaming, i just used the ```mv``` function in the terminal, see "rename_ files" tab in the metadata excel file. and then i copy and pasted that in the terminal. 

So for which ever samples you chose, find the corresponding new name and rename the samples

Example: if you chose sample A10B, then use this to rename the r1 and r2 files: 


```
## Step 3: Rename the fastq files

# sample 1 to rename
mv A11B_cleaned_S3_R1_001.fastq Drought_BulkSoil_Post_11_R1.fastq
mv A11B_cleaned_S3_R2_001.fastq Drought_BulkSoil_Post_11_R2.fastq

# sample 2 to rename
add the commands for the second sample

```

## create directories for each sample
take the sample_list.txt file (make sure you update the sample names in this file if you chose 2 different files than the two i put in there as an example) and use to loop over each sample so it will create these directories and subdirectories
- place the sample_list.txt file in the main directory: cd /scratch/alpine/$USER$@colostate.edu/DIR_NAME

## Upload the sample_list.txt file to alpine

```
## Step 4: Create directories for each sample, based on the sample_list.txt file
while read sample; do  
mkdir -p "$sample"  
done < sample_list.txt



## create subdirectories
for d in */; do  
mkdir -p "${d}raw_reads"  
done

for d in */; do  
mkdir -p "${d}processed_reads"  
done

for d in */; do  
mkdir -p "${d}assembly"  
done

#make fastq directories
for d in */; do  
mkdir -p "${d}fastqc"  
done

for d in */; do  
mkdir -p "${d}fastqc/raw" "${d}fastqc/trimmed"  
done
```


## Move raw fastq files into the raw_reads folders
```
#### Step 5: move fastq files into their sample/raw_reads folder
for f in *_R1.fastq *_R2.fastq; do  
sample=${f%_R*.fastq}  
mv "$f" "$sample/raw_reads/"  
done
```
## Step 6: Perform fastqc on the raw files (FastQC v 0.11.9)

- create a slurm directory and then create this .sh file to be able to submit this as a job
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
#SBATCH --mail-user=$USER@colostate.edu  
  
module load fastqc  
  
# Go to your working directory  
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG

module load fastqc

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


make sure this worked for all samples (so it should print # samples) 
should be 2 if thats the number of samples you picked
```
wc -l sample_list.txt
#  samples total
```

```
count=0  
while read sample; do  
compgen -G "${sample}/fastq/raw/*R1_fastqc.html" > /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count


count=0  
while read sample; do  
compgen -G "${sample}/fastq/raw/*R2_fastqc.html" > /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count

```
did fastQC ran on all of our samples?? 

## Step 7: Now check the read quality using MultiQC (v.1.19)
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


## Step 8: Install Sickle on Alpine

```
# https://github.com/najoshi/sickle
# version 1.33

acompile --ntasks=4 
module load anaconda
conda create -n sickle-trim
conda activate sickle-trim
conda install sickle-trim
```

# Step 9: Run sickle on the raw reads, we will use the sickle default of 20 quality score for trimming

the next part is where looping comes into play. this means you will have to .sh scripts, one is for loopping and one is for submitting the job. 

Create the looping file: 03_trim_sickle_loop.sh

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
```
save the file
03_trim_sickle_loop.sh


Create the job submission file: make sure the sample list is in the main directory and that you submit this file in the slurm directory. eg: 
main directory: cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG 
- the sample_list.txt should be here
slurm directory: cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/slurm
- the 03_trim_sickle_loop.sh and 03_trim_sickle.sh should be here

```
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
#SBATCH --mail-user=$USER$@colostate.edu

module load anaconda
conda activate sickle-trim
bash 03_trim_sickle_loop.sh ../sample_list.txt

#save
03_trim_sickle.sh
```
sbatch 03_trim_sickle.sh


## Step 10: Run bbduk (from the BBTools Package, bbtools/v39.01) to remove adapters
https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/

```
#this it the loop

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

## Notes on bbduk parameters: 
- **ktrim=r means k-mer based trimming**: In ktrim=r mode, once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; this is the normal mode for adapter trimming. this goes hand-in-hand with using a reference (this line ref=...adapters.fa)

- hdist is hamming distance, 1 is good, this allows one mismatch.

- flags “tbo”, which specifies to also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences), and “tpe”, which specifies to trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).


```
#count the bbduk trimmed files to make sure it ran on all of them: 
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


## Step 11: redo fastqc and make sure we pass checks and adapters were removed. 

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
#use multiqc to generate one report for all fastqc_data.txt files from trimmed reads

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG
module load anaconda
conda activate multiqc


multiqc \
/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/*/fastqc/trimmed/*bbduktrimmed*_fastqc.zip \
--outdir /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_trimmed \
--filename trimmed_multiqc_report.html

```

### Output 
 100% 176/176 | fastqc | Found 176 reports| multiqc | Report : ../../../../../scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_trimmed/trimmed_multiqc_report.html| 
 multiqc | 
 Data : ../../../../../scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_trimmed/trimmed_multiqc_report_data| 
 multiqc | Flat-image plots used. Disable with '--interactive'. See docs.| 
 multiqc | MultiQC complete

### Questions/Notes:
1. How many reads did we lose?
	1. check the multiqc_fastqc.txt from each multiqc run and compare the number of "Total Sequences" to find how many reads were removed per file.
2. were adapters actually removed? 
3. Does the number of R1 and R2 reads match? 
4. Note the quality metrics- have they improved? 
5. NOTE ABOUTE PER BASE SEQUENCES: the per base sequence count is still a little wonky for the first few base pairs; sounds like its a normal result of library prep methods which use tagmentation. see [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html) and your kit with tagmentations protocol [here](). 

# zip or delete raw reads
At this point, we will only proceed with the trimmed reads. As such, let's either zip the raw reads to save space or delete them from the working directory.

