Part 1: Data upload through trimming reads

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






