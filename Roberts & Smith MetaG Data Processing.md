
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

### zip or delete raw reads
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

### Check all samples assembled
```
#check they were all run
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/

count=0  
while read sample; do  
compgen -G "${sample}/assembly/megahit_out/final.contigs.fa" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#83.....

while read sample; do  
if ! compgen -G "${sample}/assembly/megahit_out/final.contigs.fa" > /dev/null; then  
echo "$sample"  
fi  
done < sample_list.txt

#these samples did not get assembled. 
Control_BulkSoil_Post_6
Control_Rhizo_Post_6
Drought_BulkSoil_Post_7
Drought_Rhizo_Post_7
DroughtDeluge_BulkSoil_Post_8

#why...these were the 5 that were in the middle of assembly when the first assembly job timed out... confirmed in the megahit_25178897.out file. 

#delete the megahit dirs
rm -rf /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG//Control_BulkSoil_Post_6/assembly/megahit_out  
rm -rf /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG//Control_Rhizo_Post_6/assembly/megahit_out  
rm -rf /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG//Drought_BulkSoil_Post_7/assembly/megahit_out  
rm -rf /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG//Drought_Rhizo_Post_7/assembly/megahit_out  
rm -rf /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG//DroughtDeluge_BulkSoil_Post_8/assembly/megahit_out

```


#### Rerun MEGAHIT on the 5 samples that need it

```
#!/bin/bash
#SBATCH --job-name=megahit_rerun_on_5_samples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=55
#SBATCH --qos=normal
#SBATCH --time=23:00:00
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/megahit_rerun_on_5_samples_%j.out
#SBATCH --error=slurm_output/megahit_rerun_on_5_samples_%j.err

module load anaconda
conda activate megahit

BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

megahit \
  -1 "$BASE_DIR/Control_BulkSoil_Post_6/processed_reads/Control_BulkSoil_Post_6_R1_bbduktrimmed.fastq" \
  -2 "$BASE_DIR/Control_BulkSoil_Post_6/processed_reads/Control_BulkSoil_Post_6_R2_bbduktrimmed.fastq" \
  --k-min 31 --k-max 121 --k-step 10 \
  -m 0.4 \
  -t 10 \
  -o "$BASE_DIR/Control_BulkSoil_Post_6/assembly/megahit_out" &

megahit \
  -1 "$BASE_DIR/Control_Rhizo_Post_6/processed_reads/Control_Rhizo_Post_6_R1_bbduktrimmed.fastq" \
  -2 "$BASE_DIR/Control_Rhizo_Post_6/processed_reads/Control_Rhizo_Post_6_R2_bbduktrimmed.fastq" \
  --k-min 31 --k-max 121 --k-step 10 \
  -m 0.4 \
  -t 10 \
  -o "$BASE_DIR/Control_Rhizo_Post_6/assembly/megahit_out" &

megahit \
  -1 "$BASE_DIR/Drought_BulkSoil_Post_7/processed_reads/Drought_BulkSoil_Post_7_R1_bbduktrimmed.fastq" \
  -2 "$BASE_DIR/Drought_BulkSoil_Post_7/processed_reads/Drought_BulkSoil_Post_7_R2_bbduktrimmed.fastq" \
  --k-min 31 --k-max 121 --k-step 10 \
  -m 0.4 \
  -t 10 \
  -o "$BASE_DIR/Drought_BulkSoil_Post_7/assembly/megahit_out" &

megahit \
  -1 "$BASE_DIR/Drought_Rhizo_Post_7/processed_reads/Drought_Rhizo_Post_7_R1_bbduktrimmed.fastq" \
  -2 "$BASE_DIR/Drought_Rhizo_Post_7/processed_reads/Drought_Rhizo_Post_7_R2_bbduktrimmed.fastq" \
  --k-min 31 --k-max 121 --k-step 10 \
  -m 0.4 \
  -t 10 \
  -o "$BASE_DIR/Drought_Rhizo_Post_7/assembly/megahit_out" &

megahit \
  -1 "$BASE_DIR/DroughtDeluge_BulkSoil_Post_8/processed_reads/DroughtDeluge_BulkSoil_Post_8_R1_bbduktrimmed.fastq" \
  -2 "$BASE_DIR/DroughtDeluge_BulkSoil_Post_8/processed_reads/DroughtDeluge_BulkSoil_Post_8_R2_bbduktrimmed.fastq" \
  --k-min 31 --k-max 121 --k-step 10 \
  -m 0.4 \
  -t 10 \
  -o "$BASE_DIR/DroughtDeluge_BulkSoil_Post_8/assembly/megahit_out" &

wait


```
Submitted batch job 25281224

### contigs stats .pl file, save this in a new directory called custom_scripts
```
#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift or die "Usage: $0 <fasta>\n";
open(my $fh, "<", $file) or die "Cannot open $file\n";

my @seqs;
my @headers;
my $seq = "";

# Read FASTA
while (<$fh>) {
    chomp;
    if (/^>/) {
        if ($seq ne "") {
            push @seqs, $seq;
        }
        $seq = "";
        push @headers, $_;  # store full header
    } else {
        $seq .= $_;
    }
}
# push last sequence
if ($seq ne "") {
    push @seqs, $seq;
}

close $fh;

# Lengths
my @lengths = map { length($_) } @seqs;

# Total sequences & bp
my $total_seqs = scalar(@seqs);
my $total_bp = 0;
$total_bp += $_ for @lengths;
my $avg = $total_seqs ? $total_bp / $total_seqs : 0;

# N50 calculation
my @sorted_lengths = sort { $b <=> $a } @lengths;
my $cum = 0;
my $n50 = 0;
for my $len (@sorted_lengths) {
    $cum += $len;
    if ($cum >= $total_bp / 2) {
        $n50 = $len;
        last;
    }
}

# Length bins
my %bins = (
    "0-100"        => [0,100],
    "100-500"      => [100,500],
    "500-1000"     => [500,1000],
    "1000-5000"    => [1000,5000],
    "5000-10000"   => [5000,10000],
    "10000-20000"  => [10000,20000],
    "20000-50000"  => [20000,50000],
    "50000-100000" => [50000,100000],
    "100000-500000"=> [100000,500000],
    "500000+"      => [500000,1e12],
);

# Ordered bins
my @bin_order = (
    "0-100",
    "100-500",
    "500-1000",
    "1000-5000",
    "5000-10000",
    "10000-20000",
    "20000-50000",
    "50000-100000",
    "100000-500000",
    "500000+"
);

# Count sequences in bins
my %counts;
my %bps;
foreach my $len (@lengths) {
    foreach my $bin (keys %bins) {
        my ($min,$max) = @{$bins{$bin}};
        if ($len >= $min && $len < $max) {
            $counts{$bin}++;
            $bps{$bin} += $len;
            last;
        }
    }
}

# Print length distribution
print "Length distribution\n";
print "===================\n\n";
print "Range\t# contigs (%)\t# bps (%)\n";

foreach my $bin (@bin_order) {
    my $c = $counts{$bin} // 0;
    my $b = $bps{$bin} // 0;
    my $c_pct = $total_seqs ? sprintf("%.2f", $c/$total_seqs*100) : 0;
    my $b_pct = $total_bp ? sprintf("%.2f", $b/$total_bp*100) : 0;

    print "$bin:\t$c ($c_pct%)\t$b ($b_pct%)\n";
}

# General info
print "\nGeneral Information\n";
print "==================\n\n";
print "Total number of sequences: $total_seqs\n";
print "Total number of bps:       $total_bp\n";
print "Average sequence length:   " . sprintf("%.2f", $avg) . " bps\n";
print "N50:                       $n50 bps\n";

# Build contig objects
my @contigs;
for (my $i = 0; $i < @seqs; $i++) {
    my $s = $seqs[$i];
    my $len = length($s);

    push @contigs, {
        header => $headers[$i],
        seq    => $s,
        len    => $len,
        gc     => ($s =~ tr/GCgc//),
        nonN   => ($s =~ tr/ACGTacgt//)
    };
}

# Sort by length descending
@contigs = sort { $b->{len} <=> $a->{len} } @contigs;

# Print sequence parameters
print "\nSequence parameters\n";
print "===================\n\n";
print "Sequence\tlength\tG+C\tNon-Ns\tdescription\n";

for (my $i = 0; $i < @contigs; $i++) {

    my $c = $contigs[$i];
    my $header = $c->{header};

    # Extract ONLY contig ID (k121_xxx)
    my ($id) = $header =~ /^>(\S+)/;
    $id = ">$id";

    # Full description (no >)
    my $desc = $header =~ s/^>//r;

    my $len = $c->{len};
    my $gc_pct = $len ? sprintf("%.2f", $c->{gc}/$len*100) : 0;
    my $nonN_pct = $len ? sprintf("%.2f", $c->{nonN}/$len*100) : 0;

    print $i+1, "\t", $id, "\t", $len, "\t", $gc_pct, "\t", $nonN_pct, "\t", $desc, "\n";
}
```

```
chmod +x contig_stats_full.pl
```

### try the new contig stats code - works!

```
#test assembly stats on one metaG

perl /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts/contig_stats_full.pl /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/megahit_out_testing_on_Drought_Rhizo_Post_11/final.contigs.fa > /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/megahit_out_testing_on_Drought_Rhizo_Post_11/Drought_Rhizo_Post_11_final.contigs_STATS_new.txt

```


## run contig_stats on all assemblies

```
#!/bin/bash

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

CONTIG_SCRIPT="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts/contig_stats_full.pl"

# number of samples to run at once
MAX_JOBS=5

while read SAMPLE; do
  (
    CONTIGS="${BASE_DIR}/${SAMPLE}/assembly/megahit_out/final.contigs.fa"
    OUTFILE="${BASE_DIR}/${SAMPLE}/assembly/megahit_out/${SAMPLE}_final.contigs_STATS.txt"

    if [[ -f "$CONTIGS" ]]; then
      echo "Running contig stats for $SAMPLE"

      perl "$CONTIG_SCRIPT" "$CONTIGS" > "$OUTFILE"

    else
      echo "Missing contigs file for $SAMPLE" >&2
    fi
  ) &

  # limit number of concurrent jobs
  if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
    wait -n
  fi

done < "$SAMPLE_LIST"

wait
```

08_contig_stats_loop.sh

```
#!/bin/bash
#SBATCH --job-name=contig_stats_all_samples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --qos=normal
#SBATCH --time=04:00:00
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/contig_stats_%j.out
#SBATCH --error=slurm_output/contig_stats_%j.err

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/slurm

bash 08_contig_stats_loop.sh
```
08_contig_stats.sh

Submitted batch job 25359926

### Check this ran for all samples

```
#check 

count=0  
while read sample; do  
compgen -G "${sample}/assembly/megahit_out/*_final.contigs_STATS.txt" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#88

#good!
```


### Combine all contig stats files

```
#!/bin/bash

BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"
SAMPLE_LIST="${BASE_DIR}/sample_list.txt"
OUTFILE="${BASE_DIR}/all_samples_contig_stats_summary.txt"

# write header
echo -e "Sample\
\t0-100_reads\t0-100_reads_pct\t0-100_bps\t0-100_bps_pct\
\t100-500_reads\t100-500_reads_pct\t100-500_bps\t100-500_bps_pct\
\t500-1000_reads\t500-1000_reads_pct\t500-1000_bps\t500-1000_bps_pct\
\t1000-5000_reads\t1000-5000_reads_pct\t1000-5000_bps\t1000-5000_bps_pct\
\t5000-10000_reads\t5000-10000_reads_pct\t5000-10000_bps\t5000-10000_bps_pct\
\t10000-20000_reads\t10000-20000_reads_pct\t10000-20000_bps\t10000-20000_bps_pct\
\t20000-50000_reads\t20000-50000_reads_pct\t20000-50000_bps\t20000-50000_bps_pct\
\t50000-100000_reads\t50000-100000_reads_pct\t50000-100000_bps\t50000-100000_bps_pct\
\t100000-500000_reads\t100000-500000_reads_pct\t100000-500000_bps\t100000-500000_bps_pct\
\t500000+_reads\t500000+_reads_pct\t500000+_bps\t500000+_bps_pct\
\tTotal_sequences\tTotal_bps\tAvg_length\tN50" > "$OUTFILE"


while read SAMPLE; do

  FILE="${BASE_DIR}/${SAMPLE}/assembly/megahit_out/${SAMPLE}_final.contigs_STATS.txt"

  if [[ ! -f "$FILE" ]]; then
    echo "Missing stats for $SAMPLE" >&2
    continue
  fi

  awk -v sample="$SAMPLE" '
  BEGIN { OFS="\t" }

  /Length distribution/ {in_dist=1; next}
  /General Information/ {in_dist=0; in_gen=1; next}

  # parse distribution lines (robust to whitespace + formatting)
  in_dist && /^[[:space:]]*[0-9]/ {

    # extract range (e.g., 0-100, 100-500, etc.)
    if (match($0, /([0-9]+-[0-9]+|\+):/, m)) {
      range=m[1]
      gsub(":", "", range)
    } else {
      next
    }

    # first match = reads + %
    if (match($0, /([0-9]+)[[:space:]]+\(([0-9.]+)%\)/, r)) {
      reads=r[1]
      reads_pct=r[2]
    } else {
      reads="NA"; reads_pct="NA"
    }

    # second match = bps + %
    rest=substr($0, RSTART + RLENGTH)
    if (match(rest, /([0-9]+)[[:space:]]+\(([0-9.]+)%\)/, b)) {
      bps=b[1]
      bps_pct=b[2]
    } else {
      bps="NA"; bps_pct="NA"
    }

    data[range]=reads"\t"reads_pct"\t"bps"\t"bps_pct
  }

  # general info
  in_gen && /Total number of sequences/ {
    total_seq=$5
  }
  in_gen && /Total number of bps/ {
    total_bps=$5
  }
  in_gen && /Average sequence length/ {
    avg_len=$4
  }
  in_gen && /^N50/ {
    n50=$2
  }

  END {
    printf sample

    ordered_ranges[1]="0-100"
    ordered_ranges[2]="100-500"
    ordered_ranges[3]="500-1000"
    ordered_ranges[4]="1000-5000"
    ordered_ranges[5]="5000-10000"
    ordered_ranges[6]="10000-20000"
    ordered_ranges[7]="20000-50000"
    ordered_ranges[8]="50000-100000"
    ordered_ranges[9]="100000-500000"
    ordered_ranges[10]="500000+"

    for (i=1; i<=10; i++) {
      range_key = ordered_ranges[i]
      if (range_key in data) {
        printf "\t%s", data[range_key]
      } else {
        printf "\tNA\tNA\tNA\tNA"
      }
    }

    printf "\t%s\t%s\t%s\t%s\n", total_seq, total_bps, avg_len, n50
  }

  ' "$FILE" >> "$OUTFILE"

done < "$SAMPLE_LIST"

echo "Done! Output written to: $OUTFILE"
```
bash 08a_combine_stats.sh

## Run Co-Assmebly

#### We will coassemble by treatment and soil type (rhizo versus bulik): 
- DroughtRhizo (n=10; so this includes the pre and post plots)
- DroughtBulk (n=10)
- DelugeRhizo (n=10)
- DelugeBulk (n=10)
- ControlRhizo (n=10)
- ControlBulk (n=10)
- DroughtDelugeRhizo (n=10)
- DroughtDeligeBulk (n=10)
- Control (n=8)

```
#make a new directory for coassembly
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG

mkdir coassembly
cd coassembly

#make subdirectories for each coassembly
mkdir DroughtRhizo
mkdir DroughtBulk
mkdir DelugeRhizo
mkdir DelugeBulk
mkdir ControlRhizo
mkdir ControlBulk
mkdir DroughtDelugeRhizo
mkdir DroughtDelugeBulk
mkdir Control
```

#### Create a sample list the the coassembly 
```
nano coA_sample_list.txt

#paste in the list
DroughtRhizo
DroughtBulk
DelugeRhizo
DelugeBulk
ControlRhizo
ControlBulk
DroughtDelugeRhizo
DroughtDelugeBulk
Control
```
### Create subdirectories
- concat_reads
- assembly
```
for d in */; do  
mkdir -p "${d}concat_reads"  
done

for d in */; do  
mkdir -p "${d}assembly"  
done

```

## Combine the bbduk trimmed reads into r1 and r2 for each treatment
##### ControlBulk
```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/ControlBulk/concat_reads"
samples=(
Control_BulkSoil_Post_10
Control_BulkSoil_Post_24
Control_BulkSoil_Post_27
Control_BulkSoil_Post_38
Control_BulkSoil_Post_6
Control_BulkSoil_Pre_10
Control_BulkSoil_Pre_24
Control_BulkSoil_Pre_27
Control_BulkSoil_Pre_38
Control_BulkSoil_Pre_6
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/ControlBulk_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/ControlBulk_R2_test.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09a_concat_reads_for_coA_ControlBulk.sh

##### ControlRhizo
```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/ControlRhizo/concat_reads"
samples=(
Control_Rhizo_Post_10
Control_Rhizo_Post_24
Control_Rhizo_Post_27
Control_Rhizo_Post_38
Control_Rhizo_Post_6
Control_Rhizo_Pre_10
Control_Rhizo_Pre_24
Control_Rhizo_Pre_27
Control_Rhizo_Pre_38
Control_Rhizo_Pre_6
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/ControlRhizo_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/ControlRhizo_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09b_concat_reads_for_CoA_ControlRhizo.sh


##### Controls
```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/Control/concat_reads"
samples=(
Control1_Control_Pre_NA
Control2_Control_Pre_NA
Control3_Control_Pre_NA
Control4_Control_Pre_NA
Control5_Control_Post_NA
Control6_Control_Post_NA
Control7_Control_Post_NA
Control8_Control_Post_NA
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/Control_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/Control_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09c_concat_reads_for_CoA_Control.sh

##### DelugeBulk

```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/DelugeBulk/concat_reads"
samples=(
Deluge_BulkSoil_Post_23
Deluge_BulkSoil_Post_28
Deluge_BulkSoil_Post_37
Deluge_BulkSoil_Post_5
Deluge_BulkSoil_Post_9
Deluge_BulkSoil_Pre_23
Deluge_BulkSoil_Pre_28
Deluge_BulkSoil_Pre_37
Deluge_BulkSoil_Pre_5
Deluge_BulkSoil_Pre_9
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/DelugeBulk_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/DelugeBulk_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09d_concat_reads_for_CoA_DelugeBulk.sh

##### DelugeRhizo

```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/DelugeRhizo/concat_reads"
samples=(
Deluge_Rhizo_Post_23
Deluge_Rhizo_Post_28
Deluge_Rhizo_Post_37
Deluge_Rhizo_Post_5
Deluge_Rhizo_Post_9
Deluge_Rhizo_Pre_23
Deluge_Rhizo_Pre_28
Deluge_Rhizo_Pre_37
Deluge_Rhizo_Pre_5
Deluge_Rhizo_Pre_9
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/DelugeRhizo_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/DelugeRhizo_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09e_concat_reads_for_CoA_DelugeRhizo.sh

##### DroughtBulk
```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/DroughtBulk/concat_reads"
samples=(
Drought_BulkSoil_Post_11
Drought_BulkSoil_Post_21
Drought_BulkSoil_Post_25
Drought_BulkSoil_Post_40
Drought_BulkSoil_Post_7
Drought_BulkSoil_Pre_11
Drought_BulkSoil_Pre_21
Drought_BulkSoil_Pre_25
Drought_BulkSoil_Pre_40
Drought_BulkSoil_Pre_7
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtBulk_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtBulk_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09f_concat_reads_for_CoA_DroughtBulk.sh

##### DroughtRhizo

```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/DroughtRhizo/concat_reads"
samples=(
Drought_Rhizo_Post_11
Drought_Rhizo_Post_21
Drought_Rhizo_Post_25
Drought_Rhizo_Post_40
Drought_Rhizo_Post_7
Drought_Rhizo_Pre_11
Drought_Rhizo_Pre_21
Drought_Rhizo_Pre_25
Drought_Rhizo_Pre_40
Drought_Rhizo_Pre_7
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtRhizo_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtRhizo_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09g_concat_reads_for_CoA_DroughtRhizo.sh

##### DroughtDelugeBulk
```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/DroughtDelugeBulk/concat_reads"
samples=(
DroughtDeluge_BulkSoil_Post_12
DroughtDeluge_BulkSoil_Post_22
DroughtDeluge_BulkSoil_Post_26
DroughtDeluge_BulkSoil_Post_39
DroughtDeluge_BulkSoil_Post_8
DroughtDeluge_BulkSoil_Pre_12
DroughtDeluge_BulkSoil_Pre_22
DroughtDeluge_BulkSoil_Pre_26
DroughtDeluge_BulkSoil_Pre_39
DroughtDeluge_BulkSoil_Pre_8
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtDelugeBulk_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtDelugeBulk_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```
bash 09h_concat_reads_for_CoA_DroughtDelugeBulk.sh

##### DroughtDelugeRhizo

```
#!/bin/bash

# output directory
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/DroughtDelugeRhizo/concat_reads"
samples=(
DroughtDeluge_Rhizo_Post_12
DroughtDeluge_Rhizo_Post_22
DroughtDeluge_Rhizo_Post_26
DroughtDeluge_Rhizo_Post_39
DroughtDeluge_Rhizo_Post_8
DroughtDeluge_Rhizo_Pre_12
DroughtDeluge_Rhizo_Pre_22
DroughtDeluge_Rhizo_Pre_26
DroughtDeluge_Rhizo_Pre_39
DroughtDeluge_Rhizo_Pre_8
)

BASEDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# concatenate all R1 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R1_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtDelugeRhizo_R1.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done

# concatenate all R2 reads
for sample in "${samples[@]}"; do
    infile="$BASEDIR/$sample/processed_reads/${sample}_R2_bbduktrimmed.fastq"
    outfile="$OUTDIR/DroughtDelugeRhizo_R2.fastq"

    if [[ -f "$infile" ]]; then
        cat "$infile" >> "$outfile"
    else
        echo "Missing: $infile"
    fi
done
```

bash 09i_concat_reads_for_CoA_DroughtDelugeRhizo.sh

### Check file sizes after concatenation - looks good
```
#check file sizes and make sure they make sense (should be around like 100GB each)

while IFS= read -r sample; do  
echo "=== $sample ==="  
  
dir="$sample/concat_reads"  
  
if [[ -d "$dir" ]]; then  
for f in "$dir"/*; do  
[[ -f "$f" ]] || continue  
printf "%s %s\n" "$(basename "$f")" "$(du -h "$f" | cut -f1)"  
done  
else  
echo "Missing directory: $dir"  
fi  
  
echo  
done < coA_sample_list.txt

##OUTPUT
=== DroughtRhizo ===
DroughtRhizo_R1.fastq 91G
DroughtRhizo_R2.fastq 90G
=== DroughtBulk ===
DroughtBulk_R1.fastq 86G
DroughtBulk_R2.fastq 85G
=== DelugeRhizo ===
DelugeRhizo_R1.fastq 90G
DelugeRhizo_R2.fastq 88G
=== DelugeBulk ===
DelugeBulk_R1.fastq 93G
DelugeBulk_R2.fastq 92G
=== ControlRhizo ===
ControlRhizo_R1.fastq 89G
ControlRhizo_R2.fastq 88G
=== ControlBulk ===
ControlBulk_R1.fastq 102G
ControlBulk_R2_test.fastq 101G
=== DroughtDelugeRhizo ===
DroughtDelugeRhizo_R1.fastq 100G
DroughtDelugeRhizo_R2.fastq 98G
=== DroughtDelugeBulk ===
DroughtDelugeBulk_R1.fastq 89G
DroughtDelugeBulk_R2.fastq 88G
=== Control ===
Control_R1.fastq 111M
Control_R2.fastq 110M
```

### Run CoAssembly

```
#!/bin/bash

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"

# number of samples to run at once
MAX_JOBS=6

while read SAMPLE; do
  (
    R1="${BASE_DIR}/${SAMPLE}/concat_reads/${SAMPLE}_R1.fastq"
    R2="${BASE_DIR}/${SAMPLE}/concat_reads/${SAMPLE}_R2.fastq"
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
10_megahit_coassembly_loop.sh

```

#!/bin/bash
#SBATCH --job-name=megahit_coAssembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=65
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/megahitcoAssembly%j.out
#SBATCH --error=slurm_output/megahit_coAssembly%j.err


module load anaconda
conda activate megahit

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/slurm

bash 10_megahit_coassembly_loop.sh 
```
10_megahit_coassembly.sh
i submitted this tuesday afternoon, and its on PD, so i think it should start once the amem partitions are back online.
all good here.
## Coassembly stats

```
check that all coassemblies ran

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly

count=0  
while read sample; do  
compgen -G "${sample}/assembly/megahit_out/final.contigs.fa" > /dev/null &&((count++))  
done < coA_sample_list.txt  
  
echo $count #9, so were good!
```

### Run stats

```
#!/bin/bash

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"

CONTIG_SCRIPT="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts/contig_stats_full.pl"

# number of samples to run at once
MAX_JOBS=5

while read SAMPLE; do
  (
    CONTIGS="${BASE_DIR}/${SAMPLE}/assembly/megahit_out/final.contigs.fa"
    OUTFILE="${BASE_DIR}/${SAMPLE}/assembly/megahit_out/${SAMPLE}_final.contigs_STATS.txt"

    if [[ -f "$CONTIGS" ]]; then
      echo "Running contig stats for $SAMPLE"

      perl "$CONTIG_SCRIPT" "$CONTIGS" > "$OUTFILE"

    else
      echo "Missing contigs file for $SAMPLE" >&2
    fi
  ) &

  # limit number of concurrent jobs
  if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
    wait -n
  fi

done < "$SAMPLE_LIST"

wait
```

11_contig_stats_CoAssembly_loop.sh

```
#!/bin/bash
#SBATCH --job-name=contig_stats_coAssembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --qos=normal
#SBATCH --time=04:00:00
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/contig_stats_coAssembly%j.out
#SBATCH --error=slurm_output/contig_stats_coAssembly%j.err

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/slurm

bash 11_contig_stats_CoAssembly_loop.sh
```
11_contig_stats_coAssembly.sh

Submitted batch job 27468806 (may 22 1015a)

## Compile pullseqs

```
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts

git clone https://github.com/bcthomas/pullseq.git
cd pullseq
mkdir build
cd build

module load cmake #version cmake version 4.2.3
cmake ..  
make

# This will build binaries in build/src/
  > build/src/pullseq
  > build/src/seqdiff

#check its there using help page
./src/pullseq -h

#pullseq - a bioinformatics tool for manipulating fasta and fastq filesVersion: #1.0.2 Name lookup method: UTHASH(Written by bct - copyright 2012-2015)
# ....


#good!
```

# extract contigs >2.5kb using pullseqs

### pullseqs for individual assembly 

```
#!/bin/bash
#SBATCH --job-name=pullseq_filter_indiv_assembly
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=23:00:00
#SBATCH --mem=50gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/pullseqs%j.out
#SBATCH --error=slurm_output/pullseqs%j.err

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"
PULLSEQ="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts/pullseq/build/src/pullseq"

while read SAMPLE; do

  OUTDIR="${BASE_DIR}/${SAMPLE}/assembly/megahit_out"
  INPUT="${OUTDIR}/final.contigs.fa"
  OUTPUT="${OUTDIR}/${SAMPLE}_final.contigs_2500.fa"

  echo "Processing sample: $SAMPLE"

  if [ -f "$INPUT" ]; then
    "$PULLSEQ" -i "$INPUT" -m 2500 > "$OUTPUT"
    echo " Output written to: $OUTPUT"
  else
    echo " WARNING: $INPUT not found, skipping"
  fi

done < "$SAMPLE_LIST"

echo "All samples processed."
```
12_pullseqs_2500_individual_assembly.sh

Submitted batch job 27468972

### Pullseqs for coassembly
```
#!/bin/bash
#SBATCH --job-name=pullseq_filter_Coassembly
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=23:00:00
#SBATCH --mem=50gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/pullseqs_CoA%j.out
#SBATCH --error=slurm_output/pullseqs_CoA%j.err

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"
PULLSEQ="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/custom_scripts/pullseq/build/src/pullseq"

while read SAMPLE; do

  OUTDIR="${BASE_DIR}/${SAMPLE}/assembly/megahit_out"
  INPUT="${OUTDIR}/final.contigs.fa"
  OUTPUT="${OUTDIR}/${SAMPLE}_final.contigs_2500.fa"

  echo "Processing sample: $SAMPLE"

  if [ -f "$INPUT" ]; then
    "$PULLSEQ" -i "$INPUT" -m 2500 > "$OUTPUT"
    echo " Output written to: $OUTPUT"
  else
    echo " WARNING: $INPUT not found, skipping"
  fi

done < "$SAMPLE_LIST"

echo "All samples processed."
```
12_pullseqs_2500_CoAssembly.sh

Submitted batch job 27469069

## Install bbmap (version xz)

```
acompile --ntasks=4 
module load anaconda
conda create -n bbmap
conda activate bbmap
conda install bioconda::bbmap

#get version
bbmap.sh -v
#version 39.81
```


## Map trimmed paired-end reads back to ≥2500 bp assembled contigs to generate coverage information for downstream MAG binning and abundance estimation using bbmap

### Make folder for mapped reads

```
while read SAMPLE; do  
mkdir -p /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/${SAMPLE}/mapped_reads  
done < /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt

while read SAMPLE; do  
mkdir -p /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/${SAMPLE}/mapped_reads  
done < /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt
```

## Individual assembly mapping

```
#!/bin/bash
#SBATCH --job-name=bbmap_indiv_Assembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=23:00:00
#SBATCH --mem=50gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/bbmap_indivAssembly%j.out
#SBATCH --error=slurm_output/bbmap_indivAssembly%j.err

module load anaconda
conda activate bbmap

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"  
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"   

while read SAMPLE; do

    OUTDIR="${BASE_DIR}/${SAMPLE}/assembly/megahit_out"
    REF="${OUTDIR}/${SAMPLE}_final.contigs_2500.fa"
    R1="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R1_bbduktrimmed.fastq"
    R2="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R2_bbduktrimmed.fastq"
    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    OUTPUT="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sam"

    echo "Processing sample: $SAMPLE"

    if [[ -f "$REF" && -f "$R1" && -f "$R2" ]]; then

        bbmap.sh \
            -Xmx48G \
            threads=20 \
            overwrite=t \
            ref="$REF" \
            in1="$R1" \
            in2="$R2" \
            out="$OUTPUT"

        echo "Mapping complete for: $SAMPLE"

    else
        echo "WARNING: Missing files for $SAMPLE"
    fi

done < "$SAMPLE_LIST"

echo "All samples processed."

```
sbatch 13_bbmap_coA.sh
Submitted batch job 27469328
## CoAssembly mapping

```
#!/bin/bash
#SBATCH --job-name=bbmap_CoAssembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=23:00:00
#SBATCH --mem=50gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/bbmap_CoAssembly%j.out
#SBATCH --error=slurm_output/bbmap_CoAssembly%j.err

module load anaconda
conda activate bbmap


SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"  
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"   
while read SAMPLE; do

    OUTDIR="${BASE_DIR}/${SAMPLE}/assembly/megahit_out"
    REF="${OUTDIR}/${SAMPLE}_final.contigs_2500.fa"
    R1="${BASE_DIR}/${SAMPLE}/concat_reads/${SAMPLE}_R1.fastq"
    R2="${BASE_DIR}/${SAMPLE}/concat_reads/${SAMPLE}_R2.fastq"
    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    OUTPUT="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sam"

    echo "Processing sample: $SAMPLE"

    if [[ -f "$REF" && -f "$R1" && -f "$R2" ]]; then

        bbmap.sh \
            -Xmx48G \
            threads=20 \
            overwrite=t \
            ref="$REF" \
            in1="$R1" \
            in2="$R2" \
            out="$OUTPUT"

        echo "Mapping complete for: $SAMPLE"

    else
        echo "WARNING: Missing files for $SAMPLE"
    fi

done < "$SAMPLE_LIST"

echo "All samples processed."

```
sbatch 13_bbmap_coA.sh
Submitted batch job 27469284


### Convert SAM to BAM files, sort, filter

