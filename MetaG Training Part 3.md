
Part 3: 

----------------------------------------------------------------
## Steps
Step 1: Calculate MAG abundances using coverM

Step 2: 
 
Step 3: 
 
Step 4: 

The Other Stuff...

----------------------------------------------------------------

# Step 1: Calculate MAG abundances using coverM
### Build a Bowtie2 database for mapping

```
#!/bin/bash
#SBATCH --job-name=bowtie2_build_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --output=slurm_output/bowtie2_build_db_%j.out
#SBATCH --error=slurm_output/bowtie2_build_db_%j.err

module purge
module load anaconda
module load bowtie2

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes

# make directory for bowtie database and build it using scaffolds.fna dram file 
mkdir bowtie_DB
cd bowtie_DB

#build a database of scaffolds from the dram scaffolds file
bowtie2-build ../DRAM_1.5_09092026/scaffolds.fna 30_MAG_DB --threads 16
```
sbatch bowtie_build_db.sh
Submitted batch job 32974242

#### map trimmed metagenome reads to bowtie database
```
#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/bowtie2_align_db_%j.out
#SBATCH --error=slurm_output/bowtie2_align_db_%j.err


module purge
module load anaconda
module load bowtie2

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"
BOWTIE_DB="${BASE_DIR}/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/bowtie_DB"

while read -r SAMPLE
do
    R1="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R1_bbduktrimmed.fastq"
    R2="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R2_bbduktrimmed.fastq"

    bowtie2 \
        -D 10 \
        -R 2 \
        -N 1 \
        -L 22 \
        -i S,0,2.50 \
        -p 64 \
        -x "$BOWTIE_DB" \
        -S "${BOWTIE_DB}/${SAMPLE}_mapped_99perMAGs.sam" \
        -1 "$R1" \
        -2 "$R2"

done < "$SAMPLE_LIST"
```

sbatch bowtie_align.sh
Submitted batch job 32974304

## Run coverM for abundances 
#### install coverM
```
### install coverM (which also installs samtools)
# https://github.com/wwood/CoverM?tab=readme-ov-file
conda create -n coverm 
conda activate coverm
conda install --channel bioconda coverm
coverm --version 

module load anaconda
conda activate coverm
```


# all the other stuff.... 
### Run DRAM on just genes

```
#!/bin/bash
#SBATCH --job-name=DRAM_genes
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=20
#SBATCH --time=23:30:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/DRAM_genes_%j.out
#SBATCH --error=slurm_output/DRAM_genes_%j.err 


module load anaconda
conda activate test_dram_again_sept232026

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes

DRAM.py annotate_genes \
-i /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/DRAM_1.5_09092026/genes.faa \
-o ../DRAM_1.5_total_genes \
--threads 20

DRAM.py distill \
-i DRAM_1.5_total_genes/annotations.tsv \
-o DRAM_1.5_total_genes/distill
```
DRAM_genes.sh
Submitted batch job 32973703







## Assembly with IDBA-UD
### Individual assembly with IDBA-UD

```
#!/bin/bash
#SBATCH --job-name=idba_indiv_assembly
#SBATCH --nodes=1  
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=32
#SBATCH --mem=600G
#SBATCH --time=160:00:00
#SBATCH --qos=mem-long
#SBATCH --mail-type=ALL
#SBATCH --partition=amem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/idba_indiv_assembly%j.out
#SBATCH --error=slurm_output/idba_indiv_assembly%j.err

module load miniforge
mamba activate idba

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/files_to_transfer_june12"
BASE_DIR2="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

# number of samples to run at once
MAX_JOBS=3
THREADS=16

while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
(
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && exit
    R1="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R1_bbduktrimmed.fastq"
    R2="${BASE_DIR}/${SAMPLE}/processed_reads/${SAMPLE}_R2_bbduktrimmed.fastq"
    OUTDIR="${BASE_DIR2}/${SAMPLE}/assembly/idba_ud"
    mkdir -p "$OUTDIR"
    FA="${OUTDIR}/${SAMPLE}_interleaved.fa"
    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "[$SAMPLE] Converting FASTQ to FASTA..."
        fq2fa --merge --filter "$R1" "$R2" "$FA"
        echo "[$SAMPLE] Running IDBA-UD..."
        /usr/bin/time -v idba_ud \
            -r "$FA" \
            --pre_correction \
            --num_threads "$THREADS" \
            -o "$OUTDIR" \
            2> "${OUTDIR}/${SAMPLE}_time.log"
        echo "[$SAMPLE] Finished."
        # Optional: remove intermediate FASTA to save space
        rm -f "$FA"
    else
        echo "[$SAMPLE] Missing input reads."
    fi
) &
while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
    wait -n
done
done < "$SAMPLE_LIST"
wait
echo "All samples complete."
```

10b_idba_indiv_assembly.sh
Submitted batch job 29406486

