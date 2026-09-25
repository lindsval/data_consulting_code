### check amem status here https://xdmod.rc.colorado.edu/
### acompile

```
acompile --ntasks=4 
```
## new alpine slurm partition names

```

#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
```


## loading qiime 2024.10
```

ainteractive --ntasks=4 --time=01:00:00 --partition=acpu --qos=cpu-normal
module purge
module load qiime2/2024.10_amplicon
```

## loading newest qiime install on alpine
```
ainteractive --ntasks=4 --time=03:00:00 --partition=acpu --qos=cpu-normal
module purge
module load qiime2/2026.1_amplicon
```



```
ls -lh /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/Control_BulkSoil_Post_38/assembly/idba_ud/contig.fa
```


```
count=0  
while read sample; do  
compgen -G files_to_transfer_june12/"${sample}/fastqc/trimmed/*R1_bbduktrimmed_fastqc.html" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
#88


count=0  
while read sample; do  
compgen -G files_to_transfer_june12/"${sample}/fastqc/trimmed/*R2_bbduktrimmed_fastqc.html" > /dev/null && ((count++))  
done < sample_list.txt  
  
echo $count
#88

#use multiqc to generate one report for all fastqc_data.txt files
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/files_to_transfer_june12
module load anaconda
conda activate multiqc


multiqc \
/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/files_to_transfer_june12/*/fastqc/raw/*_fastqc.zip \
--outdir /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/multiqc_raw_reads \
--filename raw_multiqc_report.html

```
since i cant get DRAM to be reinstalled in order to share with VT and Chance, i am running VT's for her

copy her files into my scratch
cp 

``
```
#!/bin/bash
#SBATCH --job-name=DRAM_VT
#SBATCH --partition=acpu
#SBATCH --qos=cpu-long
#SBATCH --ntasks=20
#SBATCH --time=160:30:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/DRAM_VT_%j.out
#SBATCH --error=slurm_output/DRAM_VT_%j.err 


module load anaconda
conda activate DRAM_v1.5.0_use

cd /pl/active/courses/2026_summer/CSU_2026/VT_well_water_metaG/dereplicated_genomes

DRAM.py annotate -i '*fa' -o  DRAM_1.5_09092026 --min_contig_size 2500 --threads 20
DRAM.py distill -i DRAM_1.5_09092026/annotations.tsv -o DRAM_1.5_09092026/distill

```
sbatch DRAM_1.5.sh
Submitted batch job 32336021

```
cd /projects/lindsval@colostate.edu
acompile --ntasks=4 
module load anaconda
conda config --set solver libmamba
conda config --show solver
time conda env create -f environment.yaml -n test_install_DRAM_v1.5.0_sept2026
```