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

```

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/

count=0  
while read sample; do  
compgen -G "${sample}/processed_reads/*R1_bbduktrimmed.fastq" > /dev/null &&((count++))  
done < sample_list.txt  
  
echo $count
# only31 finished, ran out of storage on alpine..... 

count=0

while read -r sample; do
    for file in "${sample}"/processed_reads/*_bbduktrimmed.fastq; do
        if [[ -f "$file" ]]; then
            echo "$file"
            ((count++))
        fi
    done
done < sample_list.txt

echo "Total files: $count"