
### received pre-demuxed samples, 95 samples
used the demux.qza to start with- see if primers are in sequence since there are 515f/806r on a 2x300 cycle kit mi seq i100


### Files: 
- /Volumes/RSTOR-Sones_Lab/16S/qiime2_files/01_Mouse_fecal_16S_Jun_2026_demux.qza
- metadata= '/Volumes/RSTOR-Sones_Lab/16S/Metadata Mouse Fecal 16S June 2026.xlsx'
- moved these files to alpine and to my computer: 
	- My computer: /Users/valerielindstrom/Documents/PostDoc/data_consulting/sones_lab_16s
	- Alpine: /scratch/alpine/lindsval@colostate.edu/sones_16S



## combine the lanes of the demultiplexed data

```
#raw data
cd /Volumes/RSTOR-Sones_Lab/16S/Microbiome Workshop 2026 (2)/BCLConvert_07_06_2026_14_41_39Z-938987049

#create sample list file to combine the lanes into 1 r1 and 1 r2 file.
while read SAMPLE; do
    L1=$(find . -maxdepth 1 -type d -name "${SAMPLE}_L1-ds.*")
    L2=$(find . -maxdepth 1 -type d -name "${SAMPLE}_L2-ds.*")

    if [[ -z "$L1" || -z "$L2" ]]; then
        echo "MISSING: $SAMPLE"
    fi
done < sample_list.txt

mkdir combined_fastq

while read SAMPLE; do

    L1=$(find . -maxdepth 1 -type d -name "${SAMPLE}_L1-ds.*" | head -1)
    L2=$(find . -maxdepth 1 -type d -name "${SAMPLE}_L2-ds.*" | head -1)

    echo "Processing ${SAMPLE}..."

    cat "${L1}"/*_L001_R1_001.fastq.gz \
        "${L2}"/*_L002_R1_001.fastq.gz \
        > "combined_fastq/${SAMPLE}_R1.fastq.gz"

    cat "${L1}"/*_L001_R2_001.fastq.gz \
        "${L2}"/*_L002_R2_001.fastq.gz \
        > "combined_fastq/${SAMPLE}_R2.fastq.gz"

done < sample_list.txt

#check they all concatenated (should get 97)
ls combined_fastq/*.fastq.gz | wc -l

```


## move these files to alpine - done


## Import the reads and plot read quality

```
mkdir -p /scratch/alpine/lindsval@colostate.edu/sones_16S/raw_reads/demux_reads/NP_F_only_samples
cat > /scratch/alpine/lindsval@colostate.edu/sones_16S/raw_reads/demux_reads/NP_F_only_samples/sample_list.txt << 'EOF'
Sones_BPH5_TK_32
Sones_BPH5_TK_33
Sones_BPH5_TK_35
Sones_BPH5_TK_41
Sones_BPH5_TK_43
Sones_BPH5_TK_44
Sones_BPH5_LD_26
Sones_BPH5_LD_27
Sones_BPH5_LD_29
Sones_BPH5_LD_30
Sones_BPH5_LD_31
Sones_BPH5_LD_32
Sones_C57_TK_8
Sones_C57_TK_10
Sones_C57_TK_11
Sones_C57_TK_12
Sones_C57_TK_13
Sones_C57_TK_34
Sones_C57_LD_33
Sones_C57_LD_30
Sones_C57_LD_36
Sones_C57_LD_37
Sones_C57_LD_39
Sones_C57_LD_42
EOF

wc -l /scratch/alpine/lindsval@colostate.edu/sones_16S/raw_reads/demux_reads/NP_F_only_samples/sample_list.txt

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/raw_reads/demux_reads/combined_fastq

while read SAMPLE; do
    cp "${SAMPLE}_R1.fastq.gz" ../NP_F_only_samples/
    cp "${SAMPLE}_R2.fastq.gz" ../NP_F_only_samples/
done < ../NP_F_only_samples/sample_list.txt

#filenames
ls -1 *.fastq.gz > filename_list.txt

```

```
#create manifest
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv

for r1 in "$PWD"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${PWD}/${sample}_R2.fastq.gz"
    echo -e "${sample}\t${r1}\t${r2}"
done >> manifest.tsv
```

## make dirs

```
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/
mkdir demux
mkdir taxonomy
mkdir metadata
mkdir tree
mkdir taxaplots
mkdir dada2
```


### import reads

```
ainteractive --ntasks=4 --time=03:00:00 --partition=acpu --qos=cpu-normal
module purge
module load qiime2/2026.1_amplicon
```

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /scratch/alpine/lindsval@colostate.edu/sones_16S/raw_reads/demux_reads/NP_F_only_samples/manifest.tsv \
  --output-path demux/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```


```
qiime demux summarize \
--i-data demux/demux.qza \
--o-visualization demux/demux.qzv
```

# original analysis before i realized the Sones_BPH5_LD_26 sample was missing.
## Cutadapt to remove primers: 

```
#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=amilan
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

#Activate qiime
module purge
module load qiime2/2024.10_amplicon

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/demux

qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-adapter-f ATTAGAWACCCVNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
--p-match-adapter-wildcards \
--p-match-read-wildcards \
--o-trimmed-sequences filtered_reads_cutadapt.qza \
--p-discard-untrimmed \
--verbose
```
Submitted batch job 30685241

## Denoise with dada2 

```
#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=amilan
#SBATCH --time=05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

#Activate qiime
module purge
module load qiime2/2024.10_amplicon

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/dada2

# dada2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../demux/filtered_reads_cutadapt.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--o-table table_dada2.qza \
--o-representative-sequences rep_seqs_dada2.qza \
--o-denoising-stats denoising_stats_dada2.qza

# visualize outputs
qiime feature-table summarize \
  --i-table table_dada2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization table_dada2.qzv 

qiime feature-table tabulate-seqs \
--i-data rep_seqs_dada2.qza \
--o-visualization rep_seqs_dada2.qzv

qiime metadata tabulate \
--m-input-file denoising_stats_dada2.qza \
--o-visualization denoising_stats_dada2.qzv
```

```
sbatch dada2.sh
```
Submitted batch job 30685230

## Results of denoising:

- number of ASVs =  876
- median reads per sample =  110,951
- Avg reads retained after denoising =  ~85-90%
- no positive controls
- 1 negative control and it only has 32 reads!
- were there long amplicons that need to be removed? No. longest read was 255. 

### taxonomy w/ GG2 2024.10
```
# get the classifier
wget --no-check-certificate https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.v4.nb.qza 

#classify
qiime feature-classifier classify-sklearn \
  --i-reads ../dada2/rep_seqs_dada2.qza \
  --i-classifier 2024.09.backbone.v4.nb.qza \
  --o-classification taxonomy_nb_gg2.qza

# filter tables (also remove the additional mito genome - sp004296775)
qiime taxa filter-table \
  --i-table ../dada2/table_dada2.qza \
  --i-taxonomy taxonomy_nb_gg2.qza \
  --p-exclude mitochondria,chloroplast,sp004296775 \
  --o-filtered-table ../dada2/table_noMitoChloro_nb_GG2.qza
```


## Filter tables
```

#check table to see if any samples were lost due to mito and chloro filtering
qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_GG2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization ../dada2/table_noMitoChloro_nb_GG2.qzv 

# remove all features with a total abundance of less than 10 from GG2 table
qiime feature-table filter-features \
--i-table ../dada2/table_noMitoChloro_nb_GG2.qza \
--p-min-frequency 10 \
--o-filtered-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qza

#check table to see if any samples were lost due to low abundance features
qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qzv 

# remove features that show up in only a single sample
qiime feature-table filter-features \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qza \
--p-min-samples 2 \
--o-filtered-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza

qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qzv 
  
```
What is lost? anything important? 
Very little lost due to mito/chloro filtering (total = 871)
692 features remain after filtering for min samples, min freq

## taxa plot of non-rarefied table
```

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxaplot_noMitoChloro_nb_GG2_minfreq10_minsample2.qzv
```

## Alpha rarefaction
```
qiime diversity alpha-rarefaction \
--i-table dada2/table_noMitoChloro_nb_GG2.qza \
--m-metadata-file metadata/metadata.txt \
--o-visualization alpha_rarefaction_curve.qzv \
--p-min-depth 10 \
--p-max-depth 100000
```


## Generate taxa plots of rarefied tables 
```
cd ../
cd taxaplots

qiime feature-table rarefy \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
--p-sampling-depth 40000 \
--o-rarefied-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qza

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxaplot_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qzv
```


## filter out control samples from tables
```
cd ../dada2

#filter out controls from the unrarefied table (use for core metrics)

qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_GG2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[sampleID] != 'NTC'" \
  --o-filtered-table table_noMitoChloro_nb_GG2_noControls.qza  

qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[sampleID] != 'NTC'" \
  --o-filtered-table table_noMitoChloro_nb_GG2_minfreq10_minsample2_noControls.qza  

#filter out controls from the rarefied table
qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[sampleID] != 'NTC'" \
  --o-filtered-table table_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k_noControls.qza
```

## Generate a phylogenetic tree (SEPP tree, gg2)
```
#!/bin/bash
#SBATCH --job-name=gg2_tree_sones
#SBATCH --nodes=1
#SBATCH --ntasks=23
#SBATCH --partition=amilan
#SBATCH --time=22:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

   
#### SEPP tree w/ gg2
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/tree

wget --no-check-certificate https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.sepp-reference.qza 

qiime fragment-insertion sepp \
--i-representative-sequences ../dada2/rep_seqs_dada2.qza \
--i-reference-database 2022.10.backbone.sepp-reference.qza \
--o-tree tree_gg2.qza \
--o-placements tree_placements_gg2.qza \
--p-threads 4
```

```
sbatch tree.sh
```
Submitted batch job 30816208


## Core Metrics (all samples)
```
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/

#core metrics 
qiime diversity core-metrics-phylogenetic \
--i-phylogeny tree/tree_gg2.qza \
--i-table dada2/table_noMitoChloro_nb_GG2.qza \
--p-sampling-depth 40000 \
--m-metadata-file metadata/metadata.txt \
--output-dir core_metrics_rare40k_gg2
```



## After meeting with the Sones Lab, all these sampels are not necessary for the manuscript revision. so need to fileter and redo core metrics.

- C57 = controls, BH5 (hypertension) strains of mice
- Right now, they want to know about **non-pregnant females between diet**. NP is not pregnant.
- Diet = lab chow (control = TD)
- Treatment= Hypertension is the TK it’s a plants(soy) based diet. Soybean oil replaces animal fat
	- Possibly soy pased diet had neg effects on adiposity
- mostly interested in SCFA producing bacteria
- They have SCFA data - do they have it for the NP females??

- f__Helicobacteraceae promote hypotentions
- Possibly acetate is lacking in the Hypertension group, gram neg obligate aerobes.




## Filter tables to keep only the non-pregnant females
```
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/dada2
qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_GG2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[Day_of_pregnancy] = 'NP' AND [Sex] = 'F'" \
  --o-filtered-table table_noMitoChloro_nb_GG2_NP_F_only.qza 

qiime feature-table summarize \
  --i-table table_noMitoChloro_nb_GG2_NP_F_only.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization table_noMitoChloro_nb_GG2_NP_F_only.qzv 
  
# remove all features with a total abundance of less than 10 from GG2 table
qiime feature-table filter-features \
--i-table table_noMitoChloro_nb_GG2_NP_F_only.qza \
--p-min-frequency 10 \
--o-filtered-table table_noMitoChloro_nb_GG2_NP_F_only_minfreq10.qza

#check table to see if any samples were lost due to low abundance features
qiime feature-table summarize \
  --i-table table_noMitoChloro_nb_GG2_NP_F_only_minfreq10.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization table_noMitoChloro_nb_GG2_NP_F_only_minfreq10.qzv 

# remove features that show up in only a single sample
qiime feature-table filter-features \
--i-table table_noMitoChloro_nb_GG2_NP_F_only_minfreq10.qza \
--p-min-samples 2 \
--o-filtered-table table_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2.qza

qiime feature-table summarize \
  --i-table table_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization table_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2.qzv 
  
```


there are now 23 samples: 
	for some reason the following sample was not demultiplex (i had used david's demux file..) = Sones_BPH5_LD_26
- 12 control phenotype females, NP
	- 6 were on the treatment diet
	- 6 were on regular lab chow
- 11 treatment phenotype females, NP
	- 6 were on the treatment diet
	- 6 were on regular lab chow


## Generate taxa plots of rarefied tables 
```
cd ../
cd taxaplots

qiime feature-table rarefy \
--i-table ../dada2/table_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2.qza \
--p-sampling-depth 40000 \
--o-rarefied-table ../dada2/table_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2_rare40k.qza

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2_rare40k.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxaplot_noMitoChloro_nb_GG2_NP_F_only_minfreq10_minsample2_rare40k.qzv
```


## Core Metrics (NP/F samples)
```

qiime diversity core-metrics-phylogenetic \
--i-phylogeny tree/tree_gg2.qza \
--i-table dada2/table_noMitoChloro_nb_GG2_NP_F_only.qza \
--p-sampling-depth 40000 \
--m-metadata-file metadata/metadata.txt \
--output-dir core_metrics_rare40k_gg2_NP_F_only




```
## ANCOM-BC2 (level 7)
```
module purge
module load qiime2/2026.1_amplicon

mkdir ancombc2
cd ancombc2

```

```
qiime feature-table filter-samples \
--i-table ../dada2/table_noMitoChloro_nb_GG2_NP_F_only.qza \
--p-min-frequency 40000 \
--o-filtered-table table_noMitoChloro_nb_GG2_NP_F_only_40000.qza
```

```
# set to 3 so tht a feature must be observed in 3 sampeles (and our groups are sizes of 4)
qiime feature-table filter-features \
--i-table table_noMitoChloro_nb_GG2_NP_F_only_40000.qza \
--p-min-frequency 50 \
--p-min-samples 3 \
--o-filtered-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund.qza

qiime taxa collapse \
--i-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--p-level 7 \
--o-collapsed-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund_L7.qza

```

## identify features associated with:
- Diet
- Strain
- Diet × Strain

```
#added a new column in metadata called Group that has all 4 groups designated in 1 column, since ancom wasnt able tod o interaction term

qiime composition ancombc2 \
--i-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund_L7.qza \
--m-metadata-file ../metadata/metadata_for_ancom.txt \
--p-fixed-effects-formula 'Group' \
--p-reference-levels 'Group::C57_LD' \
--o-ancombc2-output ancombc2_full_L7.qza

qiime composition ancombc2-visualizer \
--i-data ancombc2_full_L7.qza \
--o-visualization ancombc2_full_L7.qzv

#Overall diet effect

qiime composition ancombc2 \
  --i-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund_L7.qza \
  --m-metadata-file ../metadata/metadata_for_ancom.txt \
  --p-fixed-effects-formula 'Diet + Strain' \
  --p-reference-levels 'Diet::LD' \
  --o-ancombc2-output ancombc2_Diet_L7.qza

qiime composition ancombc2-visualizer \
--i-data ancombc2_Diet_L7.qza \
--o-visualization ancombc2_Diet_L7.qzv


# Overall strain effect
qiime composition ancombc2 \
  --i-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund_L7.qza \
  --m-metadata-file ../metadata/metadata_for_ancom.txt \
  --p-fixed-effects-formula 'Strain + Diet' \
  --p-reference-levels 'Strain::C57' \
  --o-ancombc2-output ancombc2_Strain_L7.qza

qiime composition ancombc2-visualizer \
--i-data ancombc2_Strain_L7.qza \
--o-visualization ancombc2_Strain_L7.qzv

### export
qiime tools export \
  --input-path ancombc2_full_L7.qza \
  --output-path ancombc2_full_L7_export

ls -lh ancombc2_full_L7_export/



```
## ANCOM-BC2 (level 6)
```

module purge
module load qiime2/2026.1_amplicon


```

```

cd ancombc2
qiime taxa collapse \
--i-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--p-level 6 \
--o-collapsed-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund_L6.qza

qiime composition ancombc2 \
--i-table table_noMitoChloro_nb_GG2_NP_F_only_40000_abund_L6.qza \
--m-metadata-file ../metadata/metadata_for_ancom.txt \
--p-fixed-effects-formula 'Group' \
--p-reference-levels 'Group::C57_LD' \
--o-ancombc2-output ancombc2_full_L6.qza

qiime composition ancombc2-visualizer \
--i-data ancombc2_full_L6.qza \
--o-visualization ancombc2_full_L6.qzv

qiime tools export \
  --input-path ancombc2_full_L6.qza \
  --output-path ancombc2_full_L6


## Export alpha and beta diversity files to then run in R
cd /scratch/alpine/lindsval@colostate.edu/sones_16S

mkdir export

#shannon
unzip core_metrics_rare40k_gg2_NP_F_only/shannon_vector.qza -d export/shannon

# Observed Features  
unzip core_metrics_rare40k_gg2_NP_F_only/observed_features_vector.qza -d export/observed_features  
  
# Faith's PD  
unzip core_metrics_rare40k_gg2_NP_F_only/faith_pd_vector.qza -d export/faith_pd  
  
# Pielou's evenness  
unzip core_metrics_rare40k_gg2_NP_F_only/evenness_vector.qza -d export/evenness

# Bray Curtis  
unzip core_metrics_rare40k_gg2_NP_F_only/bray_curtis_pcoa_results.qza -d export/bray_curtis
 
  
# Jaccard  
unzip core_metrics_rare40k_gg2_NP_F_only/jaccard_pcoa_results.qza -d export/jaccard  
  
# Unweighted Unifrac  
unzip core_metrics_rare40k_gg2_NP_F_only/unweighted_unifrac_pcoa_results.qza -d export/unweighted_unifrac  
  
# Weighted Unifrac  
unzip core_metrics_rare40k_gg2_NP_F_only/weighted_unifrac_pcoa_results.qza -d export/weighted_unifrac

# define alpha metrics  
metrics=("shannon" "evenness" "faith_pd" "observed_features")  
  
# copy their tsv files into export/  
for metric in "${metrics[@]}"; do  
 cp $metric/*/data/alpha-diversity.tsv ${metric}.tsv  
done

# define beta metrics  
metrics=("bray_curtis" "jaccard" "unweighted_unifrac" "weighted_unifrac")  
# copy their txt files into export  
for metric in "${metrics[@]}"; do  
 cp $metric/*/data/ordination.txt ${metric}.txt  
done

```
```

```

```

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/core_metrics_rare40k_gg2_NP_F_only/
unzip rarefied_table.qza
mv 8762d88b-3f2d-4fd0-9dda-7da29a7b6ce0/ rarefied_table
cd rarefied_table
biom convert -i feature-table.biom -o feature_table_rare40k.tsv --to-tsv

unzip bray_curtis_distance_matrix.qza
mv c9a154c6-bd8a-456f-8c4f-dbeacbb831fe/ bray_curtis_distance_matrix
cd bray_curtis_distance_matrix


```

#### now download these files and put them here

/Users/valerielindstrom/Documents/PostDoc/data_consulting/sones_lab_16s/export

#### then use R script for plotting and stats
/Users/valerielindstrom/Documents/PostDoc/data_consulting/sones_lab_16s/figures/alpha_beta_taxa_plots.R

# This new section is for the re-do analysis with all 24 samples
### note that i didn't bring over the NTC (already checked it in the previous analysis)
### basically takes all the code from above and re-runs it using same params in one script

```
#!/bin/bash
#SBATCH --job-name=submit_all_commands_sones16S
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=acpu
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=cpu-normal

#Activate qiime
module purge
module load qiime2/2026.1_amplicon

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/demux

qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-adapter-f ATTAGAWACCCVNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \
--p-match-adapter-wildcards \
--p-match-read-wildcards \
--o-trimmed-sequences filtered_reads_cutadapt.qza \
--p-discard-untrimmed \
--verbose

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/dada2

# dada2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../demux/filtered_reads_cutadapt.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--o-table table_dada2.qza \
--o-representative-sequences rep_seqs_dada2.qza \
--o-denoising-stats denoising_stats_dada2.qza \
--o-base-transition-stats base_stats_dada2.qza

# visualize outputs
qiime feature-table summarize \
  --i-table table_dada2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --o-summary table_dada2.qzv \
  --o-feature-frequencies feature_frequencies.qza \
  --o-sample-frequencies sample_frequencies.qza

qiime feature-table tabulate-seqs \
--i-data rep_seqs_dada2.qza \
--o-visualization rep_seqs_dada2.qzv

qiime metadata tabulate \
--m-input-file denoising_stats_dada2.qza \
--o-visualization denoising_stats_dada2.qzv

cd /scratch/alpine/lindsval@colostate.edu/sones_16S/taxonomy
# get the classifier
wget --no-check-certificate https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.v4.nb.qza 

#classify
qiime feature-classifier classify-sklearn \
  --i-reads ../dada2/rep_seqs_dada2.qza \
  --i-classifier 2024.09.backbone.v4.nb.qza \
  --o-classification taxonomy_nb_gg2.qza

# filter tables (also remove the additional mito genome - sp004296775)
qiime taxa filter-table \
  --i-table ../dada2/table_dada2.qza \
  --i-taxonomy taxonomy_nb_gg2.qza \
  --p-exclude mitochondria,chloroplast,sp004296775 \
  --o-filtered-table ../dada2/table_noMitoChloro_nb_GG2.qza

#check table to see if any samples were lost due to mito and chloro filtering
qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_GG2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --o-summary ../dada2/table_noMitoChloro_nb_GG2.qzv \
  --o-feature-frequencies ../dada2/feature_frequencies_table_noMitoChloro_nb_GG2.qza
  --o-sample-frequencies ../dada2/sample_frequenciestable_noMitoChloro_nb_GG2.qza

# remove all features with a total abundance of less than 10 from GG2 table
qiime feature-table filter-features \
--i-table ../dada2/table_noMitoChloro_nb_GG2.qza \
--p-min-frequency 10 \
--o-filtered-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qza

#check table to see if any samples were lost due to low abundance features
qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --o-summary ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qzv \
  --o-feature-frequencies ../dada2/feature_frequencies_table_noMitoChloro_nb_GG2minfreq10.qza
  --o-sample-frequencies ../dada2/sample_frequenciestable_noMitoChloro_nb_GG2minfreq10.qza

# remove features that show up in only a single sample
qiime feature-table filter-features \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10.qza \
--p-min-samples 2 \
--o-filtered-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza

qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --o-summary ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qzv \
  --o-feature-frequencies ../dada2/feature_frequencies_table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza
  --o-sample-frequencies ../dada2/sample_frequenciestable_noMitoChloro_nb_GG2_minfreq10_minsample2.qza
 
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/taxaplots

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxaplot_noMitoChloro_nb_GG2_minfreq10_minsample2.qzv

cd ../
qiime diversity alpha-rarefaction \
--i-table dada2/table_noMitoChloro_nb_GG2.qza \
--m-metadata-file metadata/metadata.txt \
--o-visualization alpha_rarefaction_curve.qzv \
--p-min-depth 10 \
--p-max-depth 100000
  
cd taxaplots
qiime feature-table rarefy \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
--p-sampling-depth 40000 \
--o-rarefied-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qza

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxaplot_noMitoChloro_nb_GG2_minfreq10_minsample2_rare40k.qzv


#### SEPP tree w/ gg2
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/tree

wget --no-check-certificate https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.sepp-reference.qza 

qiime fragment-insertion sepp \
--i-representative-sequences ../dada2/rep_seqs_dada2.qza \
--i-reference-database 2022.10.backbone.sepp-reference.qza \
--o-tree tree_gg2.qza \
--o-placements tree_placements_gg2.qza \
--p-threads 4

#core metrics 
qiime diversity core-metrics-phylogenetic \
--i-phylogeny tree/tree_gg2.qza \
--i-table dada2/table_noMitoChloro_nb_GG2.qza \
--p-sampling-depth 40000 \
--m-metadata-file metadata/metadata.txt \
--output-dir core_metrics_rare40k_gg2

```
cd /scratch/alpine/lindsval@colostate.edu/sones_16S/slurm
sbatch `submit_all.sh`
Submitted batch job 31768393


~={red}### then redo the ancombc stuff.=~
