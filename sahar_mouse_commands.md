uses qiime v2024.10
analysis by val seitz
feb 2026 
Note: these data were generated on the new MiSeq i100 instrument with 515f/926r primers
- this means that we are targeting the V4–V5 region of the 16S rRNA gene. 

```
ainteractive --ntasks=4 --time=01:00:00
module purge
module load qiime2/2024.10_amplicon
```

## Import and Demultiplex reads 
```
#!/bin/bash
#SBATCH --job-name=import_demux_sahar_mouse
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=amilan
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm-%j.out
#SBATCH --qos=normal

module purge
module load qiime2/2024.10_amplicon

cd /scratch/alpine/lindsval@colostate.edu/sahar_data_2026

#import forward and reverse reads 
qiime tools import \
--type EMPPairedEndSequences \
--input-path raw_reads \
--output-path mouse_reads.qza

qiime demux emp-paired \
--m-barcodes-file metadata/mouse_barcodes.txt \
--m-barcodes-column barcode \
--p-rev-comp-mapping-barcodes \
--p-rev-comp-barcodes \
--i-seqs mouse_reads.qza \
--o-per-sample-sequences demux/demux_mouse.qza \
--o-error-correction-details demux/demux_demux_mouse_details.qza

qiime demux summarize \
--i-data demux/demux_mouse.qza \
--o-visualization demux/demux_mouse.qzv
```

```
sbatch import_demux_sahar_mouse.sh
```
Submitted batch job 24302503

Data Quality: Reads looks wonderful! these are MiSeq i100 data, so quality scores are binned. tuncated both f and r at 250:
## Denoise 

```
#!/bin/bash
#SBATCH --job-name=dada2_sahar_mouse
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

cd /scratch/alpine/lindsval@colostate.edu/sahar_data_2026/dada2

# dada2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../demux/demux_mouse.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 300 \
--p-trunc-len-r 300 \
--o-table mouse_table_dada2.qza \
--o-representative-sequences mouse_rep_seqs_dada2.qza \
--o-denoising-stats mouse_denoising_stats_dada2.qza

# visualize outputs
qiime feature-table summarize \
  --i-table mouse_table_dada2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization mouse_table_dada2.qzv 

qiime feature-table tabulate-seqs \
--i-data mouse_rep_seqs_dada2.qza \
--o-visualization mouse_rep_seqs_dada2.qzv

qiime metadata tabulate \
--m-input-file mouse_denoising_stats_dada2.qza \
--o-visualization mouse_denoising_stats_dada2.qzv
```

```
sbatch dada2.sh
```
Submitted batch job 24447466

## Results of denoising:

- number of ASVs =  1551
- mean reads per sample =  63k
- Avg reads retained after denoising =  ~85%
- positive controls have lots of reads - good!
- negative controls have really low reads <200, good!
- were there long amplicons that need to be removed? Yes, longest is 558, anything greataer than 400 needs to be removed

### Remove long (400+ base pair) amplicons from the representative sequences file and the feature table

```
# filter out any large amplicons from the seqs and table (because they are contaminates)

cd dada2

qiime feature-table filter-seqs \
--i-data mouse_rep_seqs_dada2.qza \
--m-metadata-file mouse_rep_seqs_dada2.qza \
--p-where 'length(sequence) < 400' \
--o-filtered-data mouse_rep_seqs_dada2_filtered.qza

qiime feature-table tabulate-seqs \
--i-data mouse_rep_seqs_dada2_filtered.qza \
--o-visualization mouse_rep_seqs_dada2_filtered.qzv

qiime feature-table filter-features \
--i-table mouse_table_dada2.qza \
--m-metadata-file mouse_rep_seqs_dada2_filtered.qza \
--o-filtered-table mouse_table_dada2_filtered.qza
  
qiime feature-table summarize \
--i-table mouse_table_dada2_filtered.qza \
--m-sample-metadata-file ../metadata/metadata.txt \
--o-visualization mouse_table_dada2_filtered.qzv
    
```

- removed 87 bad features
- only keeps the reads which are primarily 374nt long, great!
- mean frequency per sample is still around 64k so that means not many reads when towards those bad ASVs anyways. 

## Taxonomy w/ GG2 2024.10 using the full-length 2024.09.backbone.full-length.nb.qza 

```
cd taxonomy
module purge
module load qiime2/2024.10_amplicon

# get the classifier
wget --no-check-certificate https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.nb.qza

#classify
qiime feature-classifier classify-sklearn \
  --i-reads ../dada2/mouse_rep_seqs_dada2_filtered.qza \
  --i-classifier 2024.09.backbone.full-length.nb.qza \
  --o-classification taxonomy_nb_gg2.qza
```


## Filter tables: remove mito and chloro
```
# filter tables (also remove the additional mito genome - sp004296775)
qiime taxa filter-table \
  --i-table ../dada2/mouse_table_dada2_filtered.qza \
  --i-taxonomy taxonomy_nb_gg2.qza \
  --p-exclude mitochondria,chloroplast,sp004296775 \
  --o-filtered-table ../dada2/table_noMitoChloro_nb_GG2.qza
  
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
~={red}What is lost? anything important? =~
~={red}- Very little lost due to mito/chloro filtering (total = 1458)=~
~={red}- over half the features were lost from freq 10 filtering (remaining 731)=~
~={red}- 656 features left after keeping only features in 2 samples...=~

~={red}this seems really low...jessica said fine=~

## Taxa plots 
## all samples, filtered but not rarefied

```
cd ../taxaplots

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_GG2_minfreq10_minsample2.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxa_plot_mouse_noMitoChloro_nb_GG2_minfreq10_minsample2.qzv

```
## Alpha rarefaction
```
qiime diversity alpha-rarefaction \
--i-table dada2/table_noMitoChloro_nb_GG2.qza \
--m-metadata-file metadata/metadata.txt \
--o-visualization alpha_rarefaction_curve.qzv \
--p-min-depth 10 \
--p-max-depth 75000
```

~={red}Where to rarefy: 40k=~

# train a nb clasifier on the v4-v5 region from the gg2 full-length backbone

Get the files from greengenes2 server: https://ftp.microbio.me/greengenes_release/current/
get these files: 
- 2024.09.backbone.full-length.fna.qza
- 2024.09.backbone.tax.qza

```
mkdir gg2_v4v5_nb_classifier
cd gg2_v4v5_nb_classifier

#extract v4-v5
qiime feature-classifier extract-reads \
--i-sequences 2024.09.backbone.full-length.fna.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer CCGYCAATTYMTTTRAGTTT \
--p-min-length 300 \
--p-max-length 400 \
--o-reads gg2_v4v5_seqs.qza


#train the nb classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads gg2_v4v5_seqs.qza \
--i-reference-taxonomy 2024.09.backbone.tax.qza \
--o-classifier gg2_2024.09_v4v5_nb_classifier.qza
```


### classify taxonomy using the new classifier

```
cd taxonomy

#classify
qiime feature-classifier classify-sklearn \
--i-reads ../dada2/mouse_rep_seqs_dada2_filtered.qza \
--i-classifier ../gg2_v4v5_nb_classifier/gg2_2024.09_v4v5_nb_classifier.qza \
--o-classification taxonomy_nb_v4v5_gg2.qza
```


## Filter tables: remove mito and chloro from v4v5 tables
```
# filter tables (also remove the additional mito genome - sp004296775)
qiime taxa filter-table \
  --i-table ../dada2/mouse_table_dada2_filtered.qza \
  --i-taxonomy taxonomy_nb_v4v5_gg2.qza \
  --p-exclude mitochondria,chloroplast,sp004296775 \
  --o-filtered-table ../dada2/table_noMitoChloro_nb_v4v5_GG2.qza
  
#check table to see if any samples were lost due to mito and chloro filtering
qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization ../dada2/table_noMitoChloro_nb_v4v5_GG2.qzv 

# remove all features with a total abundance of less than 10 from GG2 table
qiime feature-table filter-features \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2.qza \
--p-min-frequency 10 \
--o-filtered-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10.qza

#check table to see if any samples were lost due to low abundance features
qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10.qzv 

# remove features that show up in only a single sample
qiime feature-table filter-features \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10.qza \
--p-min-samples 2 \
--o-filtered-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2.qza

qiime feature-table summarize \
  --i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2.qza \
  --m-sample-metadata-file ../metadata/metadata.txt \
  --o-visualization ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2.qzv 
  
```
~={red}What is lost? anything important? =~
~={red}- Very little lost due to mito/chloro filtering (total = 1460)=~
~={red}- over half the features were lost from freq 10 filtering (remaining 732)=~
~={red}- 657 features left after keeping only features in 2 samples.=~

## Generate taxa plots of rarefied tables from the v4v5 taxonomy
```
cd ../
cd taxaplots

qiime feature-table rarefy \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2.qza \
--p-sampling-depth 40000 \
--o-rarefied-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k.qza

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_v4v5_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxaplot_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k.qzv
```

Looks great, and way better with the v4v5 classifier 

filter out control samples from tables
```
cd ../dada2

#filter out controls from the unrarefied table (use for core metrics)

qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_v4v5_GG2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[Tissue] != 'positive_controls'" \
  --o-filtered-table table_noMitoChloro_nb_v4v5_GG2_noControls.qza  

qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[Tissue] != 'positive_controls'" \
  --o-filtered-table table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_noControls.qza  

#filter out controls from the rarefied table
qiime feature-table filter-samples \
  --i-table table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k.qza \
  --m-metadata-file ../metadata/metadata.txt \
  --p-where "[Tissue] != 'positive_controls'" \
  --o-filtered-table table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k_noControls.qza

#### taxa plot rarefied, no controls 
cd ../taxaplots

qiime taxa barplot \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k_noControls.qza \
--i-taxonomy ../taxonomy/taxonomy_nb_v4v5_gg2.qza \
--m-metadata-file ../metadata/metadata.txt \
--o-visualization taxa_plot_table_noMitoChloro_nb_v4v5_GG2_minfreq10_minsample2_rare40k_noControls.qzv
```

## Generate a phylogenetic tree (SEPP tree, gg2)
```
#!/bin/bash
#SBATCH --job-name=gg2_tree_mouse
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
cd /scratch/alpine/lindsval@colostate.edu/sahar_data_2026/tree

wget --no-check-certificate https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.sepp-reference.qza 

qiime fragment-insertion sepp \
--i-representative-sequences ../dada2/mouse_rep_seqs_dada2_filtered.qza \
--i-reference-database 2022.10.backbone.sepp-reference.qza \
--o-tree tree_gg2.qza \
--o-placements tree_placements_gg2.qza \
--p-threads 4
```

```
sbatch tree.sh
```
Submitted batch job 24519162


~={orange}last thing to submit!=~
## Core Metrics
```
cd /scratch/alpine/lindsval@colostate.edu/sahar_data_2026/core_metrics

#core metrics on all samples
qiime diversity core-metrics-phylogenetic \
--i-phylogeny ../tree/tree_gg2.qza \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2.qza \
--p-sampling-depth 40000 \
--m-metadata-file ../metadata/metadata.txt \
--output-dir core_metrics_rare40k_gg2_v4v5

#core metrics for only samples (no controls)
qiime diversity core-metrics-phylogenetic \
--i-phylogeny ../tree/tree_gg2.qza \
--i-table ../dada2/table_noMitoChloro_nb_v4v5_GG2_noControls.qza \
--p-sampling-depth 40000 \
--m-metadata-file ../metadata/metadata.txt \
--output-dir core_metrics_rare40k_gg2_v4v5_noControls
```


## Export alpha and beta diversity files to then run in R
```
cd /scratch/alpine/lindsval@colostate.edu/sahar_data_2026/core_metrics

mkdir export

#shannon
unzip core_metrics_rare40k_gg2_v4v5_noControls/shannon_vector.qza -d export/shannon

# Observed Features  
unzip core_metrics_rare40k_gg2_v4v5_noControls/observed_features_vector.qza -d export/observed_features  
  
# Faith's PD  
unzip core_metrics_rare40k_gg2_v4v5_noControls/faith_pd_vector.qza -d export/faith_pd  
  
# Pielou's evenness  
unzip core_metrics_rare40k_gg2_v4v5_noControls/evenness_vector.qza -d export/evenness

# Bray Curtis  
unzip core_metrics_rare40k_gg2_v4v5_noControls/bray_curtis_pcoa_results.qza -d export/bray_curtis  
  
# Jaccard  
unzip core_metrics_rare40k_gg2_v4v5_noControls/jaccard_pcoa_results.qza -d export/jaccard  
  
# Unweighted Unifrac  
unzip core_metrics_rare40k_gg2_v4v5_noControls/unweighted_unifrac_pcoa_results.qza -d export/unweighted_unifrac  
  
# Weighted Unifrac  
unzip core_metrics_rare40k_gg2_v4v5_noControls/weighted_unifrac_pcoa_results.qza -d export/weighted_unifrac

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

now download and then go to R script

/Users/valerielindstrom/Documents/PostDoc/data_consulting/sahar_mouse_2026/data/core_metrics_export

```
cd /scratch/alpine/lindsval@colostate.edu/sahar_data_2026/core_metrics/core_metrics_rare40k_gg2_v4v5_noControls/rarefied_table/data

biom convert -i feature-table.biom -o feature_table_rare40k.tsv --to-tsv

```