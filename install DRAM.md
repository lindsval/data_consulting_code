

### Install DRAM 1.5 on Alpine

```

cd /projects/lindsval@colostate.edu
ainteractive --ntasks=4 --time=01:00:00 --partition=acpu --qos=cpu-normal
module purge

module load anaconda

wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/environment.yaml

conda config --set solver libmamba

#create env, this took about 15 mins
time conda env create -f environment.yaml -n test_dram_again_sept232026

conda activate test_dram_again_sept232026

```

### once the conda env is created, some manual changes need to happen before we can set up: 
```
# go here (projects/lindsval@colostate.edu/software/anaconda/envs/test_dram_again_sept232026/bin/DRAM-setup.py) and change the name of the file from DRAM-setup.py to DRAM-setup-original.py

# then upload the version of DRAM-setup.py i gave you, put it in your respective bin directory (mine looks like: projects/lindsval@colostate.edu/software/anaconda/envs/test_install_DRAM_v1.5.0_sept2026/bin/DRAM-setup.py)

# change the path in the first line to YOUR projects path
## eg. change from #!/projects/lindsval@colostate.edu/software/anaconda/envs/DRAM_v1.5.0_use/bin/python to #!/projects/YOUR USERNAME/software/anaconda/envs/YOUR ENV/bin/python

chmod +x /projects/lindsval@colostate.edu/software/anaconda/envs/test_dram_again_sept232026/bin/DRAM-setup.py

#then also upload this file ("database_processing_vog_fixed.py") to your environment's mag_annotator folder: (eg. /projects/$USER/software/anaconda/envs/ENVIRONMENT/lib/python3.10/site-packages/mag_annotator/)

#this file is from my original DRAM install /projects/lindsval@colostate.edu/software/anaconda/envs/DRAM_v1.5.0_use/lib/python3.10/site-packages/mag_annotator/

chmod +x /projects/lindsval@colostate.edu/software/anaconda/envs/test_dram_again_sept232026/lib/python3.10/site-packages/mag_annotator/database_processing_vog_fixed.py

#then run 
module load anaconda
conda activate test_dram_again_sept232026

#need to downgrade the setuptools version within the DRAM environment
pip install "setuptools<81"

#copy the folder "preformatted_databases_from_kayla" to your scratch directory, then run set up
DRAM-setup.py import_config --config_loc  /scratch/alpine/lindsval@colostate.edu/preformatted_databases_from_kayla/CONFIG

#check set up worked
DRAM-setup.py print_config
#this should list out your database paths

# run DRAM (example)
DRAM.py annotate -i '*fa' -o  DRAM_1.5_09082026 --min_contig_size 2500 --threads 20
DRAM.py distill -i DRAM_1.5_09082026/annotations.tsv -o DRAM_1.5_09082026/distill


```


### and if you want to be able to use annotate_genes do

```
conda activate test_dram_again_sept232026

cp "$CONDA_PREFIX/lib/python3.10/site-packages/mag_annotator/annotate_bins.py" \
"$CONDA_PREFIX/lib/python3.10/site-packages/mag_annotator/annotate_bins.py.bak" && \
sed -i "s/keep_tmp_dir=True, low_mem_mode=False, threads=10, verbose=True):/keep_tmp_dir=True, low_mem_mode=False, threads=10, verbose=True, config_loc=None):/" \
"$CONDA_PREFIX/lib/python3.10/site-packages/mag_annotator/annotate_bins.py" && \
sed -i "s/rename_genes, keep_tmp_dir, low_mem_mode, threads, verbose)/rename_genes, keep_tmp_dir, low_mem_mode, threads, verbose, log_file_path, config_loc)/" \
"$CONDA_PREFIX/lib/python3.10/site-packages/mag_annotator/annotate_bins.py"

sed -i "s/custom_hmm_cutoffs_loc, use_uniref, use_camper, use_fegenie, /custom_hmm_cutoffs_loc, use_uniref, use_camper, /; s/use_sulphur, use_vogdb/use_vogdb/" "$CONDA_PREFIX/lib/python3.10/site-packages/mag_annotator/annotate_bins.py"

```