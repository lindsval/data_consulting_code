figure out what all these files mean...

by unzipping the fastq file and coverting just to txt, that gives us: 

### What the FASTQ header tells us

|Header field|Meaning|
|---|---|
|@28c5a133-c869-4dfb-a972-c83be72722fc|Nanopore read ID|
|RG:Z:8019e23e-..._dna_r10.4.1_e8.2_400bps_hac@v5.2.0_barcode15|Sequenced on an ONT R10.4.1 flow cell, high-accuracy (HAC) basecalling, assigned to barcode 15|
|PU:Z:FBG13555|Run or instrument identifier|
|LB:Z:260609_AFK_18S|Library name; this was an 18S run|
|SM:Z:barcode15|Sample/barcode assignment|

The sequence itself starts with your forward primer motif:

AACCTGGTTGATCCTGCCAGT

which is the EukA primer sequence, confirming the amplicon target.

### Why there are 25 FASTQ files

1. Nanopore software writes reads in chunks: it does not usually create one giant FASTQ per sample, so files ending in _0.fastq.gz, _1.fastq.gz, and so on are just consecutive chunks of the same sample.
    
2. You should treat them as one dataset: concatenate them before analysis or point your software at all files together.
    

Example:

Concatenate all chunks for a sample

That gives you one FASTQ containing all reads for that sample.

### Amplicon length

- Your primers are expected to amplify ~1750–1800 bp of the 18S rRNA gene.
    
- The read shown is roughly that size and ends right before the reverse primer motif AGGTGAACCTG from EukB.
    
- So your sequencing facility likely generated full-length 18S Nanopore reads.

### Concatenate reads: 

```
ainteractive --ntasks=6 --time=03:00:00


```


```

mkdir concatenated_fastqs

for sample in */; do
    sample_name=$(basename "$sample")
    cat "$sample"/*.fastq.gz > "concatenated_fastqs/${sample_name}.fastq.gz"
done

```

After concatenation, verify the files:

```
ls -lh concatenated_fastqs/
```

and check read counts:

```
for f in concatenated_fastqs/*.fastq.gzdo    echo "$f"    echo "Reads: $(zcat "$f" | wc -l | awk '{print $1/4}')"done
```

```
qiime tools import \ 
--type 'SampleData[SequencesWithQuality]' \
--input-path manifest.tsv \
--output-path nanopore18S.qza \
--input-format SingleEndFastqManifestPhred33V2
```


```
qiime demux summarize \
--i-data nanopore18S.qza \
--o-visualization nanopore18S.qzv
```