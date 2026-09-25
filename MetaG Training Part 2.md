Part 2: Zip data, assembly, binning, QC, dreplicate to build a MAG database, annotate bins for taxonomy (GTDB) and metabolisms (DRAM)

----------------------------------------------------------------

## Steps
 Step 1: zip raw reads after quality checking
 
 Step 2: install megahit and assemble reads (individual assembly)
 
 Step 3: Create contigs stats file & run it to generate contig stats
 
 Step 4: Combine all contig stats files
 
Step 5: Run Co-Assembly (includes concatenating the R1 and R2 files and running megahit for those concat files)

Step 6: Get co-assembly stats

Step 7:  Pull out contigs >2.5kb for binning

Step 8: Map trimmed paired-end reads back to ≥2500 bp assembled contigs to generate coverage information for MAG binning and abundance estimation using bbmap

Step 9: Binning with Metabat (v. 2:2.18)

Step 10: Run checkM on the bins to get quality and completness information 

Step 11: Dereplicate the M/HQ MAGs

Step 12: Make genome db of the 95%id mapped MAGs and map reads back to DB to see % reads mapped

Step 13: Run DRAM on the M/HQ MAGs

Step 14:  Run GTDB-tk for MAG taxonomy


----------------------------------------------------------------
##  Step 1: zip or delete raw reads
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

##  Step 2: install MEGAHIT 

```
acompile --ntasks=4 --time=03:00:00
module load anaconda
conda create -n megahit
conda activate megahit
conda install -c bioconda megahit
megahit -v
#this should print: MEGAHIT v1.2.9
```

##  Step 2b: Assemble trimmed reads using MEGAHIT (v1.2.9)

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

```


## Step 3: Create contigs stats .pl file, save this in a new directory called custom_scripts
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

## Step 3b: run contig_stats on all assemblies

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

#### Check this ran for all samples

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


## Step 4: Combine all contig stats files
Copy this into a new .sh file and then run it 
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

after running the contig stats, you then need to export that txt file and in excel, calculate the % reads in each range. in particular we are interesting the % contigs in the following categories:
|100-500_reads_pct|
|500-1000_reads_pct|
|1000-5000_reads_pct|

NOTE: 
- remember that we care about those contigs that are in the 1000-5000bp bucket because we only bin/annotate contigs that are greater than 2.5kb (since we want to use good-quality contigs and 2.5kb will account for about 2-3 genes, so we want to use contigs with more than just 1 gene assembled, this will help with high quality bins)
- Assembly is bad if there is no assembly at all, sometimes you will just get low % in the 1-5kb bucket, you can try altnerative assemblers like IDBA-UD if you have a small contigs problem (see Part 3 of metaG training)
- you also want to look at the N50, this is a metric of contig legnth
- See also how many assembled contigs you have: do this by counting the sequences in the final contigs file `grep -c ">" final.contigs.fa.`
- see also what your longest contig is: 
```
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG
while read sample; do
    contigs="${sample}/assembly/megahit_out2/final.contigs.fa"
    if [ -f "$contigs" ]; then
        longest=$(awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq $0} END {if (seq) print length(seq)}' "$contigs" | sort -nr | head -1)
        echo -e "${sample}\t${longest}"
    else
        echo -e "${sample}\tFILE_NOT_FOUND"
    fi
done < sample_list.txt > longest_contigs.txt
```

## Step 5: Run Co-Assembly
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
mkdir DroughtDeligeBulk
mkdir Control
```
#### Create a sample list for the coassembly samples
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

### Create subdirectories (for each coassembly we will create new subdirectories)
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
This took about 3 days, but the assembly wasn't amazing, so a better dataset may need a lot more time/memory

## Step 6: Get co-assembly stats


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
MAX_JOBS=1
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
#SBATCH --cpus-per-task=1
#SBATCH --mem=80gb 
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
Submitted batch job 27469417 
#### Combine all contig stats files

```
#!/bin/bash

BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"
SAMPLE_LIST="${BASE_DIR}/coA_sample_list.txt"
OUTFILE="${BASE_DIR}/all_coassembly_contig_stats_summary.txt"

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
bash 11a_combine_stats_coA.sh

after running the contig stats, you then need to export that txt file and in excel, calculate the % reads in each range. in particular we are interesting the % contigs in the following categories:
|100-500_reads_pct|
|500-1000_reads_pct|
|1000-5000_reads_pct|


## Step 7: extract contigs >2.5kb using pullseqs

#### Compile pullseqs program

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
```
### Run pullseqs for individual assembly 

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

### Run pullseqs for coassembly
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


## Step 8: Map trimmed paired-end reads back to ≥2500 bp assembled contigs to generate coverage information for MAG binning and abundance estimation using bbmap

### Make folder for mapped reads

```
while read SAMPLE; do  
mkdir -p /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/${SAMPLE}/mapped_reads  
done < /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt

while read SAMPLE; do  
mkdir -p /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/${SAMPLE}/mapped_reads  
done < /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt
```

### Individual assembly mapping
- uses the contigs as the reference to get coverage
- compares every paired end read from every sample to see where the reads map to the contigs (gives position information)

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
sbatch 13_bbmap_indiv.sh
Submitted batch job 27469328 (took about 4 hours)

### CoAssembly mapping

```
#!/bin/bash
#SBATCH --job-name=bbmap_CoAssembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=240gb
#SBATCH --partition=amem
#SBATCH --qos=mem-long
#SBATCH --time=48:00:00
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
Submitted batch job 27548717

### Convert SAM to BAM files, sort, filter
- BAM file is sorted based on its position in the reference, as determined by its alignment

```
#!/bin/bash  
#SBATCH --job-name=indiv_sort  
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20  
#SBATCH --time=23:00:00  
#SBATCH --mem=20gb  
#SBATCH --qos=normal  
#SBATCH --partition=amilan  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=lindsval@colostate.edu  
#SBATCH --output=slurm_output/indiv_sort%j.out  
#SBATCH --error=slurm_output/indiv_sort%j.err  

module load samtools  
  
SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"  
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"  
  
while read SAMPLE; do

    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    OUTPUT="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sam"

    echo "Processing sample: $SAMPLE"

    if [[ -f "$OUTPUT" ]]; then

        samtools view -@ 20 -bS "$OUTPUT" \
            > "${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.bam"

        samtools sort -@ 20 \
            -T "${MAPPED_DIR}/${SAMPLE}_tmp_sort" \
            -o "${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sorted.bam" \
            "${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.bam"

        echo "Finished: $SAMPLE"

    else
        echo "Missing SAM file: $OUTPUT"
    fi

done < "$SAMPLE_LIST"

echo "All samples complete."

```

14_sort_indiv_assembly.sh
Submitted batch job 27550775

```
#!/bin/bash
#SBATCH --job-name=indiv_sort_remaining_samples
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=23:00:00
#SBATCH --mem=20gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/indiv_sort%j.out
#SBATCH --error=slurm_output/indiv_sort%j.err

module load samtools

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

while read SAMPLE; do

    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    SAM_FILE="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sam"
    BAM_FILE="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.bam"
    SORTED_BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sorted.bam"

    echo "Checking sample: $SAMPLE"

    # Skip samples that already finished
    if [[ -f "$SORTED_BAM" ]]; then
        echo "Sorted BAM already exists. Skipping $SAMPLE"
        continue
    fi

    # Make sure SAM file exists
    if [[ ! -f "$SAM_FILE" ]]; then
        echo "Missing SAM file: $SAM_FILE"
        continue
    fi

    echo "Processing $SAMPLE"

    samtools view -@ 20 -bS "$SAM_FILE" \
        > "$BAM_FILE"

    samtools sort -@ 20 \
        -T "${MAPPED_DIR}/${SAMPLE}_tmp_sort" \
        -o "$SORTED_BAM" \
        "$BAM_FILE"

    echo "Finished: $SAMPLE"

done < "$SAMPLE_LIST"

echo "All samples complete."
```
14_sort_indiv_assembly_remaining_samples.sh

```
#!/bin/bash  
#SBATCH --job-name=coA_sort  
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20 
#SBATCH --partition=amem
#SBATCH --qos=mem-long
#SBATCH --time=48:00:00 
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=lindsval@colostate.edu  
#SBATCH --output=slurm_output/coA_sort%j.out  
#SBATCH --error=slurm_output/coA_sort%j.err  
  
module load samtools  
  
SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"  
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"  
  
while read SAMPLE; do

    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    OUTPUT="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sam"

    echo "Processing sample: $SAMPLE"

    if [[ -f "$OUTPUT" ]]; then

        samtools view -@ 20 -bS "$OUTPUT" \
            > "${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.bam"

        samtools sort -@ 20 \
            -T "${MAPPED_DIR}/${SAMPLE}_tmp_sort" \
            -o "${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sorted.bam" \
            "${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.bam"

        echo "Finished: $SAMPLE"

    else
        echo "Missing SAM file: $OUTPUT"
    fi

done < "$SAMPLE_LIST"

echo "All samples complete."
```
14_sort_CoAssembly.sh
Submitted batch job 28213051

### Reformat individual assembly files
BBTools reformat.sh script will filter to keep only the highest quality matches

```

#!/bin/bash
#SBATCH --job-name=indiv_reformat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=04:00:00
#SBATCH --mem=120gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/indiv_reformat%j.out
#SBATCH --error=slurm_output/indiv_reformat%j.err

module load anaconda  
module load samtools
module load bbtools
conda activate bbmap
echo "reformat: $(which reformat.sh)"
echo "samtools: $(which samtools)"
reformat.sh --version
samtools --version | head -n 1

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

while read SAMPLE; do
    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    INPUT_BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sorted.bam"
OUTPUT_BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped99per.sorted.bam"
    echo "Processing sample: $SAMPLE"
    if [[ -f "$INPUT_BAM" ]]; then
        reformat.sh \
            -Xmx100g \
            threads=20 \
            minidfilter=0.99 \
            in="$INPUT_BAM" \
            out="$OUTPUT_BAM" \
            pairedonly=t \
            primaryonly=t \
            overwrite=true
        echo "Finished: $SAMPLE"
    else
        echo "Missing BAM file: $INPUT_BAM"
    fi
done < "$SAMPLE_LIST"
echo "All samples complete."

```
15_reformat_indiv.sh
Submitted batch job 28402102, rerunning june 16th, since the `idfilter` flag was updated to `minidfilter`, DONE

### Reformat coassembly files

```

#!/bin/bash
#SBATCH --job-name=coA_reformat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=04:00:00
#SBATCH --mem=120gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/coA_reformat%j.out
#SBATCH --error=slurm_output/coA_reformat%j.err

module load anaconda  
module load samtools
module load bbtools
conda activate bbmap
echo "reformat: $(which reformat.sh)"
echo "samtools: $(which samtools)"
reformat.sh --version
samtools --version | head -n 1

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"

while read SAMPLE; do
    MAPPED_DIR="${BASE_DIR}/${SAMPLE}/mapped_reads"
    INPUT_BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped.sorted.bam"
OUTPUT_BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped99per.sorted.bam"
    echo "Processing sample: $SAMPLE"
    if [[ -f "$INPUT_BAM" ]]; then
        reformat.sh \
            -Xmx100g \
            threads=20 \
            minidfilter=0.99 \
            in="$INPUT_BAM" \
            out="$OUTPUT_BAM" \
            pairedonly=t \
            primaryonly=t \
            overwrite=true
        echo "Finished: $SAMPLE"
    else
        echo "Missing BAM file: $INPUT_BAM"
    fi
done < "$SAMPLE_LIST"
echo "All samples complete."

```
15_reformat_coA.sh


## Step 9: Binning with Metabat (v. 2:2.18)

### Install metabat

```
acompile --ntasks=4 
module load anaconda
conda config --set channel_priority strict
conda create -n metabat2 \
    -c conda-forge \
    -c bioconda \
    boost-cpp=1.85 \
    boost=1.85 \
    metabat2
conda activate metabat2

```


### Co assembly binning
```
#!/bin/bash
#SBATCH --job-name=coA_metabat_bin
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=23:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/coA_metabat_bin%j.out
#SBATCH --error=slurm_output/coA_metabat_bin%j.err

module load anaconda
conda activate metabat2

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"

while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do

    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue

    SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
    ASSEMBLY_DIR="${SAMPLE_DIR}/assembly/megahit_out"
    MAPPED_DIR="${SAMPLE_DIR}/mapped_reads"

    CONTIGS="${ASSEMBLY_DIR}/${SAMPLE}_final.contigs_2500.fa"
    BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped99per.sorted.bam"

    OUT_DIR="${SAMPLE_DIR}/metabat_bins"
    mkdir -p "$OUT_DIR"
    if [[ ! -f "$CONTIGS" || ! -f "$BAM" ]]; then
        echo "Skipping $SAMPLE (missing input)"
        continue
    fi
    echo "Processing $SAMPLE"
    DEPTH="${ASSEMBLY_DIR}/depth.txt"
    jgi_summarize_bam_contig_depths \
        --outputDepth "$DEPTH" \
        "$BAM"
    metabat2 \
        -i "$CONTIGS" \
        -a "$DEPTH" \
        -o "$OUT_DIR/${SAMPLE}_bin" \
        -t 6
done < "$SAMPLE_LIST"
echo "All samples complete."
```
16_coA_binning.sh
Submitted batch job 28414634, done
### Individual assembly binning
```
#!/bin/bash
#SBATCH --job-name=indiv_metabat_bin
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=23:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/indiv_metabat_bin%j.out
#SBATCH --error=slurm_output/indiv_metabat_bin%j.err

module load anaconda
conda activate metabat2

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do

    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue

    SAMPLE_DIR="${BASE_DIR}/${SAMPLE}"
    ASSEMBLY_DIR="${SAMPLE_DIR}/assembly/megahit_out"
    MAPPED_DIR="${SAMPLE_DIR}/mapped_reads"

    CONTIGS="${ASSEMBLY_DIR}/${SAMPLE}_final.contigs_2500.fa"
    BAM="${MAPPED_DIR}/${SAMPLE}_final.contigs_2500_mapped99per.sorted.bam"

    OUT_DIR="${SAMPLE_DIR}/metabat_bins"
    mkdir -p "$OUT_DIR"
    if [[ ! -f "$CONTIGS" || ! -f "$BAM" ]]; then
        echo "Skipping $SAMPLE (missing input)"
        continue
    fi
    echo "Processing $SAMPLE"
    DEPTH="${ASSEMBLY_DIR}/depth.txt"
    jgi_summarize_bam_contig_depths \
        --outputDepth "$DEPTH" \
        "$BAM"
    metabat2 \
        -i "$CONTIGS" \
        -a "$DEPTH" \
        -o "$OUT_DIR/${SAMPLE}_bin" \
        -t 6
done < "$SAMPLE_LIST"
echo "All samples complete."

```

16_indiv_bining.sh
Submitted batch job 28414617, done

### Check the number of bins generated and record in a new file

```
SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"
OUTFILE="${BASE_DIR}/metabat_bin_counts.txt"

echo -e "Sample\tBin_count" > "$OUTFILE"
total=0
while read sample; do
    BIN_DIR="${BASE_DIR}/${sample}/metabat_bins"
    if [[ -d "$BIN_DIR" ]]; then
        count=$(find "$BIN_DIR" -maxdepth 1 -name "*.fa" | wc -l)
        total=$((total + count))
    else
        count="NA"
    fi
    echo -e "${sample}\t${count}" >> "$OUTFILE"
done < "$SAMPLE_LIST"
echo -e "TOTAL\t${total}" >> "$OUTFILE"
echo "Results written to $OUTFILE"

# 310 individual assembly bins


SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"
OUTFILE="${BASE_DIR}/metabat_bin_counts_coassembly.txt"

echo -e "Sample\tBin_count" > "$OUTFILE"
total=0
while read sample; do
    BIN_DIR="${BASE_DIR}/${sample}/metabat_bins"
    if [[ -d "$BIN_DIR" ]]; then
        count=$(find "$BIN_DIR" -maxdepth 1 -name "*.fa" | wc -l)
        total=$((total + count))
    else
        count="NA"
    fi
    echo -e "${sample}\t${count}" >> "$OUTFILE"
done < "$SAMPLE_LIST"
echo -e "TOTAL\t${total}" >> "$OUTFILE"
echo "Results written to $OUTFILE"
```

## Step 10: Run checkM on the bins to get quality and completness information
### checkM v1.2.3 install

```

acompile --ntasks=4 
module load anaconda
conda create -n checkm
conda activate checkm
conda install -c bioconda -c conda-forge checkm-genome
checkm -h


```
### Run checkM on individual assembly bins
```
#!/bin/bash
#SBATCH --job-name=indiv_checkM
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=23:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/indiv_checkM%j.out
#SBATCH --error=slurm_output/indiv_checkM%j.err

module load anaconda
conda activate checkm

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"

while read sample; do
    echo "Processing ${sample}..."
    BIN_DIR="${BASE_DIR}/${sample}/metabat_bins"
    if [[ ! -d "${BIN_DIR}" ]]; then
        echo "Skipping ${sample}: metabat_bins not found."
        continue
    fi
    cd "${BIN_DIR}"
    # Run CheckM
    checkm lineage_wf \
        -t 12 \
        -x fa \
        . \
        checkm
    # Generate QA table
    checkm qa \
        -o 2 \
        -f checkm/results.txt \
        --tab_table \
        -t 12 \
        checkm/lineage.ms \
        checkm
done < "${SAMPLE_LIST}"


```
sbatch 17_indiv_checkM.sh
Submitted batch job 28449358, done

### pull out the HQ/MQ bins from the checkM results on individual assembly 

```
SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"
OUTDIR="${BASE_DIR}/checkm_filtered_bins_individual"
mkdir -p "$OUTDIR"
while read sample; do
    echo "Processing ${sample}..."
    CHECKM_FILE="${BASE_DIR}/${sample}/metabat_bins/checkm/results.txt"
    OUTFILE="${OUTDIR}/${sample}_HQ_MQ_bins.txt"
    if [[ ! -f "$CHECKM_FILE" ]]; then
        echo -e "${sample}\tNO_CHECKM_FILE" > "$OUTFILE"
        continue
    fi
    echo -e "bin\tquality" > "$OUTFILE"
    awk -F "\t" '
    NR>1 {
        if ($6 >= 90 && $7 <= 5) {
            print $1 "\tHIGH"
        }
        else if ($6 >= 50 && $7 <= 10) {
            print $1 "\tMEDIUM"
        }
    }' "$CHECKM_FILE" >> "$OUTFILE"
done < "$SAMPLE_LIST"
```

```
nano extract_mqhq_bins_indiv.sh
chmod +x extract_mqhq_bins_indiv.sh 
./extract_mqhq_bins_indiv.sh 
```

### Count the number of bins from the individual assembly / checkM that were medium- or high-quality
```
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/checkm_filtered_bins_individual" 
grep -h -E "HIGH|MEDIUM" ${OUTDIR}/*_HQ_MQ_bins.txt | wc -l
 
```

### Run checkM on coassembly assembly bins
```
#!/bin/bash
#SBATCH --job-name=coA_checkM
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=23:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/coA_checkM%j.out
#SBATCH --error=slurm_output/coA_checkM%j.err

module load anaconda
conda activate checkm

SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"

while read sample; do
    echo "Processing ${sample}..."
    BIN_DIR="${BASE_DIR}/${sample}/metabat_bins"
    if [[ ! -d "${BIN_DIR}" ]]; then
        echo "Skipping ${sample}: metabat_bins not found."
        continue
    fi
    cd "${BIN_DIR}"
    # Run CheckM
    checkm lineage_wf \
        -t 12 \
        -x fa \
        . \
        checkm
    # Generate QA table
    checkm qa \
        -o 2 \
        -f checkm/results.txt \
        --tab_table \
        -t 12 \
        checkm/lineage.ms \
        checkm
done < "${SAMPLE_LIST}"


```
sbatch 17_coA_checkM.sh
### pull out the HQ/MQ bins

```
SAMPLE_LIST="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/coA_sample_list.txt"
BASE_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly"
OUTDIR="${BASE_DIR}/checkm_filtered_bins"
mkdir -p "$OUTDIR"
while read sample; do
    echo "Processing ${sample}..."
    CHECKM_FILE="${BASE_DIR}/${sample}/metabat_bins/checkm/results.txt"
    OUTFILE="${OUTDIR}/${sample}_HQ_MQ_bins.txt"
    if [[ ! -f "$CHECKM_FILE" ]]; then
        echo -e "${sample}\tNO_CHECKM_FILE" > "$OUTFILE"
        continue
    fi
    echo -e "bin\tquality" > "$OUTFILE"
    awk -F "\t" '
    NR>1 {
        if ($6 >= 90 && $7 <= 5) {
            print $1 "\tHIGH"
        }
        else if ($6 >= 50 && $7 <= 10) {
            print $1 "\tMEDIUM"
        }
    }' "$CHECKM_FILE" >> "$OUTFILE"
done < "$SAMPLE_LIST"
```

```
nano extract_mqhq_bins_coA.sh
chmod +x extract_mqhq_bins_coA.sh 
./extract_mqhq_bins_coA.sh 
```

### Count the number of bins from the coassembly / checkM that were medium- or high-quality
```
OUTDIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/coassembly/checkm_filtered_bins"  
grep -h -E "HIGH|MEDIUM" ${OUTDIR}/*_HQ_MQ_bins.txt | wc -l

```


## Step 11: dereplicate the M/HQ MAGs

### first move all of the medium and high quality bins to a new directory (eg. MedHighQualityMAGs)
#### Install dRep
```
acompile --ntasks=4 
module load miniforge
mamba create -n drep -c conda-forge -c bioconda drep
mamba activate
dRep -h #version 3.7.1
dRep check_dependencies 

```

#### Run dRep on the M/HQ database 
```
#!/bin/bash
#SBATCH --job-name=dRep
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/dRep_%j.out
#SBATCH --error=slurm_output/dRep_%j.err

module load miniforge
mamba activate drep

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs
dRep dereplicate dRep_bins -p 7 -comp 50 -con 10 -g ./*fa
```

dRep.sh 
Submitted batch job 31041133
#### Count number of bins before and after dereplication - 29 left after dRep
```

ls -dq *fa | wc -l
# there were 121 M and HQ bins before dRep

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins/dereplicated_genomes
#29 left afte dRep
```


## Step 12: Make genome db of the 95%id mapped MAGs and map reads back to DB to see % reads mapped


### First, rename the bin headers to match dram naming style where the scaffolds are now named with the bin name appended at the beginning. the bin headers must match. 
```
#!/bin/bash
INPUT_DIR="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins/dereplicated_genomes"
OUTPUT_DIR="${INPUT_DIR}/genomes_renamed"
mkdir -p "$OUTPUT_DIR"
for fasta in "$INPUT_DIR"/*.fa; do
    [ -e "$fasta" ] || continue
    filename=$(basename "$fasta")
    bin="${filename%.fa}"
    echo "Renaming $filename..."
    awk -v prefix="$bin" '
    /^>/ {
        sub(/^>/, "")
        print ">" prefix "_" $0
        next
    }
    { print }
    ' "$fasta" > "${OUTPUT_DIR}/${filename}"
done
echo "Done."
```


```
bash rename_bins_like_dram.sh
```


Concatenate MAG database and map concatenated reads metagenomic reads
```
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins/dereplicated_genomes/genome_renamed
cat *fa > cat_MAGs_dRep.fa
```

```
#loop over trimmed files and concat them in 1 r1 file and 1 r2 file:

#!/bin/bash
#SBATCH --job-name=concat_reads
#SBATCH --partition=acpu
#SBATCH --qos=c
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/concat_reads_%j.out
#SBATCH --error=slurm_output/concat_reads_%j.err 

BASE="/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG"
SAMPLE_LIST="${BASE}/sample_list.txt"
OUTDIR="${BASE}/concat_reads"

mkdir -p "$OUTDIR"

while read -r sample; do
    cat "${BASE}/${sample}/processed_reads/${sample}_R1_bbduktrimmed.fastq"
done < "$SAMPLE_LIST" > "${OUTDIR}/concat_R1_bbduktrimmed.fastq"

while read -r sample; do
    cat "${BASE}/${sample}/processed_reads/${sample}_R2_bbduktrimmed.fastq"
done < "$SAMPLE_LIST" > "${OUTDIR}/concat_R2_bbduktrimmed.fastq"
```
sbatch concat_reads.sh
Submitted batch job 31041887, done
### mapping reads to MAG db
```
#!/bin/bash
#SBATCH --job-name=map_reads
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=23:30:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/map_reads_%j.out
#SBATCH --error=slurm_output/map_reads_%j.err 

module load anaconda
conda activate bbmap

cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins/dereplicated_genomes/genome_renamed

bbmap.sh -Xmx48G threads=16 overwrite=t ref=cat_MAGs_dRep.fa in1=../../../../concat_reads/concat_R1_bbduktrimmed.fastq in2=../../../../concat_reads/concat_R2_bbduktrimmed.fastq outm=mapped_reads_interleaved.fastq #reads that mapped to the genome database

```

sbatch `map_reads_to_MAGdb.sh`, Submitted batch job 31044297

 calculate % reads mapped
```
#this will count the number of lines in the file we created above for those that mapped to the genome database, then divide by 4 to get the number of reads that mapped
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins/dereplicated_genomes/genome_renamed
wc -l mapped_reads_interleaved.fastq

##this will count the number of lines in the R1 file, then divide by 4 to get the number of reads you started with
wc -l concat_R1_bbduktrimmed.fastq

#then in excel, divide the # of reads mapped by the # of reads total to get the % of reads mapping 
```

see example calculation:  "/Users/valerielindstrom/Documents/PostDoc/data_consulting/roberts_metagenomics/misc_files/map_to_MAGs__reads_mapped.xlsx"

looking for >20% at the bare minimum 

## Step 13: Run DRAM on the M/HQ MAGs
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

### Run DRAM on the dreplicated M/HQ MAGs 

```
#!/bin/bash
#SBATCH --job-name=DRAM_50mags_use
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=20
#SBATCH --time=23:30:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --output=slurm_output/DRAM_50mags_use_%j.out
#SBATCH --error=slurm_output/DRAM_50mags_use_%j.err 


module load anaconda
conda activate DRAM_v1.5.0_use
cd /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes

DRAM.py annotate -i '*fa' -o  DRAM_1.5_09092026 --min_contig_size 2500 --threads 20
DRAM.py distill -i DRAM_1.5_09092026/annotations.tsv -o DRAM_1.5_09092026/distill

```
sbatch DRAM_50mags_use.sh
Submitted batch job 32335874



## Step 14: Run GTDB-tk for MAG taxonomy

#### Install gtdbk-tk - requires large db so installing into scratch - Release 11-RS232 (15th April 2026)
```
acompile --ntasks=4 --time=03:00:00
module load miniforge
mamba create -n gtdbtk-2.7.2 -c conda-forge -c bioconda gtdbtk=2.7.2
mamba activate gtdbtk-2.7.2
gtdbtk --version

mkdir -p /scratch/alpine/lindsval@colostate.edu/gtdbtk_r232
```

```
#!/bin/bash
#SBATCH --job-name=gtdbtk_db
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --nodes=1
#SBATCH --output=gtdbtk_db_%j.out
#SBATCH --error=gtdbtk_db_%j.err

module load miniforge
mamba activate gtdbtk-2.7.2

cd /scratch/alpine/lindsval@colostate.edu/gtdbtk_r232
wget https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
tar -xvzf gtdbtk_r232_data.tar.gz \
    --strip 1 > /dev/null

#rm gtdbtk_r232_data.tar.gz
```
sbatch `download_gtdbtk.sh`
Submitted batch job 31046012, done

After it finishes, set the database path
```
acompile --ntasks=4 --time=03:00:00
module load miniforge
mamba activate gtdbtk-2.7.2
mamba env config vars set \
GTDBTK_DATA_PATH="/scratch/alpine/lindsval@colostate.edu/gtdbtk_r232"
```
Then **deactivate and reactivate**:
```
mamba deactivate
mamba activate gtdbtk-2.7.2
```
Check:

```
echo $GTDBTK_DATA_PATH
```

Verify
```
gtdbtk check_install #everything looks good
```


### run gtdb 
```
#!/bin/bash
#SBATCH --job-name=gtdbtk_30bins
#SBATCH --partition=amem
#SBATCH --qos=mem-normal
#SBATCH --ntasks=48
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lindsval@colostate.edu
#SBATCH --nodes=1
#SBATCH --output=slurm_output/gtdbtk_30bins_%j.out
#SBATCH --error=slurm_output/gtdbtk_30bins_%j.err

module load miniforge
mamba activate gtdbtk-2.7.2

gtdbtk classify_wf \
    --genome_dir /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes \
    --extension fa \
    --out_dir /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/gtdbtk_30bins \
    --cpus 48
```

sbatch gtdb.sh
Submitted batch job 32373200


make it readable for excel
```
cp /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/gtdbtk_30bins/gtdbtk.ar53.summary.tsv \
/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/gtdbtk_30bins/gtdbtk_ar53_summary_for_excel.tsv

cp /scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/gtdbtk_30bins/gtdbtk.bac120.summary.tsv \
/scratch/alpine/lindsval@colostate.edu/roberts_soils_metaG/MedHighQualityMAGs/dRep_bins_151/dereplicated_genomes/gtdbtk_30bins/gtdbtk_bac120_summary_for_excel.tsv
```

