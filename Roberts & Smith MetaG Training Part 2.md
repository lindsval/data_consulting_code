Part 2: Zip data, assembly, binning

these training docs are just modified (for simplicity) from the data analysis script "Roberts & Smith MetaG Data Processing"

Training Doc: https://docs.google.com/presentation/d/1hXJD7P0WjOE-iaKu9HhdtIVDVLCY8hIjD71KNBJ04WY/edit?slide=id.g3cd05613a6c_0_62#slide=id.g3cd05613a6c_0_62

----------------------------------------------------------------

## Steps
 Step 1: zip raw reads after quality checking
 
 Step 2: install megahit and assemble reads (individual assembly)
 
 Step 3: Create contigs stats file & run it to generate contig stats
 
 Step 4: Combine all contig stats files
 
Step 5: Run Co-Assembly

Step 6: 

Step 7: 

Step 8:

Step 9: 

Step 10: 

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

## Step 5: Run Co-Assembly