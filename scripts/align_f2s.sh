#!/bin/bash
set -euo pipefail

# Create shell scripts to trim reads and align to genome


# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 sample_names DIR genome_dir genome partition"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
genome_dir="$3"
genome="$4"
partition="$5"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

mkdir -p "$DIR"

# sample_names file format:
#  <barcode>    <assigned_F2_name>    <external_F2_name>
# Example:
#  GGAAGAA    F2.01    1A
#  CAGAGAA    F2.02    6E
#
# - first column is the barcode sequence
# - second column is the F2 name assigned in this project
# - third column is the F2 name used by the other lab
# - fields are tab-separated

# Main loop
while IFS=$'\t' read -r a b c; do
  script_file="$b.align.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --job-name=$b"
    echo "#SBATCH --output=${b}.align.out"
    echo "#SBATCH --error=${b}.align.err"
    echo '#SBATCH --cpus-per-task=4'
    echo "#SBATCH -p $partition"
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=4G'
    echo ''
    echo "trimmomatic PE  -threads 4 -phred33 -trimlog $b.trimlog -quiet -validatePairs $b.1.fq.gz $b.2.fq.gz $b.PE.R1.fq.gz $b.U.R1.fq.gz $b.PE.R2.fq.gz $b.U.R2.fq.gz ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 AVGQUAL:25"
    echo "bwa mem -t 4 $genome_dir/$genome $b.PE.R1.fq.gz $b.PE.R2.fq.gz | samtools view -Shb | samtools sort -T $b.sort -o $b.sort.bam"
  } > "$script_file"
done < "$file_name"


