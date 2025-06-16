#!/bin/bash

# Create shell scripts to trim reads and align to genome


# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 sample_names DIR genome_dir genome"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
genome_dir="$3"
genome="$4"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r a b c d; do
  script_file="$b.align.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$b"
    echo "#SBATCH --output=${b}.align.out"
    echo "#SBATCH --error=${b}.align.err"
    echo '#SBATCH --cpus-per-task=4'
    echo '#SBATCH --mail-type=ALL'
    echo '#SBATCH --mail-user=eo22@duke.edu'
    echo '#SBATCH -p common,scavenger'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=4G'
    echo ''
    echo "trimmomatic PE  -threads 4 -phred33 -trimlog $b.trimlog -quiet -validatePairs $b.1.fq.gz $b.2.fq.gz $b.PE.R1.fq.gz $b.U.R1.fq.gz $b.PE.R2.fq.gz $b.U.R2.fq.gz ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 AVGQUAL:25"
    echo "bwa mem -t 4 $genome_dir/$genome $b.PE.R1.fq.gz $b.PE.R2.fq.gz | samtools view -Shb | samtools sort -T $b.sort -o $b.sort.bam"
    echo "sbatch ../$b.rg.sh"
  } > "$script_file"
done < "$file_name"
``

