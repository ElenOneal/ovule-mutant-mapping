#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 sequences.txt DIR read_dir genome_dir genome"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
read_dir="$3"
genome_dir="$4"
genome="$5"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r f c a b d; do
  script_file="$a.align.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.align.out"
    echo "#SBATCH --error=${a}.align.err"
    echo '#SBATCH --cpus-per-task=6'
    echo '#SBATCH --account=biodept'
    echo '#SBATCH -p common,biodept'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=24G'
    echo ''
    echo "bwa mem -t 6 -M $genome_dir/$genome $read_dir/$a.PE.R1.fq.gz $read_dir/$a.PE.R2.fq.gz | samtools view -Shb | samtools sort -T $a.sort -o $a.sort.bam"
    echo "samtools index $a.sort.bam"
    echo "sbatch $a.filter.sh"
  } > "$script_file"
done < "$file_name"
``

