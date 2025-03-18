#!/bin/bash

# Create shell scripts to trim first 2 reads from un-aggregated ddRAD library files
# Uses trimmomatic
# Uses process_radtags From Stacks
# http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php


# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 library_filenames DIR"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r a; do
  script_file="$a.trim.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.trim.out"
    echo "#SBATCH --error=${a}.trim.err"
    echo '#SBATCH --cpus-per-task=6'
    echo '#SBATCH --account=biodept'
    echo '#SBATCH -p common,biodept'
    echo "#SBATCH --chdir=$DIR/$a"
    echo '#SBATCH --mem=6G'
    echo ''
    echo "trimmomatic PE -threads 6 -phred33 -trimlog $a.trimlog filtered_forward.fastq.gz filtered_reverse.fastq.gz filtered.PE.R1.fq.gz filtered.U.R1.fq.gz filtered.PE.R2.fq.gz filtered.U.R2.fq.gz HEADCROP:2"
    echo "process_radtags --barcode-dist-1 0 --paired -1 filtered.PE.R1.fq.gz -2 filtered.PE.R2.fq.gz -b ../../$a.barcodes.txt -e pstI -i gzfastq -y gzfastq -r -c -q -E phred33 --inline_null"
  } > "$script_file"
done < "$file_name"
