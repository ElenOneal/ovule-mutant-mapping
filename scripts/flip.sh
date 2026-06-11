#!/bin/bash

# Create shell scripts to flip reads using flip2BeRad.py


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
  script_file="$a.flip.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.flip.out"
    echo "#SBATCH --error=${a}.flip.err"
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH --account=biodept'
    echo '#SBATCH -p common'
    echo "#SBATCH --chdir=$DIR/$a"
    echo '#SBATCH --mem=2G'
    echo ''
    echo "gunzip $a.rmdup.1.fq.gz"
    echo "gunzip $a.rmdup.2.fq.gz"
    echo "python ../../flip2BeRAD.py -c TGCAG -f $a.rmdup.1.fq -r $a.rmdup.2.fq -b ../../ddrad_barcode_seqs.txt -m 1 -o 2"
    echo "gzip *fq"
    echo "gzip *fastq"
  } > "$script_file"
done < "$file_name"
``
