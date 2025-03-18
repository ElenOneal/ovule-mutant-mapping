#!/bin/bash

# Create shell scripts to remove duplicate reads using rmdup_molbarcodes.py file


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
  script_file="$a.rmdup.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.rmdup.out"
    echo "#SBATCH --error=${a}.rmdup.err"
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH -p common'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=5G'
    echo ''
    echo "python rmdup_molbarcodes.py -p $a -s fq.gz"
  } > "$script_file"
done < "$file_name"
``
