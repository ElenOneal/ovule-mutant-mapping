#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 sequences.txt DIR clean_dir execute_dir"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
clean_dir="$3"
execute_dir="$4"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r f c a b; do
  script_file="$a.trim.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.out"
    echo "#SBATCH --error=${a}.err"
    echo '#SBATCH --cpus-per-task=6'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=48G'
    echo ''
    echo "trimmomatic PE -threads 6 -phred33 -trimlog $a.trimlog -quiet -validatePairs $f $c $clean_dir/$a.PE.R1.fq.gz $clean_dir/$a.U.R1.fq.gz $clean_dir/$a.PE.R2.fq.gz $clean_dir/$a.U.R2.fq.gz ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 AVGQUAL:25"
  } > "$script_file"
done < "$file_name"
``
