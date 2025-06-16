#!/bin/bash

# Create shell scripts to trim reads and align to genome


# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 sample_names DIR picard genome_dir genome"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
picard="$3"
genome_dir="$4"
genome="$5"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Main loop
while IFS=$'\t' read -r a b c d; do
  script_file="$b.rg.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$b"
    echo "#SBATCH --output=${b}.out"
    echo "#SBATCH --error=${b}.err"
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH --mail-type=ALL'
    echo '#SBATCH --mail-user=eo22@duke.edu'
    echo '#SBATCH -p common,scavenger'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=8G'
    echo ''
    echo "mkdir $DIR/temp/$b"
    echo "samtools index $b.sort.bam"
    echo "samtools view -h $b.sort.bam | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o $b.filtered.bam"
    echo "samtools index $b.filtered.bam"
    echo "samtools flagstat $b.filtered.bam"
    echo "java -Xmx6g -jar $picard AddOrReplaceReadGroups RGLB=$b RGPL=ILLUMINA RGPU=$b RGSM=$b I=$b.filtered.bam O=$b.RG.bam VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=$DIR/temp/$b"
    echo "samtools index $b.RG.bam"
    echo "rm -r $DIR/temp/$b "
  } > "$script_file"
done < "$file_name"
``

