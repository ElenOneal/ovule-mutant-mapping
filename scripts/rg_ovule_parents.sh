#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 sequences.txt DIR picard genome_dir genome"
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
while IFS=$'\t' read -r f c a b d; do
  script_file="$a.filter.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=$a"
    echo "#SBATCH --output=${a}.filter.out"
    echo "#SBATCH --error=${a}.filter.err"
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH --account=biodept'
    echo '#SBATCH -p common,biodept'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=24G'
    echo ''
    echo "mkdir $DIR/temp/$a"
    echo "#java -Xmx14g -jar $picard MarkDuplicates INPUT=$a.sort.bam OUTPUT=$a.MD.bam M=$a.metrics_file VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true"
    echo "#samtools index $a.MD.bam" 
    echo "#samtools view -h $a.MD.bam | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o $a.filtered.bam"
    echo "#samtools index $a.filtered.bam"
    echo "#java -Xmx14g -jar $picard FixMateInformation INPUT=$a.filtered.bam OUTPUT=$a.FM.bam SORT_ORDER=coordinate TMP_DIR=$DIR/temp/$a VALIDATION_STRINGENCY=LENIENT"
    echo "#samtools index $a.FM.bam"
    echo "java -Xmx14g -jar $picard AddOrReplaceReadGroups RGLB=$a RGPL=illumina RGPU=run RGSM=$b I=$a.FM.bam O=$a.RG.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT TMP_DIR=$DIR/temp/$a"
    echo "samtools index $a.RG.bam"
  } > "$script_file"
done < "$file_name"
``
