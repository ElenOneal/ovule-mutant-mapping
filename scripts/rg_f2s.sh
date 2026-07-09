#!/bin/bash
set -euo pipefail

# Create shell scripts to trim reads and align to genome

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 sample_names DIR picard genome_dir genome partition"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
picard="$3"
genome_dir="$4"
genome="$5"
partition="$6"

# Check if the input file exists
if [ ! -f "$file_name" ]; then
    echo "Error: File '$file_name' not found."
    exit 2
fi

# Check if the output directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' not found."
    exit 3
fi

# sample_names file format:
#  <barcode>    <assigned_F2_name>    <external_F2_name> <cross?
# Example:
#  GGAAGAA    F2.01    1A crossname
#  CAGAGAA    F2.02    6E crossname
#
# - first column is the barcode sequence
# - second column is the F2 name assigned in this project
# - third column is the F2 name used by the other lab
# - fourth column is the cross name
# - fields are tab-separated

# Main loop
while IFS=$'\t' read -r a b c rest; do
  script_file="${b}.rg.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo "#SBATCH --job-name=${b}"
    echo "#SBATCH --output=${b}.out"
    echo "#SBATCH --error=${b}.err"
    echo '#SBATCH --cpus-per-task=1'
    echo "#SBATCH -p ${partition}"
    echo "#SBATCH --chdir=${DIR}"
    echo '#SBATCH --mem=8G'
    echo ''
    echo "mkdir -p \"${DIR}/temp/${b}\""
    echo "samtools index \"${b}.sort.bam\""
    echo "samtools view -h \"${b}.sort.bam\" | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o \"${b}.filtered.bam\""
    echo "samtools index \"${b}.filtered.bam\""
    echo "samtools flagstat \"${b}.filtered.bam\""
    echo "java -Xmx6g -jar \"${picard}\" AddOrReplaceReadGroups RGLB=${b} RGPL=ILLUMINA RGPU=${b} RGSM=${b} I=\"${b}.filtered.bam\" O=\"${b}.RG.bam\" VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=\"${DIR}/temp/${b}\""
    echo "samtools index \"${b}.RG.bam\""
    echo "rm -rf \"${DIR}/temp/${b}\""
  } > "$script_file"
done < "$file_name"


