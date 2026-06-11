#!/bin/bash
set -euo pipefail

# Add read groups to pseudo-F1 parent BAMs in the 03_align directory.

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 sequences.txt DIR picard  partition"
    exit 1
fi

# Assign arguments to variables
file_name="$1"
DIR="$2"
picard="$3"
partition="$4"

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

# Check if the Picard JAR exists
if [ ! -f "$picard" ]; then
    echo "Error: Picard JAR '$picard' not found."
    exit 4
fi

mkdir -p "$DIR"

# sequences.txt file format:
#  <R1_fastq>    <R2_fastq>    <sample_id>    <parent_group>    <barcode>
# Example:
#  Mn_DHRO22_S8_L001_R1_001.fastq.gz    Mn_DHRO22_S8_L001_R2_001.fastq.gz    DHRO22_S8    DHRO22    ACTTGA
#  Mn_GMR2_S7_L001_R1_001.fastq.gz    Mn_GMR2_S7_L001_R2_001.fastq.gz    GMR2_S7    GMR2    CAGATC
#
# - first column is the forward read FASTQ filename
# - second column is the reverse read FASTQ filename
# - third column is the sample identifier used for output naming
# - fourth column is the parent or group name
# - fifth column is the barcode sequence
# - fields are tab-separated


# Main loop
while IFS=$'\t' read -r f c a b d; do
  script_file="${DIR}/${a}.rg.sh"
  {
    echo '#!/bin/bash'
    echo 'set -euo pipefail'
    echo '#'
    echo "#SBATCH --job-name=${a}"
    echo "#SBATCH --output=${a}.rg.out"
    echo "#SBATCH --error=${a}.rg.err"
    echo '#SBATCH --cpus-per-task=1'
    echo "#SBATCH -p ${partition}"
    echo "#SBATCH --chdir=${DIR}"
    echo '#SBATCH --mem=24G'
    echo ''
    echo "#make temp directory for Picard"
    echo "#mark duplicates with Picard, removing duplicates and indexing the filtered BAM"
    echo "mkdir -p \"${DIR}/temp/${a}\""
    echo "java -Xmx14g -jar \"${picard}\" MarkDuplicates INPUT=\"${a}.sort.bam\" OUTPUT=\"${a}.MD.bam\" M=\"${a}.metrics_file\" VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true"
    echo "samtools index \"${a}.MD.bam\""
    echo "#filter the marked BAM to retain only properly paired reads with mapping quality >= 29, excluding secondary and supplementary alignments, and indexing the filtered BAM"
    echo "samtools view -h \"${a}.MD.bam\" | awk 'BEGIN {OFS=\"\\t\"} {if(\$1 ~ /^@/) {print \$0; next;} if(\$7 == \"=\" || \$7 == \$3) {print \$0;}}' | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -f 2 -F 8 -q 29 -b -o \"${a}.filtered.bam\""
    echo "samtools index \"${a}.filtered.bam\""
    echo "#add read groups with Picard and index the final BAM"
    echo "java -Xmx14g -jar \"${picard}\" FixMateInformation INPUT=\"${a}.filtered.bam\" OUTPUT=\"${a}.FM.bam\" SORT_ORDER=coordinate TMP_DIR=\"${DIR}/temp/${a}\" VALIDATION_STRINGENCY=LENIENT"
    echo "samtools index \"${a}.FM.bam\""
    echo "java -Xmx14g -jar \"${picard}\" AddOrReplaceReadGroups RGLB=${a} RGPL=illumina RGPU=run RGSM=${b} I=\"${a}.FM.bam\" O=\"${a}.RG.bam\" SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT TMP_DIR=\"${DIR}/temp/${a}\""
    echo "samtools index \"${a}.RG.bam\""
  } > "$script_file"
  chmod +x "$script_file"
done < "$file_name"

