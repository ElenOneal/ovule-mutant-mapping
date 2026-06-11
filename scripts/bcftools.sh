#!/bin/bash
set -euo pipefail

# Create shell scripts to call variants using bcftools mpileup and call, and filter variants using bcftools filter

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 listfiles bamfiles DIR genome_dir genome prefix partition"
    exit 1
fi

# Assign arguments to variables
list="$1"
bamfiles="$2"
DIR="$3"
genome_dir="$4"
genome="$5"
prefix="$6" # Added prefix as the sixth argument
partition="$7" # Added partition as the seventh argument

# Check if the list file exists
if [ ! -f "$list" ]; then
    echo "Error: File '$list' not found."
    exit 2
fi

# Check if the BAM list exists
if [ ! -f "$bamfiles" ]; then
    echo "Error: BAM list file '$bamfiles' not found."
    exit 3
fi

# Check if the genome directory exists
if [ ! -d "${genome_dir}" ]; then
    echo "Error: Genome directory '${genome_dir}' not found."
    exit 4
fi

# Check if the reference genome exists
if [ ! -f "${genome_dir}/${genome}" ]; then
    echo "Error: Reference genome '${genome_dir}/${genome}' not found."
    exit 5
fi

# List file format:
#  <region>    <chromosome_id>
# Example:
#  Chr_01:1-13923366    01
#  Chr_02:1-19757465    02
#
# - first column is the region passed to bcftools mpileup with -r
# - second column is the chromosome suffix used in output file names
# - fields are tab-separated

# Ensure the output directory exists
mkdir -p "$DIR"

# Main loop

while IFS=$'\t' read -r a b; do
  script_file="${DIR}/Chr_${b}.${prefix}.bcftools.sh"
  cat > "$script_file" <<EOF
#!/bin/bash
set -euo pipefail
#
#SBATCH --get-user-env
#SBATCH --job-name=Chr_${b}.${prefix}
#SBATCH --output=Chr_${b}.${prefix}.out
#SBATCH --error=Chr_${b}.${prefix}.err
#SBATCH --cpus-per-task=6
#SBATCH -p ${partition}
#SBATCH --chdir=${DIR}
#SBATCH --mem=30G

# call bcftools then filter with the following: 1) minimum mapping quality = 29, 2) minimum base quality = 30, 3) minimum distance between SNPs = 3, 4) strand bias filter (RPBZ and SCBZ) with a threshold of -2 to 2
bcftools mpileup --threads 6 --redo-BAQ --min-MQ 29 --min-BQ 30 --per-sample-mF --annotate FORMAT/AD,FORMAT/DP,INFO/AD -f "${genome_dir}/${genome}" -b "${bamfiles}" -I -r "${a}" | bcftools call --multiallelic-caller | bcftools filter --SnpGap 3 -e 'QUAL<40 || INFO/RPBZ<-2 || INFO/RPBZ>2 || INFO/SCBZ<-2 || INFO/SCBZ > 2' -Oz -o "Chr_${b}.${prefix}.filtered.vcf.gz"
tabix "Chr_${b}.${prefix}.filtered.vcf.gz"
# keep snps only, excluding indels and other variant types, and index the SNP-only VCF
bcftools view -v snps "Chr_${b}.${prefix}.filtered.vcf.gz" -Oz -o "Chr_${b}.${prefix}.snps_only.vcf.gz"
tabix "Chr_${b}.${prefix}.snps_only.vcf.gz"
EOF
  chmod +x "$script_file"
done < "$list"

