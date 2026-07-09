#!/bin/bash
set -euo pipefail

# Concatenate the per-chromosome SNP-only VCFs produced by bcftools.sh into a
# single genome-wide VCF, ready for filter_vcf_for_lepmap3_ovule_cross.py.
#
# Reuses the same "listfiles" file passed to bcftools.sh, so the chromosome
# list only has to be maintained in one place.
#
# List file format (same as bcftools.sh expects):
#  <region>    <chromosome_id>
# Example:
#  Chr_01:1-13923366    01
#  Chr_02:1-19757465    02
# - fields are tab-separated

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 listfiles prefix bcftools_dir partition"
    exit 1
fi

# Assign arguments to variables
list="$1"          # same listfiles file used by bcftools.sh
prefix="$2"        # same prefix used by bcftools.sh (6th arg to bcftools.sh)
bcftools_dir="$3"  # directory containing bcftools.sh's output (its DIR arg)
partition="$4"

# Check if the list file exists
if [ ! -f "$list" ]; then
    echo "Error: File '$list' not found."
    exit 2
fi

# Check if the bcftools output directory exists
if [ ! -d "$bcftools_dir" ]; then
    echo "Error: Directory '$bcftools_dir' not found."
    exit 3
fi

# Build the list of per-chromosome SNP-only VCFs to concatenate, in listfile
# order. Filenames must match bcftools.sh's naming: Chr_<id>.<prefix>.snps_only.vcf.gz
vcf_list_file="${bcftools_dir}/${prefix}.snps_only.vcffiles.txt"
> "$vcf_list_file"
while IFS=$'\t' read -r region chrom_id; do
  [ -z "$chrom_id" ] && continue
  vcf="${bcftools_dir}/Chr_${chrom_id}.${prefix}.snps_only.vcf.gz"
  if [ ! -f "$vcf" ]; then
    echo "Error: Expected VCF '$vcf' not found. Has bcftools.sh finished for chromosome ${chrom_id}?"
    exit 4
  fi
  echo "$vcf" >> "$vcf_list_file"
done < "$list"

# Write the SLURM job script that concatenates and indexes the merged VCF.
script_file="${bcftools_dir}/${prefix}.concat.sh"
cat > "$script_file" <<EOF
#!/bin/bash
set -euo pipefail
#
#SBATCH --get-user-env
#SBATCH --job-name=${prefix}.concat
#SBATCH --output=${prefix}.concat.out
#SBATCH --error=${prefix}.concat.err
#SBATCH --cpus-per-task=1
#SBATCH -p ${partition}
#SBATCH --chdir=${bcftools_dir}
#SBATCH --mem=4G

bcftools concat -f "${prefix}.snps_only.vcffiles.txt" -Oz -o "${prefix}.snps_only.vcf.gz"
tabix "${prefix}.snps_only.vcf.gz"

# filter_vcf_for_lepmap3_ovule_cross.py reads --invcf with plain open(), not
# gzip.open(), so it needs an uncompressed VCF. Keep the .gz + .tbi above for
# indexing/region queries, and unpack a plain copy here for the filtering step.
gunzip -k "${prefix}.snps_only.vcf.gz"
EOF
chmod +x "$script_file"

echo "Wrote ${script_file}"
echo "Submit with: sbatch ${script_file}"
