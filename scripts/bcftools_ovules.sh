#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 listfiles bamfiles DIR genome_dir genome prefix"
    exit 1
fi

# Assign arguments to variables
list="$1"
bamfiles="$2"
DIR="$3"
genome_dir="$4"
genome="$5"
prefix="$6" # Added prefix as the sixth argument

# Check if the list file exists
if [ ! -f "$list" ]; then
    echo "Error: File '$list' not found."
    exit 2
fi

# Main loop

while IFS=$'\t' read -r a b; do
  script_file="Chr_${b}.bcftools.sh"
  {
    echo '#!/bin/bash'
    echo '#'
    echo '#$ -S /bin/bash' 
    echo "#SBATCH --get-user-env"
    echo "#SBATCH --job-name=Chr_${b}"
    echo "#SBATCH --output=Chr_${b}.${prefix}.out"
    echo "#SBATCH --error=Chr_${b}.${prefix}.err"
    echo '#SBATCH --cpus-per-task=6'
    echo '#SBATCH --account=biodept'
    echo '#SBATCH -p common,biodept'
    echo "#SBATCH --chdir=$DIR"
    echo '#SBATCH --mem=30G'
    echo ''
    echo "bcftools mpileup --threads 6 --redo-BAQ --min-MQ 29 --min-BQ 30 --per-sample-mF --annotate FORMAT/AD,FORMAT/DP,INFO/AD -f $genome_dir/$genome -b $bamfiles -I -r $a | bcftools call --multiallelic-caller | bcftools filter --SnpGap 3 -e 'QUAL<40 || INFO/RPBZ<-2 || INFO/RPBZ>2 || INFO/SCBZ<-2 || INFO/SCBZ > 2' -Oz -o Chr_${b}.${prefix}.filtered.vcf.gz"
    echo "tabix Chr_${b}.${prefix}.filtered.vcf.gz"
    echo "bcftools view -v snps Chr_${b}.${prefix}.filtered.vcf.gz -Oz -o Chr_${b}.${prefix}.snps_only.vcf.gz"
    echo "tabix Chr_${b}.${prefix}.snps_only.vcf.gz"
  } > "$script_file"
done < "$list"
``
