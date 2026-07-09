#!/bin/bash
#
# This script builds a LepMap3 linkage mapping SLURM job script.
# It constructs a pedigree file, transposes it, then writes a .sh file
# that runs ParentCall2 -> Filtering2 -> SeparateChromosomes2 -> OrderMarkers2,
# then re-maps marker IDs to contig/position, splits by linkage group, and
# converts to phased genotype calls (see Linkage_map_assembly.Rmd for the
# underlying logic this automates).
#
# Usage: ./this_script.sh [-p PRE] [-o OUT] [-l LOD] [-t THETA] [-d dT] [-s MINORSCAFS] [-q PARTITION]
#   -p  Input VCF prefix (default: ovule_f2s.min60.max3)
#   -o  Output file prefix (default: ovule_f2s.min60.max3)
#   -l  LOD limit for SeparateChromosomes2 (default: 10)
#   -t  Theta limit for SeparateChromosomes2 (default: .20)
#   -d  Data tolerance for Filtering2 (default: 0.0001)
#   -s  Include scaffold contigs in split_lepmap.py output: True or False (default: False)
#   -q  SLURM partition (default: common)
# You can alter memory and CPU settings in the generated job script as needed.

# --- Default parameters (override with command-line flags) ---
PRE=ovule_f2s.min60.max3
OUT=ovule_f2s
lod='10'
theta='.20'
dT=0.0001
minorscafs='False'
partition='common'
DIR="$(pwd)"

# Resolve the directory this script lives in, so the generated job script can
# find split_lepmap.py regardless of how lepmap.sh itself was invoked.
script_dir="$(cd "$(dirname "$0")" && pwd)"

# --- Parse command-line options ---
usage() { echo "Usage: $0 [-p PRE] [-o OUT] [-l LOD] [-t THETA] [-d dT] [-s MINORSCAFS] [-q PARTITION]"; exit 1; }
while getopts "p:o:l:t:d:s:q:h" opt; do
case "$opt" in
p) PRE="$OPTARG" ;;
o) OUT="$OPTARG" ;;
l) lod="$OPTARG" ;;
t) theta="$OPTARG" ;;
d) dT="$OPTARG" ;;
s) minorscafs="$OPTARG" ;;
q) partition="$OPTARG" ;;
h) usage ;;
*) usage ;;
esac
done

# path to LepMap3 modules is ~/lepmap3. Adjust as needed for your system. Ensure you have Java installed and configured
LEPMAP_DIR=~/lepmap3

# --- Pedigree family/individual IDs ---
# Edit these to match your study organism and sample names
family_name="Fam"
GM="DHRO22"        # Maternal grandparent
GD="GMR2"          # Paternal grandparent
Mom="PseudoF1"     # F1 mother (offspring of GM x GD)
Dad="PseudoF1_D1"  # F1 father (offspring of GM x GD)

# --- Build pedigree file (ovule.ped.txt) ---
# LepMap3 expects a 6-column pedigree: family, individual, father, mother, sex, phenotype
# Sex coding: 1=male, 2=female, 0=unknown
rm ovule.ped.txt
rm ovule.ped.transpose.txt

# Header rows required by LepMap3
printf "CHR\tCHR\tCHR\tCHR\tCHR\tCHR\n" >> ovule.ped.txt
printf "POS\tPOS\tPOS\tPOS\tPOS\tPOS\n" >> ovule.ped.txt

# Grandparents (no parents listed, sex specified)
printf "$family_name\t$GM\t0\t0\t2\t0\n" >> ovule.ped.txt
printf "$family_name\t$GD\t0\t0\t1\t0\n" >> ovule.ped.txt

# F1 parents (offspring of grandparents)
printf "$family_name\t$Dad\t$GD\t$GM\t1\t0\n" >> ovule.ped.txt
printf "$family_name\t$Mom\t$GD\t$GM\t2\t0\n" >> ovule.ped.txt

# F2 individuals with known parentage but unknown sex
printf "$family_name\tDHRO22_D1\t$Dad\t$Mom\t0\t0\n" >> ovule.ped.txt
printf "$family_name\tGMR2_D1\t$Dad\t$Mom\t0\t0\n" >> ovule.ped.txt

# Add remaining F2 samples from file (one sample name per line)
while IFS=$'\n', read -r f ; do
printf "$family_name\t$f\t$Dad\t$Mom\t0\t0\n"  >> ovule.ped.txt
done < ovule_samples_nomissing.txt

# Transpose pedigree for LepMap3 input format (individuals as columns)
cat ovule.ped.txt | datamash transpose > ovule.ped.transpose.txt

# --- Write SLURM job script (lepmap<lod><theta>.<dT>.sh) ---
# The generated script will run the full LepMap3 pipeline as a SLURM job.
echo '#!/bin/bash' > lepmap$lod$theta.$dT.sh
echo '#' >> lepmap$lod$theta.$dT.sh

# SLURM directives
echo '#SBATCH --get-user-env' >> lepmap$lod$theta.$dT.sh
echo '#SBATCH --job-name='$lod,$theta   >> lepmap$lod$theta.$dT.sh
echo '#SBATCH --output='$lod$theta.out  >> lepmap$lod$theta.$dT.sh  # stdout and stderr go here
echo '#SBATCH --cpus-per-task=8' >> lepmap$lod$theta.$dT.sh
echo "#SBATCH -p ${partition}" >> lepmap$lod$theta.$dT.sh
echo "#SBATCH --chdir=${DIR}" >> lepmap$lod$theta.$dT.sh
echo '#SBATCH --mem=48G' >> lepmap$lod$theta.$dT.sh
echo -en '\n' >> lepmap$lod$theta.$dT.sh

# Step 1: ParentCall2 — call genotypes using pedigree and VCF
echo java -cp  ~/lepmap3/bin ParentCall2 data=ovule.ped.transpose.txt vcfFile=$PRE.thinned.vcf removeNonInformative=1 \> data.$OUT.call >> lepmap$lod$theta.$dT.sh
echo -en '\n' >> lepmap$lod$theta.$dT.sh

# Step 2: Filtering2 — filter markers by data tolerance and missing data rate
echo java -cp ~/lepmap3/bin Filtering2 data=data.$OUT.call dataTolerance=$dT removeNonInformative=1 outputHWE=0 missingLimit=0.30 \> data.$OUT.$dT.call >> lepmap$lod$theta.$dT.sh
echo -en '\n' >> lepmap$lod$theta.$dT.sh

# Step 3: SeparateChromosomes2 — cluster markers into linkage groups using LOD/theta thresholds
echo java -cp ~/lepmap3/bin/ SeparateChromosomes2 numThreads=8 grandparentPhase=1 data=data.$OUT.$dT.call lodLimit=$lod theta=$theta \> map.$lod$theta.$dT.$OUT.txt >> lepmap$lod$theta.$dT.sh
echo -en '\n' >> lepmap$lod$theta.$dT.sh

# Step 4: OrderMarkers2 — order markers within linkage groups and output phased data
echo java -cp ~/lepmap3/bin/ OrderMarkers2 numThreads=8 grandparentPhase=1 data=data.$OUT.$dT.call map=map.$lod$theta.$dT.$OUT.txt outputPhasedData=1 identicalLimit=0.01 \> order.$lod$theta.$dT.$OUT.txt >> lepmap$lod$theta.$dT.sh
echo -en '\n' >> lepmap$lod$theta.$dT.sh

# Steps 5-7: OrderMarkers2 identifies markers by an internal index, not by
# contig/position, so re-map its output using the .call file, split the
# result by linkage group, and convert to phased genotype calls.
# See Linkage_map_assembly.Rmd (LepMAP3 output processing) for the source of
# this logic.
cat >> lepmap$lod$theta.$dT.sh <<EOF

# Step 5: Re-map marker IDs to contig/position
cut -f1,2 data.$OUT.$dT.call | awk '(NR>=7)' > $OUT.$dT.snps.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=\$0}(NR!=FNR){if (\$1 in s) \$1=s[\$1];print}' $OUT.$dT.snps.txt order.$lod$theta.$dT.$OUT.txt > order.$lod$theta.$dT.$OUT.mapped

# Step 6: Split the mapped LepMap order file into per-linkage-group output files.
python ${script_dir}/split_lepmap.py order.$lod$theta.$dT.$OUT.mapped $OUT.$lod$theta $minorscafs

# Step 7: Convert the LepMap order output into genotype calls for downstream analysis.
awk -vfullData=1 -f ${LEPMAP_DIR}/map2genotypes.awk order.$lod$theta.$dT.$OUT.txt > $OUT.$lod$theta.phase.txt
EOF