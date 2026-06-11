# Ovule Mutant Mapping Workflow

This document outlines the complete analysis pipeline for QTL mapping in the ovule mutant project. For detailed code and analysis, see `Linkage_map_assembly.Rmd` and `Rqtl.Rmd`.

## Prerequisites

Create conda environments before starting:
```bash
conda env create -f environment/linkage-mapping_environment.yml
   conda env create -f environment/py2_environment.yml
   conda env create -f environment/linkage-mapping_rmd.yml
```

Activate the main environment:
```bash
conda activate linakge-mapping
```

- **Git:** Ensure you have Git installed. [Download Git](https://git-scm.com/downloads)
- **Python:** Version 3.8.12 is recommended.
- **R and RStudio:** For R Markdown scripts.
- **Rqtl** for QTL mapping

### Installation

1. **Clone the Repository using SSH:**
   ```bash
   git clone git@github.com:ElenOneal/ovule-mutant-mapping.git
   ```

2. **Create Conda Environments:**
   ```bash
   conda env create -f environment/linkage-mapping_environment.yml
   conda env create -f environment/py2_environment.yml
   conda env create -f environment/linkage-mapping_rmd.yml
   ```

3. **Activate the linkage mapping environment:**
   ```bash
   conda activate linakge-mapping
   ```


---

## Step 1: Data Preprocessing

### 1.1 Prepare Sequencing Files
- Copy raw fastq files from sequencing run to working directory
- Organize files by library plate (JW77, JW78, etc.)

### 1.2 Remove PCR Duplicates
Activate `py2` environment for Python 2.7 compatibility:
```bash
conda activate py2
bash scripts/rmdup.sh data/ovule_library_filenames.txt 01_rawreads/
for file in *.rmdup.sh; do sbatch $file; done
```
Output: Deduplicated reads with `*rmdup*.gz`

### 1.3 Flip Reads (BeRAD protocol)
```bash
conda activate py2
bash scripts/flip.sh data/ovule_library_filenames.txt 02_flipreads/
for file in *.flip.sh; do sbatch $file; done
```
Output: `filtered_forward.fastq` and `filtered_reverse.fastq` per sample

### 1.4 Trim and Extract Barcodes
```bash
bash scripts/trim_ddrad.sh data/ovule_library_filenames.txt 02_flipreads/
for file in *.trim.sh; do sbatch $file; done
```
Move trimmed reads to alignment directory:
```bash
mkdir 03_align
mv 02_flipreads/*/*.1.fq.gz 03_align/
mv 02_flipreads/*/*.2.fq.gz 03_align/
rm 03_align/*rem*.fq.gz
```

---

## Step 2: Alignment and Genotyping

### 2.1 Trim and Align F2 Samples
```bash
bash scripts/trim_align_ovules.sh data/ovule_barcode_key.txt 03_align/ /path/to/genome/ Mguttatus.IM767.v2.fa
bash scripts/rg_ovules.sh data/ovule_barcode_key.txt 03_align/ /path/to/picard.jar /path/to/genome/ Mguttatus.IM767.v2.fa
```

### 2.2 Process Parent Reads
```bash
bash scripts/trim_parents.sh data/ovule_genome_sequence_files.txt /parent_raw_reads/ /parent_clean_reads/ 03_align/
bash scripts/align_parents.sh data/ovule_genome_sequence_files.txt 03_align/ /parent_clean_reads/ /path/to/genome/ Mguttatus.IM767.v2.fa
bash scripts/rg_ovule_parents.sh data/ovule_genome_sequence_files.txt 03_align/ /path/to/picard.jar /path/to/genome/ Mguttatus.IM767.v2.fa
```

### 2.3 Call SNPs with bcftools
```bash
bash scripts/bcftools_ovules.sh im767.v2.list bamlist.txt 04_bcftools/ /path/to/genome/ Mguttatus.IM767.v2.fa ovule_f2s
cd 04_bcftools/
bcftools concat -f vcffiles_snpsonly.txt -Ov -o ovule_f2s.snps_only.vcf
```

---

## Step 3: Filter and Prepare for Linkage Mapping

### 3.1 Initial VCF Filtering
```bash
cd 05_filter/
module add VCFtools
vcftools --vcf ovule_f2s.snps_only.vcf --missing-indv
python scripts/filter_vcf_for_lepmap3_ovule_cross.py \
  --invcf ovule_f2s.snps_only.vcf \
  --ParentList ovule_parents.txt \
  --SampleList ovule_samples.txt \
  --out ovule_f2s.min60 \
  --minMapQ 20 \
  --filterInbredParents True \
  --filterforhets True \
  --minparentdepth --maxparentdepth \
  --missingList ovule_f2s.frac_missing_data.txt \
  --MaxMissingData 1 \
  --filterMissingParents True \
  --FractionF2 60 \
  --minF2depth 4
```

### 3.2 Remove High-Genotype Sites
```bash
grep -v "0/2" ovule_f2s.min60.filtered.vcf | grep -v "2/3" > ovule_f2s.min60.max3.vcf
vcftools --vcf ovule_f2s.min60.max3.vcf --thin 150 --recode --recode-INFO-all --stdout > ovule_f2s.min60.max3.thinned.vcf
```

### 3.3 Create Pedigree File
Build `ovule.ped.txt` with grandparents (DHRO22, GMR2) and F1s as cross parents.

---

## Step 4: Linkage Map Assembly with LepmAP3

### 4.1 Run LepmAP3 Pipeline
```bash
bash scripts/lepmap10.20.0.0001.sh
```
This runs (in sequence):
- `ParentCall2` – Phase parental genotypes
- `Filtering2` – Filter informative markers
- `SeparateChromosomes2` – Group markers into linkage groups (LOD 10, theta 0.20)
- `OrderMarkers2` – Order markers within groups

Output: `order.10.20.0.0001.id.ovule_f2s.txt` (mapped markers with cM positions)

### 4.2 Extract Phased Genotypes
```bash
awk -vfullData=1 -f map2genotypes.awk order.10.20.0.0001.id.ovule_f2s.txt > ovule_f2s.10.20.phase1.txt
```

---

## Step 5: Data Curation and QTL Preparation

### 5.1 Replace Missing Data Calls
Run in R (within `scripts/Ovules.Rmd`):
- Replace LepMAP3 phased missing calls with actual VCF calls
- Output: `ovule_f2s_phased_w_missing.txt`

### 5.2 Filter High-Heterozygosity Samples
Identify and remove F2s with excessive heterozygosity (likely contamination):
- Use `het_distribution.R` to visualize
- Remove samples with >80% heterozygosity

### 5.3 Create Marker Subset
Thin markers to 25 kb spacing for QTL analysis:
```bash
# In R
ovule.qtl <- ovule.geno %>%
  group_by(lg) %>%
  arrange(position) %>%
  filter(c(TRUE, diff(position) >= 25000))
```

---

## Step 6: QTL Analysis with rqtl

### 6.1 Prepare rqtl Input Files
Create from phased genotypes and phenotypes:
- `ovule_f2s_quantgen.csv` – Genotypes in rqtl format
- `ovule_f2s_rqtlphen.csv` – Phenotypes (WT=1, MUT=0)

### 6.2 Run QTL Mapping
Open and render `scripts/Ovules.Rmd`:
```bash
cd scripts/
Rscript -e "rmarkdown::render('Ovules.Rmd')"
```

Outputs:
- `ovules_rqtl_full.html` – QTL scan results and plots
- LOD scores and QTL positions
- Linkage group plots

---

## Key Data Files

| File | Description |
|------|-------------|
| `data/ovule_samples.txt` | List of F2 sample IDs |
| `data/ovule_barcode_key.txt` | Barcode-to-sample mapping |
| `data/FinalPhenotypeCalls.txt` | Phenotypic calls (WT/MUT/MOSAIC) |
| `data/ovule_sites.bed` | Final SNP positions used in QTL analysis |
| `ovule.qtl.txt` | QTL mapping results summary |

---

## Intermediate Outputs

- `01_rawreads/` – Deduplicated fastq files
- `02_flipreads/` – Flipped, protocol-corrected reads
- `03_align/` – Aligned BAM files and read groups
- `04_bcftools/` – VCF files (SNPs only)
- `05_filter/` – Filtered VCF and mapping order

---

## Notes

- **Conda environments:** Use `linkage-mapping` for main pipeline, `py2` for deaggregation only
- **LepmAP3 parameters:** LOD 10, theta 0.20 produced best map (see Ovules.Rmd for justification)
- **Marker filtering:** Thinned to 150 bp spacing before LepmAP3, then 25 kb for QTL to reduce computation
- **F2 phenotyping:** Medusa mutants segregate at ~25–30% as expected for single locus
