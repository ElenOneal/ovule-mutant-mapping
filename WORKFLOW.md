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
conda activate f1mapping
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
   conda activate f1mapping
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
./rmdup.sh data/ovule_library_filenames.txt 01_rmdup common
for file in rmdup_*.sh; do sbatch $file; done
```
Output: Deduplicated reads with `*rmdup*.gz`

### 1.3 Flip Reads (BeRAD protocol)
```bash
./flip.sh data/ovule_library_filenames.txt 02_flipreads common
for file in *.flip.sh; do sbatch $file; done
```
Output: `filtered_forward.fastq` and `filtered_reverse.fastq` per sample

### 1.4 Trim and Extract Barcodes
```bash
./process_radtags.sh data/ovule_library_filenames.txt 02_flipreads common
for file in *.radtags.sh; do sbatch $file; done
```
Move trimmed reads to alignment directory:
```bash
mkdir 03_align
mv 02_flipreads/*/*.1.fq.gz 03_align/
mv 02_flipreads/*/*.2.fq.gz 03_align/
```

---

## Step 2: Alignment and Genotyping

### 2.1 Trim and Align F2 Samples
```bash
./align_f2s.sh data/ovule_barcode_key.txt 03_align genome_dir Mguttatus.IM767.v2.fa common
./rg_f2s.sh data/ovule_barcode_key.txt 03_align picard_location genome_dir Mguttatus.IM767.v2.fa common
```

### 2.2 Process Parent Reads
```bash
./align_parents.sh data/ovule_genome_sequence_files.txt 03_align 00_parent_reads genome_dir Mguttatus.IM767.v2.fa common
./rg_parents.sh data/ovule_genome_sequence_files.txt 03_align/ picard_location common
```

### 2.2b Create F1 Read Groups
```bash
./f1_readgroups.sh
```

### 2.3 Call SNPs with bcftools
```bash
./bcftools.sh im767.v2.list bamlist.txt 04_bcftools/ genome_dir Mguttatus.IM767.v2.fa ovule_f2s common
```
This produces one SNP-only VCF per chromosome, e.g. `04_bcftools/Chr_01.ovule_f2s.snps_only.vcf.gz`.

### 2.4 Concatenate per-chromosome VCFs
```bash
./concatenate_vcfs.sh im767.v2.list ovule_f2s 04_bcftools/ common
sbatch 04_bcftools/ovule_f2s.concat.sh
```
Reuses the same `im767.v2.list` and `ovule_f2s` prefix passed to `bcftools.sh` to merge the per-chromosome VCFs into one genome-wide file and unpack a plain-text copy. `filter_vcf_for_lepmap3_ovule_cross.py` reads `--invcf` with plain `open()` (not `gzip.open()`), so the plain-text copy is required.

Output: `04_bcftools/ovule_f2s.snps_only.vcf.gz` (+ `.tbi` index) and `04_bcftools/ovule_f2s.snps_only.vcf` (plain text, used in Step 3.1)

---

## Step 3: Filter and Prepare for Linkage Mapping

### 3.1 Initial VCF Filtering
```bash
mkdir -p 05_filter
cp 04_bcftools/ovule_f2s.snps_only.vcf 05_filter/
cd 05_filter/
module add VCFtools
vcftools --vcf ovule_f2s.snps_only.vcf --missing-indv
python ../scripts/filter_vcf_for_lepmap3_ovule_cross.py \
  --invcf ovule_f2s.snps_only.vcf \
  --ParentList ../data/ovule_parents.txt \
  --SampleList ../data/ovule_f2_samples.txt \
  --out ovule_f2s.min60 \
  --minMapQ 20 \
  --filterInbredParents True \
  --filterforhets True \
  --minparentdepth True \
  --maxparentdepth True \
  --missingList ovule_f2s.min60.max3.frac_missing_data.txt \
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

## Step 4: Linkage Map Assembly with LepMAP3

### 4.1 Run LepMAP3 Pipeline
```bash
cd 05_filter/
./lepmap.sh -o ovule_f2s -l 10 -t .20 -s False
```
This writes a SLURM job script (`lepmap10.20.0.0001.sh`) that runs, in sequence, with LOD=10, theta=0.20, dT=0.0001:
- `ParentCall2` – Phase parental genotypes
- `Filtering2` – Filter informative markers by data tolerance threshold
- `SeparateChromosomes2` – Group markers into linkage groups
- `OrderMarkers2` – Order markers within groups
- Re-map marker IDs to contig/position (`OrderMarkers2` identifies markers by an internal index, not contig/position; see `Linkage_map_assembly.Rmd` for the underlying logic) and split into per-linkage-group output via `split_lepmap.py`
- Convert to phased genotype calls via `map2genotypes.awk`

Submit it, and wait for it to finish before continuing:
```bash
sbatch lepmap10.20.0.0001.sh
```

Output:
- `order.10.20.0.0001.ovule_f2s.txt` — raw `OrderMarkers2` output (markers by internal index)
- `order.10.20.0.0001.ovule_f2s.mapped` — re-mapped to contig/position
- `ovule_f2s.10.20_mappingorder.txt` — per-linkage-group mapping order (mapped markers with cM positions)
- `ovule_f2s.10.20.phase.txt` — phased genotypes

---

## Step 5: Data Curation and QTL Preparation

### 5.1 Replace Missing Data Calls
Run in R (within `Linkage_map_assembly.Rmd`, `replace_missing_phases` chunk):
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
**Not included in this public repository pending publication** (see Data Access in README.md). This step subsamples phased genotypes and phenotypes into rqtl format, producing:
- `data/ovule_f2s_quantgen.csv` – Subsampled genotypes in rqtl format (1 marker per 25 kb)
- `data/ovule_f2s_rqtlphen.csv` – Phenotypes (WT=1, MUT/MOSAIC=0)

### 6.2 Run QTL Mapping (Optional)
Open and render `Rqtl.Rmd` for full EM-based QTL mapping with permutation tests:
```bash
conda activate linkage-mapping_rmd
Rscript -e "rmarkdown::render('Rqtl.Rmd')"
```

Outputs:
- QTL scan plots and LOD scores
- Permutation-derived significance thresholds
- Two-dimensional scan results
- Model comparison and QTL refinement results

---

## Key Data Files

| File | Description |
|------|-------------|
| `data/ovule_samples.txt` | List of F2 sample IDs |
| `data/ovule_barcode_key.txt` | Barcode-to-sample mapping |
| `data/FinalPhenotypeCalls.txt` | Phenotypic calls (WT/MUT/MOSAIC) |
| `data/ovule_sites.bed` | Final SNP positions used in QTL analysis |
| `ovule.qtl.txt` | QTL mapping results summary |
| `ovule_sub_genotypes_missing_lowhet.txt` | Phased genotype data |

---

## Intermediate Outputs

- `01_rawreads/` – Deduplicated fastq files
- `02_flipreads/` – Flipped, protocol-corrected reads
- `03_align/` – Aligned BAM files and read groups
- `04_bcftools/` – Per-chromosome and concatenated genome-wide VCFs (SNPs only)
- `05_filter/` – Filtered VCF and mapping order

---

## Notes

- **Conda environments:** Use `linkage-mapping` for main pipeline, `py2` for Python 2 scripts, `linkage-mapping_rmd` for report rendering
- **LepMAP3 parameters:** LOD 10, theta 0.20, data tolerance 0.0001 produced best map (see `Linkage_map_assembly.Rmd` for justification)
- **Marker filtering:** Biallelic sites only before LepMAP3; genotype calls limited to 3 alleles max; markers thinned to 150 bp for initial scan
- **QTL marker subset:** Thinned to 25 kb spacing to reduce computation; phased genotypes tested for missing data calls
- **F2 phenotyping:** Ovule mutant phenotype segregates at ~25–30% consistent with single-locus Mendelian inheritance
