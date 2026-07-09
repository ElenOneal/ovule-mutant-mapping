# ovule-mutant-mapping

A repository for the analysis scripts and data used in the Mimulus ovule mutant mapping project. This project includes shell scripts, Python scripts, R scripts, and Rmarkdown files for data analysis and mapping.

The Medusa Mutant Disrupts Ovule Integument and Female Gametophyte Development in Mimulus.

Authors: Miguel Flores-Vergara, Albert Tucci, Elen Oneal, Alessandro Ruiu, other authors,  Robert G. Franks* (final authorship to be determined)



## Table of Contents

- [Project Overview](#project-overview)
- [Getting Started](#getting-started)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Project Overview

This mapping project was used to discover the locus underlying an ovule mutant phenotype first observed in a Mimulus nudatus accesion (GMR2). This line was at least 6 generations inbred when some individuals were observed with an unusual ovule phenotype. Instead of exhibiting normal development, the ovules took on a "Medusa"-like appearance.  In Arabidopsis, BEL1 is involved in the establishment of ovule identity. Early qPCR work of the expression of the BEL1 homolog in M. nudatus were inconclusive.

## Rationale
We did QTL mapping in an F2 mapping population to discover the loci underlying the ovule mutant in M. nudatus GMR2. The M. nudatus inbred line DHRO22 exhibits the wildtype phenotype. An F1 (DHRO22 x GMR2) was created, then selfed to create F2s. Phenotyping of F2s revealed that the ovule mutant appears to segregate at a 25-35% rate, suggesting a single locus.


A conda virtual environment (“linkage-mapping”) was employed for NGS bioinformatics and linkage map assembly. A description of the conda environment linkage-mapping, including package versions, can be found in linkage-mapping_environment.yml.  This environment used Python 3.8.12. A separate conda environment ("py2") was employed when deaggregating the ddRAD mapping data for the F2s used in this project using scripts written by Thom Nelson. "Py2" used Python 2.7.18. Finally, the conda environment "linkage-mapping_rmd" was used for data analysis in R and requires R version 4.0.5 or above.

An Rmarkdown file (Linkage_map_assembly.Rmd) was created to document the ddRAD linkage map assembly pipeline.

This project includes a series of scripts written in Python and shell, along with R Markdown to document flow, perform statistical analysis and generate reproducible figures.

##  Data Access

The full dataset and figures will be available upon publication.  For now, this repository provides code and illustrative examples only. Phenotype data and rQTL input files will be made available upon publication. Contact elen.oneal@duke.edu for reviewer access.

## Getting Started

These instructions will help you set up the project on your local machine.

### Prerequisites

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

## Repository Structure

- `scripts/`: scripts for data processing, QTL mapping, visualization.
- `data/`: Raw and processed data files.
- `figs/`: Output figures and visualizations.
- `environment/`: Conda environment YAML files for reproducibility.
- `LICENSE`: License for code reuse (MIT).

### Key Files

- `Linkage_map_assembly.Rmd`: R Markdown documenting NGS pipeline.
- `Rqtl.Rmd`: R markdown documenting rQTL analysis.
- `WORKFLOW.md`: Step-by-step pipeline notes for running the full analysis.
- `scripts/ovules_rqtl_full.Rmd`: Expanded in-repo rQTL analysis notebook.

## Usage

### Main Analysis Workflow
    ```bash
    bash scripts/rmdup.sh data/ovule_library_filenames.txt 01_rmdup common
    bash scripts/flip.sh data/ovule_library_filenames.txt 02_flipreads common
    bash scripts/process_radtags.sh data/ovule_library_filenames.txt 02_flipreads common
    bash scripts/align_f2s.sh data/ovule_barcode_key.txt 03_align /path/to/genome_dir Mguttatus.IM767.v2.fa common
    bash scripts/rg_f2s.sh data/ovule_barcode_key.txt 03_align /path/to/picard.jar /path/to/genome_dir Mguttatus.IM767.v2.fa common
    ```

2. **Process parent reads:**
    ```bash
    bash scripts/align_parents.sh data/ovule_genome_sequence_files.txt 03_align /path/to/parent_fastqs /path/to/genome_dir Mguttatus.IM767.v2.fa common
    bash scripts/rg_parents.sh data/ovule_genome_sequence_files.txt 03_align /path/to/picard.jar common
    ```

3. **Variant calling:**
    ```bash
    bash scripts/bcftools.sh data/im767.v2.list data/bamlist.txt 04_bcftools /path/to/genome_dir Mguttatus.IM767.v2.fa ovule_f2s common
    ```
    This produces one SNP-only VCF per chromosome (e.g. `04_bcftools/Chr_01.ovule_f2s.snps_only.vcf.gz`).

4. **Concatenate per-chromosome VCFs:**
    ```bash
    bash scripts/concatenate_vcfs.sh data/im767.v2.list ovule_f2s 04_bcftools common
    sbatch 04_bcftools/ovule_f2s.concat.sh
    ```
    Reuses the same `listfiles` file passed to `bcftools.sh` (must have the same `prefix`, here `ovule_f2s`) to merge the per-chromosome `snps_only.vcf.gz` files into one genome-wide VCF, then unpacks a plain-text copy: `04_bcftools/ovule_f2s.snps_only.vcf`. The uncompressed copy is required because `filter_vcf_for_lepmap3_ovule_cross.py` reads `--invcf` with plain `open()`, not `gzip.open()`.

5. **VCF filtering and LepMap3 preparation:**
    ```bash
    python scripts/filter_vcf_for_lepmap3_ovule_cross.py \
       --invcf 04_bcftools/ovule_f2s.snps_only.vcf \
       --ParentList data/ovule_parents.txt \
       --SampleList data/ovule_f2_samples.txt \
       --out ovule_f2s.min60 \
       --minMapQ 20 \
       --filterMissingParents True \
       --filterInbredParents True \
       --filterforhets True \
       --minparentdepth True \
       --maxparentdepth True \
       --FractionF2 60 \
       --minF2depth 4 \
       --missingList data/ovule_f2s.min60.max3.frac_missing_data.txt \
       --MaxMissingData 1
    ```

6. **Linkage map assembly:**
    ```bash
    bash scripts/lepmap.sh -p ovule_f2s.min60.max3 -o ovule_f2s -l 10 -t .20 -s False -q common
    sbatch lepmap10.20.0.0001.sh
    ```
    This single job runs the full LepMap3 pipeline (`ParentCall2` → `Filtering2` → `SeparateChromosomes2` → `OrderMarkers2`), then re-maps marker IDs to contig/position, splits by linkage group (`split_lepmap.py`), and converts to phased genotype calls (`map2genotypes.awk`) — see `Linkage_map_assembly.Rmd` for the underlying logic this automates.

7. **QTL analysis:**
    - Render `Rqtl.Rmd` from repository root:
       ```bash
       conda activate linkage-mapping_rmd
       Rscript -e "rmarkdown::render('Rqtl.Rmd')"
       ```
    - Alternative full analysis is in `scripts/ovules_rqtl_full.Rmd`.

### Script Arguments (Current)

- `scripts/rmdup.sh`
   - `library_filenames DIR partition [time] [mem]`
- `scripts/flip.sh`
   - `library_filenames DIR partition`
- `scripts/process_radtags.sh`
   - `library_filenames DIR partition`
- `scripts/align_f2s.sh`
   - `sample_names DIR genome_dir genome partition`
- `scripts/rg_f2s.sh`
   - `sample_names DIR picard genome_dir genome partition`
- `scripts/align_parents.sh`
   - `sequences.txt DIR read_dir genome_dir genome partition`
- `scripts/rg_parents.sh`
   - `sequences.txt DIR picard partition`
- `scripts/bcftools.sh`
   - `listfiles bamfiles DIR genome_dir genome prefix partition`
- `scripts/concatenate_vcfs.sh`
   - `listfiles prefix bcftools_dir partition`
   - `listfiles` and `prefix` must match the values passed to `scripts/bcftools.sh`.
- `scripts/filter_vcf_for_lepmap3_ovule_cross.py`
   - Required flags:
      - `--invcf --ParentList --SampleList --out --minMapQ --minF2depth --missingList --MaxMissingData`
   - Optional boolean flags:
      - `--filterMissingParents --filterInbredParents --filterforhets --minparentdepth --maxparentdepth`
   - Optional threshold:
      - `--FractionF2` (default `0`)
- `scripts/lepmap.sh`
   - `[-p PRE] [-o OUT] [-l LOD] [-t THETA] [-d dT] [-s MINORSCAFS] [-q PARTITION]`
   - `-d` was previously broken (missing from the option parser); now fixed.
   - `-s` is new: controls the `minorscafs` argument passed to `split_lepmap.py` (default `False`).
   - The generated job now also re-maps marker IDs to contig/position, splits by linkage group, and produces phased genotype calls in the same submission (previously separate manual steps).
- `scripts/split_lepmap.py`
   - `mapping_order_file outprefix minorscafs`
   - Normally invoked automatically by the job script `scripts/lepmap.sh` generates, not run directly.

### Data Files

- **Input data:** Located in `data/`
- **Genotypes:** `ovule_f2s.finalsites.genotypes.txt`, `ovule_sub_genotypes_missing_lowhet.txt`
- **Sample metadata:** `ovule_f2_samples.txt`, `ovule_barcode_key.txt`, `ovule_parents.txt`
- **Marker/site references:** `im767.v2.list`, `ovule_sites.bed`

## Dependencies

### Software
- **R:** version 4.0.5 or later
- **RStudio:** for rendering R Markdown files
- **Python:** 3.8.12 (for linkage-mapping) and 2.7.18 (for py2 environment)
- **Git:** for version control
- **Conda:** for environment management

### R Packages
- `rqtl` (QTL mapping)
- `tidyverse` (data manipulation)
- `ggplot2` (visualization)
- `rmarkdown` (document rendering)

### Python Packages
- See `environment/linkage-mapping_environment.yml` for complete linkage mapping environment specs
- See `environment/py2_environment.yml` for py2 environment specs

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

- **Thom Nelson:** ddRAD deaggregation scripts (py2 environment)
- **Robert G. Franks:** Project supervision and guidance
- **Mimulus Research Community:** For resources and support

### Funding and Resources

This research was funded by the National Science Foundation grants IOS-1557692 (RF) and IOS-1558113 (JW). This work is supported by the Research Capacity Fund (HATCH), project award no. 7005173, from the U.S. Department of Agriculture’s National Institute of Food and Agriculture.

### Related Publications


