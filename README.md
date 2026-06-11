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
   conda activate linkage-mapping
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
- `data/ovule.cross`: R/qtl cross object.
- `data/FinalPhenotypeCalls.txt`: Phenotypic data for F2 individuals.

## Usage

### Main Analysis Workflow

1. **Data Preprocessing:**
   - Run shell scripts in `scripts/` for sequence alignment and barcode processing:
     ```bash
     bash scripts/trim_ddrad.sh
     bash scripts/trim_parents.sh
     bash scripts/align_parents.sh
     bash scripts/rg_ovules.sh
     ```

2. **Genotyping and Mapping:**
   - Use `scripts/lepmap10.20.0.0001.sh` to run LepmAP3 for linkage mapping.
   - Run `scripts/bcftools_ovules.sh` to filter VCF files.

3. **QTL Analysis:**
   - Open `scripts/Ovules.Rmd` in RStudio to run the main QTL mapping analysis.
   - Render with:
     ```bash
     Rscript -e "rmarkdown::render('Ovules.Rmd')"
     ```
   - Alternative analysis in `scripts/ovules_rqtl_full.Rmd`.

4. **Results:**
   - Output figures and tables are saved to `figs/` and root directory.
   - See `ovule.qtl.txt` for QTL summary results.

### Data Files

- **Input data:** Located in `data/`
- **Genotypes:** `ovule_f2s.finalsites.genotypes.txt`, `ovule_sub_genotypes.txt`
- **Phenotypes:** `FinalPhenotypeCalls.txt`
- **Maps:** `ovule.newmap`, `ovule.cross`

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


