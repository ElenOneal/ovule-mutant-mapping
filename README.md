# ovule-mutant-mapping

A repository for the analysis scripts and data used in the ovule mutant mapping project. This project includes shell scripts, Python scripts, and R Markdown files for data analysis and mapping.

## Table of Contents

- [Project Overview](#project-overview)
- [Getting Started](#getting-started)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Project Overview

This mapping project was used to discover the locus underlying an ovule mutant phenotype first observed in Mimulus nudatus inbred line (GMR2). This line was at least 6 generations inbred when some individuals were observed with an unusual ovule phenotype. Instead of exhibiting normal development, the ovules took on a "Medusa"-like appearance. Early qPCR work suggested that changes in expression in the M. nudatus homolog of BEL1 could be linked to this developmental change. In Arabidopsis, BEL1 is involved in the establishment of ovule identity. 

Rationale
We did QTL mapping in an F2 mapping population to discover the loci underlying the ovule mutant in M. nudatus GMR2. The F2 population was created by first crossing DHRO22 x GMR2, both at least 6 generations inbred, then selfing the F1. Phenotyping of F2s revealed that the ovule mutant appears to segregate at a 25-35% rate, suggesting a single locus.

An Rmarkdown file (Ovules.Rmd) was created to document analyses. Work was performed in R version 4.0.5 (2021-03-31).

A conda virtual environment (“f1mapping”) was employed for most programs used in this project. A description of the conda environment f1mapping, including package versions, can be found in f1mapping_environment.yml.  This environment used Python 3.8.12.

A separate conda environment ("py2") was employed when deaggregating the ddRAD mapping data for the F2s used in this project. "Py2" used Python 2.7.18.

This project includes a series of scripts written in Python and shell, along with R Markdown documents to perform statistical analysis and generate reproducible figures.

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

