[![DOI](https://img.shields.io/badge/DOI-10.1073%2Fpnas.2100751118-blue.svg)](https://doi.org/10.1073/pnas.2100751118)

# GRINS: Genetic elements that recode assembly-line polyketide synthases and accelerate their diversification

This repository contains code asssociated with the following publication.

Please cite the code in this repository with the following:

    Nivina A, Herrera Paredes S, Fraser H & Khosla C. "GRINS: genetic
    elements that recode assembly-line polyketide synthases and
    accelerate their diversification". PNAS June 29, 2021 118 (26) e2100751118;
    https://doi.org/10.1073/pnas.2100751118

## Directories

* **annotations**: Code to annotate genomes.

* **checkm**: Code to run [checkM](https://github.com/Ecogenomics/CheckM).

* **detection**: Code to detect GRINS.

* **phylo**: Code for phylogenetic analysis.

* **conda_envs**: Conda environment .yml files to recreate the environments
in the publication.

## Detecting GRINS

Instructions for running the code used to predict GRINS in novel genomes are
in the [detection directory](detection). Bear in mind that the current
workflow requires unannotated fasta assemblies as input, and performs
de-novo annotation.

## Copyright & License

    (C) Copyright 2017-2021 Sur Herrera Paredes, Aleksandra Nivina

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
