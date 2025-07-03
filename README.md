<!--
SPDX-FileCopyrightText: Open Energy Transition gGmbH, Ember, and contributors to the Ember Flexibility Study
SPDX-License-Identifier: CC-BY-4.0
-->

# Ember Flexibility Study

This repository contains the code and analysis for the **Ember Flexibility Study**, conducted in collaboration with [Ember](https://ember-climate.org/) and [Open Energy Transition (OET)](https://openenergytransition.org/). The study investigates clean flexibility options for Europe's energy system, building on the PyPSA-Eur framework. All results are computed from raw data and code to ensure full reproducibility.


This repository is a soft-fork of [OET-PyPSA-Eur](https://github.com/open-energy-transition/pypsa-eur) and contains the entire project **Clean Flexibility for Europe's Energy System** supported by [Open Energy Transition (OET)](https://openenergytransition.org/)<sup>*</sup>, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

This repository is maintained using [OET's soft-fork strategy](https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy). OET's primary aim is to contribute as much as possible to the open source (OS) upstream repositories. For long-term changes that cannot be directly merged upstream, the strategy organizes and maintains OET forks, ensuring they remain up-to-date and compatible with upstream, while also supporting future contributions back to the OS repositories.

---

# Repository structure

* `benchmarks`: will store `snakemake` benchmarks (does not exist initially)
* `config`: configurations used in the study
* `cutouts`: will store raw weather data cutouts from `atlite` (does not exist initially)
* `data`: includes input data that is not produced by any `snakemake` rule
* `doc`: includes all files necessary to build the `readthedocs` documentation of PyPSA-Eur
* `envs`: includes all the `mamba` environment specifications to run the workflow
* `logs`: will store log files (does not exist initially)
* `notebooks`: includes all the `notebooks` used for ad-hoc analysis
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `rules`: includes all the `snakemake`rules loaded in the `Snakefile`
* `resources`: will store intermediate results of the workflow which can be picked up again by subsequent rules (does not exist initially)
* `results`: will store the solved PyPSA network data, summary files and plots (does not exist initially)
* `scripts`: includes all the Python scripts executed by the `snakemake` rules to build the model

# Installation and Usage

## 1. Installation

Clone the repository:

    git clone https://github.com/open-energy-transition/{{repository}}

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Users may also prefer to use [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) or [conda](https://docs.conda.io/projects/conda/en/latest/index.html). Using `mamba`, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

Activate the newly created `{{project_short_name}}` environment:

    mamba activate {{project_short_name}}

## 2. Run the analysis

    snakemake -call

This will run all analysis steps to reproduce results and build the report.

To generate a PDF of the dependency graph of all steps `resources/dag.pdf` run:

    snakemake -c1 dag

---

# Contributing and Support

We strongly welcome anyone interested in contributing to this project. If you have any ideas, suggestions or encounter problems, feel invited to file issues or make pull requests on GitHub.

- To **discuss** with other users, organise projects, share news, and get in touch with the community you can use the [Discord server](https://discord.gg/AnuJBk23FU).
- For **bugs and feature requests**, please use the [GitHub Issues page](https://github.com/open-energy-transition/{{repository}}/issues).

# Licence

The code in this repository is released as free software under the [MIT License](https://opensource.org/licenses/MIT), see [`doc/licenses.rst`](doc/licenses.rst). However, different licenses and terms of use may apply to the various input data, see [`doc/data_sources.rst`](doc/data_sources.rst).