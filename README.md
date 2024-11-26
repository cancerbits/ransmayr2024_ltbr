# Code repository for the analysis of single-cell RNA-seq data presented in Ransmayr et al., 2024

Christoph Hafemeister and Florian Halbritter

St. Anna Children's Cancer Research Institute (CCRI), Vienna, Austria

**Abstract**

Secondary lymphoid organs (SLOs) provide the confined microenvironment required
for stromal cells to interact with immune cells to initiate adaptive immune responses resulting in
B-cell differentiation. Here, we studied three patients from two families with functional
hyposplenism, absence of tonsils and complete lymph node aplasia, leading to recurrent bacterial
and viral infections. We identified biallelic loss-of-function mutations in LTBR encoding the
lymphotoxin beta receptor (LTβR), primarily expressed on stromal cells. LTβR-deficient patients
had hypogammaglobulinemia, diminished memory B cells, regulatory and follicular T-helper
cells, and dysregulated expression of several Tumor Necrosis Factor family members. B-cell
differentiation in an ex vivo co-culture system was intact, implying that the observed B-cell defects
were not intrinsic in nature, and instead result from LTβR-dependent stromal cell interaction
signaling critical for SLO formation. Collectively, we define a human inborn error of immunity
caused primarily by a stromal defect affecting the development and function of SLOs.

## Repository structure

* `project.Dockerfile` defines the environment used to run the code
* `config.yaml` is used to set paths 
* `R/` holds R function definitions and misc utility scripts
* `Rmd/` holds R markdown documents for the individual steps of the project
* `bash/` holds shell scripts to build and run the docker image, and to parse the config file

## Reproducing the results

The file `knit.R` calls all `Rmd/*.Rmd` files in order to reproduce the analysis.

Paths in the `config.yaml` file starting with "/path/to/" will have to be set.

## Links

**Paper:** Science Immunology, 22 Nov 2024, Vol 9, Issue 101, DOI: [10.1126/sciimmunol.adq8796](https://doi.org/10.1126/sciimmunol.adq8796)

**Data files:** Raw data files are available at [The European Genome-phenome Archive (EGA)](https://ega-archive.org/studies/EGAS00001007271).
