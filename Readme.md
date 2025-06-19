## Comprehensive Genomic Analysis of Somatic Mutations in Head and Neck Cancer (HNC)

This repository contains an R-based workflow for analyzing somatic mutation data from cancer patients. The project includes survival analysis, mutation landscape visualization, pathway enrichment, copy number alteration profiling, mutational signature deconvolution, and structural variant analysis.

---

### Features

- Kaplan-Meier survival analysis by clinical features (e.g., HPV status)
- Mutation landscape visualization using `maftools`
- Driver gene and protein domain analysis
- Co-mutation and mutual exclusivity analysis
- KEGG and GO enrichment analysis of mutated genes
- Copy number alteration (CNA) profiling with segmentation data
- Tumor mutational burden (TMB) comparison
- Circos plots of structural variants (gene fusions)
- Mutational signature analysis with `deconstructSigs`
- Fusion recurrence analysis in primary vs metastatic samples

---


### Input Data Files

The dataset were taken from cBioPortal - https://www.cbioportal.org/study/summary?id=hnc_mskcc_2016

Different files included in the dataset are:

- `data_clinical_patient.txt` – Clinical metadata (patient-level)
- `data_clinical_sample.txt` – Sample-level metadata
- `data_mutations.txt` – Somatic mutation data (MAF-like format)
- `data_cna_hg19.seg` – Copy number segmentation file
- `data_sv.txt` – Structural variant (fusion) data

> These are assumed to be downloaded from a resource like **cBioPortal** or **TCGA**.

---


### Workflow 

The project includes the following analytical steps:

1. Clinical Data Preprocessing and Survival Analysis
2. Mutation Landscape Visualization using `maftools`
3. Driver Mutation and Functional Domain Analysis
4. Pathway Enrichment Analysis using `clusterProfiler`
5. Copy Number Alteration Analysis
6. Tumor Mutational Burden Analysis
7. Circos Plot of Structural Variants
8. Mutational Signature Deconvolution
9. Analysis of Recurrent Fusions in Metastatic Disease

---


### Dependencies

- R ≥ 4.1
- maftools
- survival, survminer
- clusterProfiler, org.Hs.eg.db
- circlize
- deconstructSigs
- ggplot2, dplyr, data.table

---


### Acknowledgments

- cBioPortal
- The Cancer Genome Atlas (TCGA)
- Bioconductor packages used in the analysis

---

