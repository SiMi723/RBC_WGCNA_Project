WGCNA of RBC Clogging in Lung Tissue
Project Objective
This project uses Weighted Gene Co-expression Network Analysis (WGCNA) to identify key gene networks and biological pathways associated with Red Blood Cell (RBC) clogging in lung tissue samples. This repository contains the full, reproducible script and results of the analysis.

Workflow
The analysis pipeline follows these steps:

Data Pre-processing: Raw RNA-seq counts from data/expression_data.csv are filtered and normalized using DESeq2's Variance Stabilizing Transformation.

Network Construction: A co-expression network is built and clustered into gene modules using the WGCNA R package.

Module-Trait Correlation: Gene modules are correlated with the RBC clogging phenotype (from data/clinical_traits.csv) to identify the most relevant module.

Functional Enrichment: The genes from the key module are analyzed for Gene Ontology (GO) term enrichment to understand their biological function.

How to Run This Analysis
Clone this repository.

Place your expression_data.csv and clinical_traits.csv files in the data/ directory.

Open the RBC_WGCNA_Project.Rproj file in RStudio.

Open and run the scripts/01_WGCNA_Analysis.R script from top to bottom.

All outputs (plots and tables) will be automatically saved to the results/ directory.

Key Results
1. Module-Trait Relationship
The analysis identified a gene module that is highly correlated with the RBC clogging phenotype. The heatmap below shows the correlation value and p-value for each module.

2. Functional Enrichment of Key Module
The most significant module was found to be enriched for biological processes highly relevant to the disease pathology.