# scRNA_analysis
### Step 1: 10X single cell cell ranger analysis
Using cell ranger package to map the reads, choose bam file true for the downstream velocity analysis

### Step 1 : data cleaning
0. Ambient RNA Removal: Utilized the R package SoupX to remove ambient RNA, with the resulting counts stored in a new layer named `SoupX_counts`.  
1. Data Pre-processing and Concatenation: Performed initial data pre-processing and concatenated datasets using pre_processing.py.
2. Core scRNA-seq Processing: This stage included doublet removal, data normalization, scaling, dimensionality reduction, batch correction using BBKNN, and initial cell clustering.
3. In Silico Lineage Removal: Implemented a strategy to remove cell lineages. This involved two sub-steps:
i. Filtering out clusters identified as lineage-negative cells (e.g., B cells, NK/T cells, pDCs).
ii. Applied a secondary filter using specific marker genes to remove remaining non-myeloid cells:
    * B cells: ***Cd79b***, ***Pax5***
    * pDCs: ***Siglech***
    * Erythroid cells: ***Gypa***, ***Hbb-bs***, ***Alas2***
    * Mast cells: ***Mcpt8***, ***Pf4***
    * Granulocytes: ***Ly6g***
4. Re-clustering and Annotation: Following the filtering steps, cells were re-clustered, and the resulting clusters were annotated based on known gene markers.
Note: data cleaning process is modified from [**Single-cell best practices**](https://www.sc-best-practices.org/preamble.html)
