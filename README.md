# Single-cell-RNA-seq-data-of-sex-and-androgens
# Codes provided here are to define and visualize the DEGs, AASB-DEGs and significantly enriched biological pathways of each cell type of 17 tissues, as well as to calculate cell proportion of each cell type.
# 1. Cellranger_pipeline.sh: code script in Shell for cellranger pipeline which is executed on Linux System. 
# 2. Seurat_object_construction.R: code in R to construct and integrate seurat objects.
# 3. DEG_definition: code in R to define DEGs based on p-adjust < 0.05 and |Log2FC| > 0.5.
# 4. Significantly-enriched-biological-pathway: code in R to definine the significantly enriched biological pathways based on DEGs (as examplified by male-biased DEGs here, Gene Ontology: Biological Process, p < 0.01 & q < 0.01.).
# 5. AASB-DEG_definition.R: code in R to define androgen-associated sex-biased DEGs (AASB-DEGs) for each cell type across the 17 tissues, which comprises positive AASB-DEGs (the expression levels of which were male-biased and positively regulated by androgens) and negative AASB-DEGs (the expression levels of which were female-biased and negatively regulated by androgens).
