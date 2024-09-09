# RF_NMDS
 Ordination through Random Forest framework in Caret.

 This function will a 3D ordination using Random Forest for the positional coordinates. The function has a few inputs:
  1) Data matrix (i.e. B-values from a methylation or gene expression array)
  2) Limma outputs for differential analysis
  3) A pheno file that has basic metadata with a minimum of "phenotype" and "SampleID" columns for sample identification
  4) the absolute value LFC valie from the limma analysis

 Once the function completes the output will contain coordinates for later plotting.
