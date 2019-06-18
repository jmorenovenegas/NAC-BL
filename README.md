# NAC-BL

This experiment deals with a WGCNA analysis of preNAC samples. We will investigate the molecular mechanisms associated with pCR.

It is very similar to the experiment in exp2019-04-16. But now we filter the sample by BC subtuype and analyse the Basal-like subtype.


In this case we are integrating WGCNA and DIAMOnD. We will find gene modules correlated with pCR and use these genes as seeds in DIAMOnD to capture more genes. Since our expression data comes from a limited panel (ca. 700 genes), we need to expand the set of genes we find correlated with pCR to have a broader picture of the biological processes and pathways that are implicated in pCR.

Update: There is no module correlated with ChemoRes or pCR. I will redirect the experiment to check pre NAC and post NAC samples of basal-like subtype. This way I´ll get some insigth into the mechanims that underlies the changes caused by NAC.

Analysis steps:

(chemoRes_NAC-BL.R)  
- Data loading, preprocessing, normalization and cleaning.
- Network construction and module detection.
- Relate modules to clinical chemoresistance(-pCR).
- WGCNA relevant modules enrichment.
- Module expansion using DIAMOnD.py (user must set path to python source with the required packages).
- Subnetwork building and enrichment of expanded modules.
- Link communities analysis subnetworks preprocessing.  

(link_communities.R)  
- Link communities analysis

(select_biolprocess.R)
- Extract genes from a linkcomm enriched function
- Annotate which genes were present in the panel, and if it has been activated o inibited

Instructions to perform the analysis:

1. Install python required libraries (see below)
2. Set path to python source with required libraries inside chemoRes_NAC-BL.R. 
3. Run chemoRes_NAC-BL.R by console or in RStudio (recommended).
4. Run link_communities.R. You must pass .gml files obtained in the previous analysis as arguments.
5. Run select_biolprocess.R

Used R functions available in RCode directory.

Interactomes used in DIAMOnD available in Interactomes directory.

RCode/interactomes_processing.R is the script used to process the interactomes.

More information is available in the scripts.

R version 3.5.1 in RStudio
- NanoStringQCPro
- dplyr 
- WGCNA
- biomaRt
- clusterProfiler
- linkcomm
- ggplot2
- igraph

Python version used 2.7.15
- networkx
- scipy 
- numpy
