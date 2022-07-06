# Mouse spinal cord during neuropathic pain.

This is a repository for code used in the publication Niehaus et al. 2021 entitled "Spinal cord macrophages resolve nociceptive hypersensitivity after peripheral injury".

# Example script files for several components of the paper
  
  Code for the Shiny online application (zylka-lab.med.unc.edu):
  * `app.R` - This code should be straightforward, with instructions I used to generate this application available at 'https://rstudio.github.io/shiny/tutorial/'.
  * `avgbar_lattice.R` and `propbar_lattice.R` - Plotting functions
  
  A general pipeline for single-cell RNA-seq (scRNAseq) dataprocessing and the source file used to generate plots; adapted from Loo et al. 2019 and Shekhar et al.2018.
  * `class042518.R` - source file for plotting, data manipulation, testing, and clustering functions.
  * `scRNAseq` - a step by step guide to two rounds of clustering from initial genes x cells dataframe (available at GEO GSE134003). This pipeline was adapted from the above publications, whose pipelines are found at 'https://github.com/jeremymsimon/MouseCortex' and 'https://github.com/broadinstitute/BipolarCell2016'. Clustering began with with an initial round followed by subclustering of each initial cluster.
  
  
Example scripts of immunofluorescence analyses via imageJ (Fiji).
  * `autoAdjust.ijm` - autoAdjust macro for Fiji; must be installed for GFAP and IBA1 quantification
  * `GFAParea.ijm`   Quantifies GFAP area per ROI area from maximum intensity projection image (GFAP = channel 1)
  * `IBA1area_and_cells.ijm` Quantifies IBA1 area and number of IBA1 DAPI cells from maximum intensity projection image (DAPI = channel 1, IBA1 = channel 2)
  
Example Rscript for testing and visualization used in the manuscript
  * `ExampleBoxPlotCode.R`
  * `LinePlotExample.R`
  * `MacrophageTimepointANOVAs.R`
Direct any questions or concerns via email to kylius0@email.unc.edu
