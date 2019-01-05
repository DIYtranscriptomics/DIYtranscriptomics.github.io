---
layout: page
title: R scripts
image:
  feature: Rscript.png
---

The lectures will be paired with R scripts which will 'step' through the process of analyzing transcriptional profiling data from RNAseq experiments

----

|	Script Name	|	Comments	|
|---------|:-----------:|
[readMapping.sh](http://DIYtranscriptomics.github.io/Code/files/readMapping.sh) 	|	Shell script for 'walk away' mapping of reads with Kallisto
[Step1_preprocessingKallisto.R](http://DIYtranscriptomics.github.io/Code/files/Step1_preprocessingKallisto.R) 	|	Packages: [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html), [readr](https://cran.r-project.org/web/packages/readr/README.html), [organism-specific database](https://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData), [Sleuth](https://github.com/pachterlab/sleuth)
[Step2_dataExploration.R](http://DIYtranscriptomics.github.io/Code/files/Step2_dataExploration.R) 	|	Packages: [ggplot2](http://ggplot2.org/), [reshape2](http://had.co.nz/reshape/)
[Step3_dataWrangling.R](http://DIYtranscriptomics.github.io/Code/files/Step3_dataWrangling.R) 	|	Packages: [ggplot2](http://ggplot2.org/), [reshape2](http://had.co.nz/reshape/), [dplyr](http://genomicsclass.github.io/book/pages/dplyr_tutorial.html), [ggvis](http://ggvis.rstudio.com/), [scatterD3](https://github.com/juba/scatterD3)
[Step4_publicData](http://DIYtranscriptomics.github.io/Code/files/Step4_publicData.R) 	|	Packages: [rhdf5](http://bioconductor.org/packages/release/bioc/html/rhdf5.html), [slinky](https://github.com/VanAndelInstitute/slinky)
[Step5_diffGenes.R](http://DIYtranscriptomics.github.io/Code/files/Step5_diffGenes.R) 	|	Packages: [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [SVA](https://bioconductor.org/packages/release/bioc/html/sva.html)
[Step6_heatmap.R](http://DIYtranscriptomics.github.io/Code/files/Step6_heatmap.R) 	|	Packages: [gplots](https://cran.r-project.org/web/packages/gplots/index.html), [Rcolorbrewer](http://earlglynn.github.io/RNotes/package/RColorBrewer/index.html), [heatmaply](https://cran.r-project.org/web/packages/heatmaply/index.html), [d3heatmap](https://cran.r-project.org/web/packages/d3heatmap/index.html)
[Step7_geneSetAnalysis.R](http://DIYtranscriptomics.github.io/Code/files/Step7_geneSetAnalysis.R) 	|	Packages: [GSVA](http://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf), [GSEAbase](http://bioconductor.org/packages/release/bioc/html/GSEABase.html)
[Rmarkdown_example.Rmd](http://DIYtranscriptomics.github.io/Code/files/Rmarkdown_example.Rmd), [Txi_gene](http://DIYtranscriptomics.github.io/Code/files/Txi_gene), [myDGEList](http://DIYtranscriptomics.github.io/Code/files/myDGEList) 	|	Packages: [Rmarkdown](http://rmarkdown.rstudio.com/), [knitr](http://yihui.name/knitr/)
