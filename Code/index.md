---
layout: page
title: R scripts
---

The lectures will be paired with R scripts which will 'step' through the process of analyzing transcriptional profiling data from RNAseq experiments

----

|	Script Name	|	Comments	|
|---------|:-----------:|
[readMapping.sh](http://DIYtranscriptomics.github.io/Code/files/readMapping.sh) 	|	Shell script for 'walk away' mapping of reads with Kallisto
[Step1_preprocessingKallisto.R](http://DIYtranscriptomics.github.io/Code/files/Step1_preprocessingKallisto.R) 	|	Packages: [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html), [readr](https://cran.r-project.org/web/packages/readr/README.html), [organism-specific database](https://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData)
[Step2_dataExploration.R](http://DIYtranscriptomics.github.io/Code/files/Step2_dataExploration.R) 	|	Packages: [ggplot2](http://ggplot2.org/), [reshape2](http://had.co.nz/reshape/)
[Step3_dataWrangling.R](http://DIYtranscriptomics.github.io/Code/files/Step3_dataWrangling.R) 	|	Packages: [ggplot2](http://ggplot2.org/), [reshape2](http://had.co.nz/reshape/), [dplyr](http://genomicsclass.github.io/book/pages/dplyr_tutorial.html), [ggvis](http://ggvis.rstudio.com/)
[Step4_diffGenes.R](http://DIYtranscriptomics.github.io/Code/files/Step4_diffGenes.R) 	|	Packages: [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [SVA](https://bioconductor.org/packages/release/bioc/html/sva.html), [Sleuth](http://pachterlab.github.io/sleuth/)
[Step5_heatmap.R](http://DIYtranscriptomics.github.io/Code/files/Step5_heatmap.R) 	|	Packages: [gplots](https://cran.r-project.org/web/packages/gplots/index.html), [Rcolorbrewer](http://earlglynn.github.io/RNotes/package/RColorBrewer/index.html)
[Step6_geneSetAnalysis.R](http://DIYtranscriptomics.github.io/Code/files/Step6_geneSetAnalysis.R) 	|	Packages: [GSVA](http://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf), [GSEAbase](http://bioconductor.org/packages/release/bioc/html/GSEABase.html)
[markdownSummaryReport.Rmd](http://DIYtranscriptomics.github.io/Code/files/markdownSummaryReport.Rmd) 	|	Packages: [Rmarkdown](http://rmarkdown.rstudio.com/), [knitr](http://yihui.name/knitr/)
[server.R](http://DIYtranscriptomics.github.io/Code/files/server.R), [ui.R](http://DIYtranscriptomics.github.io/Code/files/ui.R), [global.R](http://DIYtranscriptomics.github.io/Code/files/global.R), [Instructions.md](http://DIYtranscriptomics.github.io/Code/files/Instructions.md), [Abstract.md](http://DIYtranscriptomics.github.io/Code/files/Abstract.md)	|	Packages: [shiny](http://shiny.rstudio.com/), [shinyapps](https://www.shinyapps.io/), [DT](https://rstudio.github.io/DT/)