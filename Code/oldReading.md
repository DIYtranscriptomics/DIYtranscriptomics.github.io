---
layout: page
title: R scripts
---

The lectures will be paired with R scripts which will 'step' through the process of analyzing transcriptional profiling data from RNAseq experiments

----

|	Script Name	|	Comments	|
|---------|:-----------:|
[readMapping.sh](http://transcriptomicsworkshop.github.io/Code/files/readMapping.sh) 	|	Shell script for 'walk away' mapping of reads with Kallisto
[Step1_preprocessingKallisto.R](http://transcriptomicsworkshop.github.io/Code/files/Step1_preprocessingKallisto.R) 	|	Packages: [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html), [readr](https://cran.r-project.org/web/packages/readr/README.html), [organism-specific database](https://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData)
[Step2_dataExploration_part1.R](http://transcriptomicsworkshop.github.io/Rscripts/files/Step2_dataExploration_part1.R) 	|	base R (dist, hclust, prcomp), ggplot2, reshape2
[Step3_dataExploration_part2.R](http://transcriptomicsworkshop.github.io/Rscripts/files/Step3_dataExploration_part2.R)  	|	dplyr (filter, arrange, select), reshape2 (melt), ggplot2, ggviz
[Step4_diffGenes.R](http://transcriptomicsworkshop.github.io/Rscripts/files/Step4_diffGenes.R)  	|	Limma (topTable and DecideTests),
[Step5_heatmap.R](http://transcriptomicsworkshop.github.io/Rscripts/files/Step5_heatmap.R) 	|	heatmap.2
[Step6_geneSetAnalysis.R](http://transcriptomicsworkshop.github.io/Rscripts/files/Step6_geneSetAnalysis.R) 	|	GSEAbase, GSVA
[Step7_Rmarkdown.rmd](http://transcriptomicsworkshop.github.io/Rscripts/files/Supplementary_code.Rmd)	|	Rmarkdown, knitr
[server.R](http://transcriptomicsworkshop.github.io/RNAseqWeb-app/server.R), [ui.R](http://transcriptomicsworkshop.github.io/RNAseqWeb-app/ui.R), [global.R](http://transcriptomicsworkshop.github.io/RNAseqWeb-app/global.R), [Instructions.md](http://transcriptomicsworkshop.github.io/RNAseqWeb-app/Instructions.md), [Abstract.md](http://transcriptomicsworkshop.github.io/RNAseqWeb-app/Abstract.md)	|	shiny, shinyapps, devtools, DT, ggvis, dplyr