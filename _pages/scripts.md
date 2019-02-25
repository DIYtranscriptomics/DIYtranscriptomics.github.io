---
title: Scripts
subtitle: Access to all the scripts and code used throughout the course.
description: Access to all the scripts and code used throughout the course.
featured_image: /images/code2.png
---

## Scripts

|	Script Name	|	Packages	|
|---------|:-----------:|
[Step1_Sleuth_TxImport.R](http://DIYtranscriptomics.github.io/Code/files/Step1_Sleuth_TxImport.R) 	|	[tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html), [tidyverse](https://www.tidyverse.org/), [Sleuth](https://github.com/pachterlab/sleuth), [Biomart](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
[Step2_dataWrangling.R](http://DIYtranscriptomics.github.io/Code/files/Step2_dataWrangling.R) 	|	[ggplot2](http://ggplot2.org/), [reshape2](http://had.co.nz/reshape/)
[Step3_multivariate.R](http://DIYtranscriptomics.github.io/Code/files/Step3_multivariate.R) 	|	[ggplot2](http://ggplot2.org/), [reshape2](http://had.co.nz/reshape/), [dplyr](http://genomicsclass.github.io/book/pages/dplyr_tutorial.html), [ggvis](http://ggvis.rstudio.com/), [scatterD3](https://github.com/juba/scatterD3)
[Step4_publicData.R](http://DIYtranscriptomics.github.io/Code/files/Step4_publicData.R) 	|	[rhdf5](http://bioconductor.org/packages/release/bioc/html/rhdf5.html), [slinky](https://github.com/VanAndelInstitute/slinky)
[Step5_diffGenes.R](http://DIYtranscriptomics.github.io/Code/files/Step5_diffGenes.R) 	|	[limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [SVA](https://bioconductor.org/packages/release/bioc/html/sva.html)
[Step6_modules.R](http://DIYtranscriptomics.github.io/Code/files/Step6_modules.R) 	|	[gplots](https://cran.r-project.org/web/packages/gplots/index.html), [Rcolorbrewer](http://earlglynn.github.io/RNotes/package/RColorBrewer/index.html), [heatmaply](https://cran.r-project.org/web/packages/heatmaply/index.html), [d3heatmap](https://cran.r-project.org/web/packages/d3heatmap/index.html)
[Step7_functionalEnrichment.R](http://DIYtranscriptomics.github.io/Code/files/Step7_functionalEnrichment.R) 	|	[GSVA](http://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf), [GSEAbase](http://bioconductor.org/packages/release/bioc/html/GSEABase.html)
[Rmarkdown_example.Rmd](http://DIYtranscriptomics.github.io/Code/files/Rmarkdown_example.Rmd), [Txi_gene](http://DIYtranscriptomics.github.io/Code/files/Txi_gene), [myDGEList](http://DIYtranscriptomics.github.io/Code/files/myDGEList) 	|	[Rmarkdown](http://rmarkdown.rstudio.com/), [knitr](http://yihui.name/knitr/)

---

## Data

<p class="message">
You'll get the most from this class if you B.Y.O.D. -- Bring your own data.  Nothing will keep you more enagaged and crystalize the content of the course quite like working on questions you *actually* care about.  That said, it is not a requirement that you come with data in hand.  Throughout the course, I will demonstrate every step in the analysis using a 'real' dataset, which you can download below and use to follow along with me in class.  If you run into serious problems following along with your own data, I strongly recommend that you use the dataset provided below.  You can always return to your data once you have a functional pipeline and understand the process.  
</p>

### The response of the intestinal epithelium to infection with the protozoan parasite, *Cryptosporidium parvum*
* The featured dataset for this iteration of the course comes courtesy of [Boris Striepen's lab](http://www.striepenlab.org/) and is unpublished, so please be respectful of this.  
* Download the [raw data](https://www.dropbox.com/sh/df58trgab010s55/AAAQ86KkKPzuqvGG-YoeISNEa?dl=0), which consists of 9 fastq files.  You will need about 30Gb of storage space on your harddrive to accomodate these file.  *please do not uncompress these files (leave them as .gz files)*
* you also need to download this [basic project directory](http://DIYtranscriptomics.github.io/Data/files/MIDAS.zip) that you can use to organize this (or any other) RNAseq project.  This project directory comes preloaded with some QC data, as well as a study design file, and a shell script for automating read mapping.  We'll discuss these in more detail during class.
* In the event that you have any problems installing or using Kallisto to map this raw data, I've already mapped this data to produce transcript-level abundance data, which you can download and start working with immediately.  These files are available as a single compressed file [here](https://drive.google.com/file/d/1KDbXsGT0EGW9qiVihesWvspYiYxSwW76/view?usp=sharing).  Unzipping this file will reveal 9 folders (each containing the Kallisto output from mapping each of the 9 fastq files above), along with 9 log files (one for each sample that was mapped).  You may notice that each folder contains several files.  Please leave these in place.  During the course, we will discuss what these files actually mean.
