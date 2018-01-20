---
layout: page
title: Reading material
---

{% include _toc.html %}


There are no required readings for the course, but I've provided links to Ebooks, primary literature, youtube videos of lectures or technical blog posts for each topic we discuss in class.

-----------------------------

## Ebooks on general R/bioconductor

* **[The Art of R programming by Norman Matloff](http://DIYtranscriptomics.github.io/Reading/files/The Art of R Programming.pdf)**
* **[R Graphics Cookbook](http://hdl.library.upenn.edu/1017.12/1675994)**
* **[Bioinformatics Data Skills: Reproducible and Robust Research with Open Source Tools](http://hdl.library.upenn.edu/1017.12/1448673)**
* **[R for Data Science](http://r4ds.had.co.nz/)**, by Hadley Wickham
* **Stay tuned for more titles available through the PennVet Library**


-----------------------------

## Introduction to RNAseq data and technology

* review Illumina's sequencing by synthesis (SBS) method and take a look at [this great overview of some file formats associated with transcriptomics](http://binf.snipcademy.com/lessons/sequence-file-formats).

* Check out the Univ. of Orgeon's [RNA-seqlopedia](http://rnaseq.uoregon.edu/), an online resource for understanding the entire RNAseq workflow

-----------------------------

## Read mapping with Kallisto

_papers, blogs posts and videos on Kallisto_

* the [April 2016 Nature Biotech paper](http://DIYtranscriptomics.github.io/Reading/files/Kallisto.pdf) from Lior Pachter's lab describing Kallisto

* [Harold Pimentel's talk on alignment](https://www.youtube.com/watch?v=b4tVokh6Law)

* [Lior Pachter's blog post on Kallisto](https://liorpachter.wordpress.com/2015/05/10/near-optimal-rna-seq-quantification-with-kallisto/)

* a great [blog post](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html) that explains the principle of 'pseudoalignments' that Kallisto uses to map reads to transcripts

* Did you notice that Kallisto is using 'Expectation Maximization (EM)' during the alignment?  You can read more about what this is [here](http://DIYtranscriptomics.github.io/Reading/files/EM.pdf)

* Follow [Kallisto discussions/questions](https://groups.google.com/forum/#!forum/kallisto-sleuth-users) and [Kallisto announcements](https://groups.google.com/forum/#!forum/kallisto-sleuth-announcements) on Google groups

_more general info about ultra lightweight methods for transcript quantification_

* [This 2014 Nature Biotech paper](http://DIYtranscriptomics.github.io/Reading/files/Sailfish.pdf) described Sailfish, which implimented the first lightweight method for quantifying transcript expression.

* Rob Patro, the first author of the Sailfish paper, wrote a nice blog post called [Not quite alignments](http://robpatro.com/blog/?p=248), comparing and contrasting alignment-free methods used by Sailfish, Salmon and Kallisto.

* Read more about [Salmon](http://biorxiv.org/content/early/2015/10/03/021592), the newest lightweight aligment tool out from Rob Patro and Carl Kinsford.  [Check out the website too](https://combine-lab.github.io/salmon/).

* Don't understand what a de Bruijn graph is? Find out more in the [this primer from Nature Biotechnology](http://DIYtranscriptomics.github.io/Reading/files/deBruijn.pdf)

* [Lior Pachter's talk at Cold Spring Harbor](http://theleadingstrand.cshl.edu/Course/Keynote/2013/INFO/135)

* [Greg Grant's recent paper comparing different aligners](http://DIYtranscriptomics.github.io/Reading/files/gregGrant_aligners_natMeth.pdf).  This should be a helpful guide in choosing alignment software outside of what we used in class.


-----------------------------

## Understanding RNAseq count data

* [What the FPKM?](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/) - Blog post by Harold Pimentel discussing within sample normalization and the meaning of RNAseq expression units

* [Between sample normalization in RNAseq](https://haroldpimentel.wordpress.com/2014/12/08/in-rna-seq-2-2-between-sample-normalization/) - another great blog post from Harold Pimentel on between-sample normalization.

* _Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples_ [Theory in Biosciences, Dec 2012](http://DIYtranscriptomics.github.io/Reading/files/wagnerTPM.pdf)

* _Revisiting global Gene Expression Analysis_ [Cell, Oct 2012](http://DIYtranscriptomics.github.io/Reading/files/Revisiting global gene express.pdf). A great example of the perils of normalizing to total read depth.

* the [original manuscript](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25) describing the Trimmed Mean of M values (TMM) method for normalizing between samples.

-----------------------------

## Starting your analysis script

* Read the [vignette for the Tximport package](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html) that we'll use to read the Kallisto mapping results into R.

*  *[Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](http://f1000research.com/articles/4-1521/v2)* F1000Research, Dec 2015. This paper describes the Tximport package and it's application for handling transcript-level expression measurments from lightweight aligners (Salmon, Sailfish, Kallisto)


-----------------------------

## Exploring, graphing and wrangling your data in R

* Take a look at some of the various ways to graph your data and the underlying R code in this [catalog of R graphs](http://shiny.stat.ubc.ca/r-graph-catalog/)

* If you end up using R to make a lot of graphs, you will find the [R Graphics Cookbook](http://hdl.library.upenn.edu/1017.12/1675994) to be an important reference. It's available free to UPenn folks as an Ebook.

* Colors palettes are an often underappreciated aspect of making beautiful and informative plots in R. You can access a suite of color palettes using the [RColorBrewer package](http://colorbrewer2.org).  These palettes can be viewed in [this cheatsheet](http://DIYtranscriptomics.github.io/Reading/files/colorbrewerPalettes.pdf).  Unfortunately, these standard palettes often don't cut it, and you'll need custom palettes. For this, I love using [Sip](https://sipapp.io/) to pick, organize and access color palettes.  

* I mentioned various unsupervised methods for dimensional reduction of your data (PCA, MDS, T-SNE).  In particular, T-SNE has become popular in representing single-cell RNAseq data, but it is also one of the more complex visualization methods to understand.  Although we didn't discuss this in class, I wanted to include a link to a [great blog post describing T-SNE](http://distill.pub/2016/misread-tsne/), as well as [the original T-SNE paper](http://DIYtranscriptomics.github.io/Reading/files/TSNE.pdf).  Please familiarize yourself with these if you plan on using this visualziation method.

* Hadley Wickham (author of Dplyr, Reshape2 and ggplot2 packages) has a nice pre-print paper on the key aspects of making ['Tidy Data'](http://vita.had.co.nz/papers/tidy-data.pdf)

* This [dplyr](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf) cheatsheet is useful resource to have on hand

* Although I love working with graphs in R, sometimes I just can figure out how to produce the final graphic exactly the way I want it.  So, I also really like the program [DataGraph](http://www.visualdatatools.com/DataGraph/).  Incredibly powerful graphing program, and very inexpensive!

-----------------------------

## Differential Gene Expression

* _Differential analysis of RNA-seq incorporating quantification uncertainty_ [Nature Methods, June, 2017](http://DIYtranscriptomics.github.io/Reading/files/sleuth.pdf)

* _voom: precision weights unlock linear model analysis tools for RNA-seq read counts_ [Genome Biology, Feb, 2014](http://DIYtranscriptomics.github.io/Reading/files/voom.pdf). Describes one of the approaches for adjusting RNAseq count data based on the mean-variance relationship.

* _RNA-Seq workflow: gene-level exploratory analysis and differential expression_ [F1000, Oct 2015](http://f1000research.com/articles/4-1070/v1)

* [Lior Pachter's blog post on Sleuth](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/)

* Harold Pimentel's [talk on differential expression with RNAseq](https://www.youtube.com/watch?v=BRWj6re9iGc)

* Good website that walks through [how to set-up a design matrix](http://genomicsclass.github.io/book/pages/expressing_design_formula.html)

* _Count-based differential expression analysis of RNA sequencing data using R and Bioconductor_ [Nature Protocols, Aug 22, 2013](http://DIYtranscriptomics.github.io/Reading/files/nprot.2013.099.pdf).  This is great overview of the edgeR and DESeq packages, their use, and explains how each one approaches differential gene expression.

* The [Limma user's guide](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

* The [EdgeR user's guide](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).  See section 3.4 and 3.5 for details about how to modify your model.matrix function for a 'blocking' design.

-----------------------------

## Functional Enrichment Analysis

* [The what, where, how and why of gene ontology -- a primer for bioinformaticians.](http://DIYtranscriptomics.github.io/Reading/files/GO.pdf)  Briefings in Bioinformatics, Feb 2011

* [A nice blog post](httP;//mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html) describing the hypergeometric test and Fisher's exact test that are at the core of many functional enrichment approaches.

* The [original 2003 Nat. Methods paper describing Gene Set Enrichment Analysis (GSEA)](http://DIYtranscriptomics.github.io/Reading/files/Mootha2003_GSEA.pdf), and the [2005 PNAS paper](http://mootha.med.harvard.edu/PubPDFs/Subramanian2005.pdf) that formally detailed its usage.

* Not a fan of using the Broad Inst. GSEA program?  You can carry out self-contained and competitive GSEA in R using [ROAST](http://DIYtranscriptomics.github.io/Reading/files/ROAST.pdf) and [CAMERA](http://DIYtranscriptomics.github.io/Reading/files/CAMERA.pdf), respectively.

* The paper describing [Gene Set VARIATION Analysis (GSVA)](http://DIYtranscriptomics.github.io/Reading/files/GSVA.pdf), which I find useful for looking GSEA-type results across a heterogeneous dataset (like a large cohort of patients).

* The Molecular Signatures Database [(MSigDB)](http://software.broadinstitute.org/gsea/msigdb)

* [2016 Immunity Paper](http://DIYtranscriptomics.github.io/Reading/files/ImmuneSigDB.pdf) describing how the 'Immunological Signatures' collection (C7) in MSigDB was created.

* A handy link if you're looking for [mouse and human specific gene signature databases](http://bioinf.wehi.edu.au/software/MSigDB/)

* If you need examples of properly formatted files for carrying out GSEA, you can check out the [GSEA file formats page](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats).

* you know how I feel about Venn diagrams, so if you're interested in exploring interactions between many groups of genes, have a look at [this Nature Methods paper](http://DIYtranscriptomics.github.io/Reading/files/upSet_plot.pdf), the accompanying R package, [UpSetR](https://cran.r-project.org/web/packages/UpSetR/README.html), as well as the [UpSet website](http://caleydo.org/tools/upset/).  Note, there's a shiny app for this as well!

* Here's [my Datagraph template](https://drive.google.com/drive/folders/1J1Fvw-73BjXYnStvAv7XIkK3FqI6dsMe?usp=sharing) for making bubble diagrams to show functional enrichment results.

-----------------------------

## Producing dynamic reports with Rmarkdown

* read [this blog post](http://kbroman.org/knitr_knutshell/) from Karl Broman on Knitr

* As with Git and dplyr, you can keep this [RMarkdown cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf) handy for quick reference

* For an example of how markdown/knitr can be used with publications, checkout [this Supplementary code](http://DIYtranscriptomics.github.io/Reading/files/supplementaryCode.pdf) in a recent [Nature Immunology paper](http://DIYtranscriptomics.github.io/Reading/files/singleCellTranscriptome.pdf) or our recent [PLOS Pathogens paper](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1005347) 

-----------------------------

## Deploying your data to the web

-----------------------------