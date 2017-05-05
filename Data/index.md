---
layout: page
title: Data
comments: false
---

<p class="message">
You'll get the most from this class if you B.Y.O.D. -- Bring your own data.  Nothing will keep you more enagaged and crystalize the content of the course quite like working on questions you *actually* care about.  That said, it is not a requirement that you come with data in hand.  Throughout the course, I will demonstrate every step in the analysis using a 'real' dataset, which you can download below and use to follow along with me in class.  If you run into serious problems following along with your own data, I strongly recommend that you use the dataset provided below.  You can always return to your data once you have a functional pipeline and understand the process.  
</p>

## The response of the intestinal epithelium to infection with the protozoan parasite, *Cryptosporidium parvum*
* The featured dataset for this iteration of the course comes courtesy of [Boris Striepen's lab](http://cellbio.uga.edu/directory/faculty/boris-striepen) and is unpublished, so please be respectful of this.  
* Download the [raw data](https://drive.google.com/drive/folders/0B-uUeUVY3YYUOVNGeDZvTXg1Q1U?usp=sharing), which consists of 9 fastq files.  You will need about 30Gb of storage space on your harddrive to accomodate these file.  *please do not uncompress these files (leave them as .gz files)*
* you also need this basic [study design file](http://DIYtranscriptomics.github.io/Software/files/Crypto_studyDesign.txt) that describes the experiment.
* In the event that you have any problems installing or using Kallisto to map this raw data, I've already mapped this data to produce transcript-level abundance data , which you can download and start working with immediately.  These files are available as a single compressed file [here](https://drive.google.com/drive/folders/0B-uUeUVY3YYUdHBidl9CTVR2N2M?usp=sharing).  Unzipping this file will reveal 9 folders (each containing the Kallisto output from mapping each of the 9 fastq files above).  You may notice that each folder contains several files.  Please leave these in place.  During the course, we will discuss what these files actually mean.

## HackDash #1
The challenge: your collaborator is interested in understanding how cells respond to infection with Respiratory Syncytial Virus (RSV).  Previous work from their lab has show that some cells infected with RSV harbor primarily full-lenght (FL) viral genomes, while other cells accumulate short 'defective' genome genomes (DVGs), while yet other cells accumulate a mix of the two.  They ask for your help in interpreting data from their recent sequencing experiment in which they profiled the host cell response to these different viral states in sorted FL-hi, DVGs-hi, intermediate, or neither (not-infected).  Using only your collaborator's Kallisto alignments and study design file, your job is to identify the main sources of variance in the data.  The winner will be the first team to email me the % variance explained by PC1 and PC2, and an explanation of what these principle components reflect in the experment. Please include images of PCA plots to help make your case.  [Download the files to get started!](https://drive.google.com/file/d/0B-uUeUVY3YYUSTl6ZmZfcElid28/view?usp=sharing).

## HackDash #2
The challenge: A colleague has asked for your help to mine data produced from a very large RNAseq study of the parasitc worm, Schistosoma mansoni.  In this experiment, male (M), female (F), juvenile (J) and mixed sex (X) worms were recovered from infected mice at various timepoints (control, 3hr, 12hr, and 24hr) following in vivo treatment with a low dose of the frontline anti-parasitic drug, praziquantel.  Experiments were carried out with three different strains of worms: NMRI, LE, and LEPZQ.  To complete this challenge, you'll need to [Download the processed data in a text file to get started!](http://DIYtranscriptomics.github.io/Data/files/data.unfiltered.txt), read it into the R environment and begin using the tools you've learned in class for wrangling dataframes (NOTE: you do NOT need to carry out a formal differential gene expression analysis here).  The first team to submit the most complete answer to the following questions by the end of class will win the challenge.  Good luck! 

* what are the dimensions of this data frame?
* select only the columns containing annotation info and the expression data for the female, LE strain worms
* add new columns that show the average expression for the triplicates at each timepoint for these samples
* add new columns that show the Log2 Fold Change for each AVG column, compared to control samples
* arrange genes in descending order based on Log2 FC of 24hr vs control
* how many genes have a log2 FC of 2 or more for 24hr vs control?
* produce an html table using the DT package that shows the top 10 most DE genes (based on Log2 FC only) for 24hr vs control, and includes all AVG columns for the female LE worms, along with annotation data.  
* Send me your html table by email, and tell me how many genes met the FC cutoff above.

** Solution: [see this script](http://DIYtranscriptomics.github.io/Data/files/hackdash2_solution.R)


## HackDash #3
The challenge: Your PI has just sent you a frantic email requesting a figure for a grant submission that will be submitted tomorrow morning. The data is from a collaborator studying calcium signaling in T cells and is key to the experiments proposed in Aim 3. Unfortunately, the collaborator has only had time to map the data and they have not had an opportunity to complete any of the analysis themselves. Luckily, you have almost completed your RNA-Seq Analysis workshop and are ready and willing to assist! You are given mapped data from 12 samples which includes, in triplicate (except for one sample), anti-CD3/28 stimulated (aCD328) and unstimulated (US) WT and STIM1/2 cDKO CD4SP thymocytes. The STIM1/2 cDKO cells have muted calcium signals and therefore serve as a model for understanding the role of calcium in T cell signaling. Your PI wants 2 heatmaps of the genes are differentially expressed due to STIM1/2 DKO in unstimulated and stimulated T cells. You are also curious what pathways are enriched within these differentially expressed genes. To complete this challenge, you will need to download the mapped files and study design, read it into the R environment, and use a combination of tools that you have learned in class over the past several weeks to answer the questions below. The first team to submit the most complete answer to the following questions by the end of class will win the challenge. Good luck!

* Complete a principal component analysis. What do PC1 and PC2 represent?
* How many genes are differentially expressed (absolute LFC > 1, FDR < 0.01) due to anti-CD3/28 stimulation of WT cells? What about for the cDKO cells? Note: filter your data prior to differential gene expression analysis to remove genes that have less than 10 counts across all samples.
* How many genes are differentially expressed (absolute LFC > 0.59, FDR < 0.05) in unstimulated and stimulated T cells due to STIM1/2 deletion?
* Please send a heatmap of the genes that are differentially expressed across genotype for i) unstimulated and ii) stimulated T cells along with answers for the above questions to corbettberry@gmail.com.
* Bonus: What are three pathways enriched in WT relative to STIM1/2 DKO in the largest cluster of differentially expressed genes for stimulated T cells?


** Solution: [see this script](http://DIYtranscriptomics.github.io/Data/files/hackdash3_solution.R)


