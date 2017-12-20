---
layout: page
title: Data
comments: false
image:
  feature: MGCscreen2.jpg
---

<p class="message">
You'll get the most from this class if you B.Y.O.D. -- Bring your own data.  Nothing will keep you more enagaged and crystalize the content of the course quite like working on questions you *actually* care about.  That said, it is not a requirement that you come with data in hand.  Throughout the course, I will demonstrate every step in the analysis using a 'real' dataset, which you can download below and use to follow along with me in class.  If you run into serious problems following along with your own data, I strongly recommend that you use the dataset provided below.  You can always return to your data once you have a functional pipeline and understand the process.  
</p>

## The response of the intestinal epithelium to infection with the protozoan parasite, *Cryptosporidium parvum*
* The featured dataset for this iteration of the course comes courtesy of [Boris Striepen's lab](http://www.striepenlab.org/) and is unpublished, so please be respectful of this.  
* Download the [raw data](https://www.dropbox.com/sh/df58trgab010s55/AAAQ86KkKPzuqvGG-YoeISNEa?dl=0), which consists of 9 fastq files.  You will need about 30Gb of storage space on your harddrive to accomodate these file.  *please do not uncompress these files (leave them as .gz files)*
* you also need this basic [study design file](http://DIYtranscriptomics.github.io/Software/files/Crypto_studyDesign.txt) that describes the experiment.
* In the event that you have any problems installing or using Kallisto to map this raw data, I've already mapped this data to produce transcript-level abundance data , which you can download and start working with immediately.  These files are available as a single compressed file [here](https://www.dropbox.com/s/iyvttrdcjs300i1/mapppedData.zip?dl=0).  Unzipping this file will reveal 9 folders (each containing the Kallisto output from mapping each of the 9 fastq files above).  You may notice that each folder contains several files.  Please leave these in place.  During the course, we will discuss what these files actually mean.

## HackDash #1
The challenge: your collaborator is interested in understanding how cells respond to infection with Respiratory Syncytial Virus (RSV).  Previous work from their lab has show that some cells infected with RSV harbor primarily full-length (FL) viral genomes, while other cells accumulate short 'defective' genome genomes (DVGs), while yet other cells accumulate a mix of the two.  They ask for your help in interpreting data from their recent sequencing experiment in which they profiled response of human cells to these different viral states in sorted FL-hi, DVGs-hi, intermediate, or neither (not infected).  Using only your collaborator's Kallisto alignments and study design file, your job is to identify the main sources of variance in the data.  The winner will be the first team to email me with an explanation of what source of variation are evident in their experiment. Please include images of PCA plots to help make your case.  [Download the files to get started!](https://drive.google.com/file/d/0B-uUeUVY3YYUSTl6ZmZfcElid28/view?usp=sharing).

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

