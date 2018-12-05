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
* In the event that you have any problems installing or using Kallisto to map this raw data, I've already mapped this data to produce transcript-level abundance data, which you can download and start working with immediately.  These files are available as a single compressed file [here](https://www.dropbox.com/s/zr33kb1z1axqfi1/mappedData.zip?dl=0).  Unzipping this file will reveal 9 folders (each containing the Kallisto output from mapping each of the 9 fastq files above), along with 9 log files (one for each sample that was mapped).  You may notice that each folder contains several files.  Please leave these in place.  During the course, we will discuss what these files actually mean.

## HackDash #1
The challenge: your collaborator is interested in understanding how cells respond to infection with Respiratory Syncytial Virus (RSV).  Previous work from their lab has show that some cells infected with RSV harbor primarily full-length (FL) viral genomes, while other cells accumulate short 'defective' genome genomes (DVGs), while yet other cells accumulate a mix of the two.  They ask for your help in interpreting data from their recent sequencing experiment in which they profiled response of human cells to these different viral states in sorted FL-hi, DVGs-hi, intermediate, or neither (not infected).  Using only your collaborator's Kallisto alignments and study design file, your job is to identify the main sources of variance in the data.  The winner will be the first team to email me with an explanation of what source of variation are evident in their experiment. Please include images of PCA plots to help make your case.  [Download the files to get started!](https://drive.google.com/file/d/0B-uUeUVY3YYUSTl6ZmZfcElid28/view?usp=sharing).

## HackDash #2
The challenge: Download [this data](https://drive.google.com/file/d/0B-uUeUVY3YYUdTZ2bTJCZ3VyNWs/view?usp=sharing) from basic experiment in which T cells were left unstimulated or stimulated with anti-CD3/28 in duplicate.  Import these files into R and use your the skills/methods we covered in class to answer the following questions.  The teams that accumulates the most points (regardless of whether your answer is submitted first or last) wins the challenge.

1.	What is the gene name of the 1015th row of the output of tximport? (2 points)
2.	What are the total number of counts for each sample within Txi_gene$counts (4 points)
3.	Create a dendrogram and include an image here. Describe how your samples are clustering. (4 points)
4.	Complete a PCA and include an image here. What do PC1 and PC2 correspond to? (5 points)
5.	How many genes are differentially expressed (FDR < 0.05, abs LFC > 1) when you compare unstimulated and stimulated T cells? (5 points)
6.	Create a heatmap of all differentially expressed genes. (5 points)

## HackDash #3

The challenge: word has spread that you are an expert at RNAseq data analysis, so a collaborator asks for your help in identifying the major transcriptional changes that occur in skin following infection with the protozoan parasite _Leishmania braziliensis_.  You have been given RNAseq data from the skin of 5 healthy volunteers and 5 patients with cutaneous leishmaniasis.  The data has already been aligned using Kallisto and you can download the Kallisto outputs [here](https://drive.google.com/open?id=1scfKwSRCpwp-2OM0xZVDpu2moeCKIGwo), but you will need to carry this data through exploratory analysis, differential expression, and functional enrichment analysis.  Your goal in the next two hours is to answer as many of the following questions as possible:

1. What is the gene name of the 12391th row of the output of tximport? **(2 points)**

2. What are the total number of counts for each sample within Txi_gene$counts? **(2 points)**

3. How many genes were filtered _out_ when you set a cutoff of >1 CPM in at least 5 samples? **(2 points)**

4. Complete a PCA and include an image here. What does PC1 correspond to? Are there any outliers in the dataset? Which sample was this? **(4 points)**

5. How many genes are differentially expressed (FDR =< 0.05, abs logFC >= 1) when you compare disease vs. healthy skin? How many are up and downregulated in disease relative to healthy skin? **(2 points)**

6. Create a heatmap of all differentially expressed genes (FDR =< 0.05, abs logFC >= 1). Note: you can earn a bonus point if you label your groups with color. **(4 points; 5 possible with bonus)**

7.  Comparing cutaneous leishmaniasis to healthy skin, and using GO BP terms and the DAVID web resource, which two _Annotation Clusters_ are most enriched  in patients with disease? **(5 points)**

8. Using GSEA software to investigate only the Biocarta portion of the C2 collection from MSigDB, answer the following questions:
a) What are the top 3 enriched pathways in _disease_, compared to healthy skin? Indicate normalized enrichment score and FDR for each. **(3 points)**
b) What are the top 3 enriched pathways in _healthy skin_, compared to disease.  Indicate normalized enrichment score and FDR for each. **(3 points)**
c) For the top most enriched pathway in disease, how many genes in the pathway make up the 'leading edge'? **(2 points)**
d) Give a brief description and interpretation of enriched signatures in healthy skin relative to disease. **(2 points)**

