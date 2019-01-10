---
layout: page
title: Data
comments: false
---


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
The challenge: Your PI has just sent you a frantic email requesting a figure for a grant submission that will be submitted tomorrow morning. The data is from a collaborator studying calcium signaling in T cells and is key to the experiments proposed in Aim 3. Unfortunately, the collaborator has only had time to align the reasd and they have not actually analyzed the results. Luckily, you have almost completed your RNA-Seq Analysis workshop and are ready and willing to assist! You are given mapped data from 12 samples which includes, in triplicate (except for one sample), anti-CD3/28 stimulated (aCD328) and unstimulated (US) WT and STIM1/2 cDKO CD4SP thymocytes. The STIM1/2 cDKO cells have muted calcium signals and therefore serve as a model for understanding the role of calcium in T cell signaling. Your PI wants 2 heatmaps of the genes are differentially expressed due to STIM1/2 DKO in unstimulated and stimulated T cells. You are also curious what pathways are enriched within these differentially expressed genes. To complete this challenge, you will need to download the mapped files and study design, read it into the R environment, and use a combination of tools that you have learned in class over the past several weeks to answer the questions below. The first team to submit the most complete answer to the following questions by the end of class will win the challenge. Good luck!

* Complete a principal component analysis. What do PC1 and PC2 represent?
* How many genes are differentially expressed (absolute LogFC > 1, FDR < 0.01) due to anti-CD3/28 stimulation of WT cells? What about for the cDKO cells? Note: filter your data prior to differential gene expression analysis to remove genes that have less than 10 counts across all samples.
* How many genes are differentially expressed (absolute LFC > 0.59, FDR < 0.05) in unstimulated and stimulated T cells due to STIM1/2 deletion?
* Please send a heatmap of the genes that are differentially expressed across genotype for i) unstimulated and ii) stimulated T cells along with answers for the above questions to corbettberry@gmail.com.
* Bonus: What are three pathways enriched in WT relative to STIM1/2 DKO in the largest cluster of differentially expressed genes for stimulated T cells?


** Solution: [see this script](http://DIYtranscriptomics.github.io/Data/files/hackdash3_solution.R)



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
