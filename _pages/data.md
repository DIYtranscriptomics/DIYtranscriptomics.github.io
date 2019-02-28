---
title: Data
subtitle: Access to fastq files and info for the dataset used for the course.
description: Access to fastq files and info for the dataset used for the course.
featured_image: /images/datapoints.jpg
---

## Bring your own data

You'll get the most from this class if you B.Y.O.D. -- Bring your own data.  Nothing will keep you more enagaged and crystalize the content of the course quite like working on questions you *actually* care about.  That said, it is not a requirement that you come with data in hand.  Throughout the course, I will demonstrate every step in the analysis using a 'real' dataset, which you can download below and use to follow along with me in class.  If you run into serious problems following along with your own data, I strongly recommend that you use the dataset provided below.  You can always return to your data once you have a functional pipeline and understand the process.  


## Course dataset

The featured dataset for this iteration of the course comes courtesy of [Boris Striepen's lab](http://www.striepenlab.org/) and is unpublished, so please be respectful of this. 

[raw data](https://www.dropbox.com/sh/df58trgab010s55/AAAQ86KkKPzuqvGG-YoeISNEa?dl=0) which consists of 9 fastq files.  You will need about 30Gb of storage space on your harddrive to accomodate these file.  *please do not uncompress these files (leave them as .gz files)*

[basic project directory](http://DIYtranscriptomics.github.io/Data/files/MIDAS.zip) that you can use to organize this (or any other) RNAseq project.  This project directory comes preloaded with some QC data, as well as a study design file, and a shell script for automating read mapping.  We'll discuss these in more detail during class.

[prealigned data](https://drive.google.com/file/d/1KDbXsGT0EGW9qiVihesWvspYiYxSwW76/view?usp=sharing) - In the event that you have any problems installing or using Kallisto to map this raw data, I've already mapped this data to produce transcript-level abundance data, which you can download and start working with immediately.  These files are available as a single compressed file.  Unzipping this file will reveal 9 folders (each containing the Kallisto output from mapping each of the 9 fastq files above), along with 9 log files (one for each sample that was mapped).  You may notice that each folder contains several files.  Please leave these in place.  During the course, we will discuss what these files actually mean.

[myDGEList](http://DIYtranscriptomics.github.io/Data/files/myDGEList) - This 'digital gene expression list' contains the raw counts for the dataset describe above and  provides a convenient structure for working with RNAseq count data in R. 
