---
title: Data
subtitle: Access to fastq files and info for the dataset used for the course.
description: Access to fastq files and info for the dataset used for the course.
featured_image: /images/corte.jpg
---

## Course dataset

The course dataset comes from a collaboration between my lab, Phil Scott's lab, and our colleagues in Salvador, Brazil (manuscript [here](https://doi.org/10.1126/scitranslmed.aax4204)).  This is an RNAseq dataset from skin biopsies obtained from patients with cutaneous leishmaniasis, a parasitic disease endemic in Brazil and other areas of South America.  You'll be working with data from 5 patients with this disease and 5 healthy endemic controls.   

This course offers multiple 'on-ramps', or entry points, for starting an analysis.  Below are a few of these options.

---

### Option 1 (preferred)- Download raw data and accessory files

[fastq files](https://drive.google.com/drive/folders/1sEk1od1MJKLjqyCExYyfHc0n7DAIy_x7?usp=sharing) – You will need about **30Gb** of storage space on your harddrive to accomodate this download.  *please do not uncompress these files (leave them as .gz files)*.  

If you're unable to download the files above, you can always get them from **[here](https://www.ebi.ac.uk/ena/browser/view/PRJNA525604?show=reads)** where these files were deposited during publication. Be sure to select only the specific files we'll be using during the course, which you can see **[here](https://diytranscriptomics.com/images/ena_screenshot.png)**.

[study design file](https://drive.google.com/file/d/1sGkMk-ekZk3pVD0aT5qzvMBo3UYnZSyU/view?usp=sharing) - A simple text file that contains one row for each sample and columns for each variable in the study.

[shell script](https://drive.google.com/file/d/1sGkMk-ekZk3pVD0aT5qzvMBo3UYnZSyU/view?usp=sharing) – For carrying out read mapping and QC analysis for multiple samples

---

### Option 2 – Download data already mapped

In the event that you don't have enough free harddrive space to download the raw fastq files in Option 1 above, or if you have any problems installing or using Kallisto to map this raw data, you can bypass the Kallisto step and download mapped data to start working immediately. 

[Kallisto output](https://drive.google.com/file/d/1gjWWTmStP1VT2TeIOTKJwByFOpyDHBkg/view?usp=sharing) -  These files are available as a single compressed zip file.  Download and unzip this file to reveal 10 folders (each containing the Kallisto output from mapping each of the 10 fastq files above), along with 10 log files – one for each sample that was mapped.  You may notice that each folder contains several files.  Please do not move or rename these files.  During the course, we will discuss these Kallisto outputs.

---

## Can I work with my own data during the course?

Yes!  We encourage students to work with their own data during the course.  Nothing will keep you more enagaged and crystalize the content of the course quite like working on questions you *actually* care about.  That said, it is not a requirement that you come with data in hand.  Throughout the course, I will demonstrate every step in the analysis using the dataset above. If you are working with your own dataset and run into problems, I strongly recommend that you switch to the course dataset so you don't fall behind.  You can always return to working on your own data once you have a functional pipeline and understand the process.  

---

## Where can I access public RNA-seq data

In the second half of the course, you'll learn how to access hundreds of thousands of samples already pre-mapped using Kallisto!  These are available through several different resources, which we'll discuss, but we'll focus on the [ARCHS4 database from Avi Ma'ayan's lab](https://amp.pharm.mssm.edu/archs4/index.html).  You'll want to download the compressed ARCHS4 databases for mouse samples [here](https://www.dropbox.com/s/crz1ctpqhbpb1ox/mouse_matrix_v8.h5?dl=0).  Note that this download is ~12GB of data.