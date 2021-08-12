[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/DIYtranscriptomics/DIYtranscriptomics.github.io)

# Overview of this repo

A semester-long 'do-it-yourself' (DIY) transcriptomics course that teaches basic principles of bioinformatics using the R programming language and RNA-seq datasets from infectious disease studies.  

# How to contribute to this course

We're always interested in interesting infectious disease datasets to use as the basis for new lab exercises.  If you think you have a great idea for a new lab we should offer, just draft as a markdown document (you can download any one of our lab markdown documents from the /Projects folder in this repo) and submit as a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) so that we can review.

We also recognize that bioinformatics is a constant learning experience, even for instructors.  So, if you spot something in our lecture content that you think should be changed, please reach out directly to the course director, Dr. Dan Beiting (beiting@upenn.edu).

# How to modify and deploy this course at your own institution

Because the entire course, even down to the website, is contained within this github repo, you can simply clone the repo and modify the contents to tailor it for your own purposes.  Of course, if you do this, we would appreciate proper attribution to DIYtranscriptomics.org and Dr. Dan Beiting (University of Pennsylvania).  The key folders and files that you'll want to pay attention to and potentially change are:

- /Code - contains all R scripts used throughout the course. 
- /Data - contains course datasets.
- /Reading - This is where you put PDFs and other files that are linked in the lecture and lab markdown files described above.
- /images - contains all the graphics used on the lecture and lab pages.  
- /\_posts - contains all markdown documents corresponding to lecture pages, which get rendered on the [course main page](https://diytranscriptomics.com/)
- /\_projects - contains markdown documents corresponding to in-person lab exercises rendered on the [course lab page](https://diytranscriptomics.com/lab/)
- /\_pages - contains markdown docs corresponding to the [about](https://diytranscriptomics.com/about), [data](https://diytranscriptomics.com/data), [scripts](https://diytranscriptomics.com/scripts), [help](https://diytranscriptomics.com/help) and [video](https://diytranscriptomics.com/video) pages linked at the top of the main course page.
- /\_includes - contains a google analytics file (analytics.html) that you'll want to delete or edit if you decide to clone this repo and launch your own course site
- /\_data - contains settings.yml file that determines the header, menu, colors, fonts and other settings for the main page.
- index.html - this is where you can set the text and image for the main page
