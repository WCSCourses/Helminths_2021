# Project

Congratulation!! You have reached the final day of the course. This is where you will gather the skills you have learnt and practice them on a new problem. Exchange the knowledge, tips and solutions with your coursemates, and ask your instructors if you get stuck.  

The project files that you have downloaded should look like this for [project 3](https://wgcadvancedcourses.github.io/Helminths_2021/manuals/module_4_project_intro/module_4_project_introduction.html#proj3)

And this for [project 4](https://wgcadvancedcourses.github.io/Helminths_2021/manuals/module_4_project_intro/module_4_project_introduction.html#proj4)


First, we need to access the content of this file
Notice that it ends with .tar.gz
To unzip this kind of file, we do 
$ tar -xvf FILENAME

Have a look inside by doing
$ ls FILENAME/

You should see the read count files and a text file containing the metadata
The metadata is there for your reference. It tells you the accession number of the original fastq files. There is no need to download the fastq files for this Project session because
1)	The fastq files are HUGE!! To process these files, we normally work on a computer cluster or a high-performance computer server, not on individual laptop
2)	We have already downloaded and processed the fastq files for you, up to the point that it can be analysed on a laptop – these are the read count files you see for each sample

So, the steps that have already been done for you are
1)	Downloaded the fastq files (see appendix XXX if you want to know how to do this)
2)	Download the reference genome and annotation file (GTF) of the worm from WormBase ParaSite (see WormBase ParaSite module)
3)	Index the reference genome using hisat2 (see transcriptomics module)
4)	Mapped the fastq to the reference genome using hisat2 (see transcriptomics module)
5)	Convert SAM to BAM using samtools view (see transcriptomics module)
6)	Sort BAM files using samtools sort (see transcriptomics module)
7)	Count reads per features (in this case, reads per gene) using htseq-count (see transcriptomics module)

The steps that you will perform by yourself in this Projects are (most of them refer to the Transcriptomics module)
1)	Starting in Rstudio
2)	Import the read counts data into Rstudio
3)	Explore the data using various data visualisation
4)	Perform pair-wise comparison of any comparisons that you are interested in and produce relevant plots
5)	Look at your results, compare your results, follow the lead of your results, explore genes etc. 
6)	Perform functional analysis using GO term enrichment
a.	Tips!!!! Notice that we have provided a GO term annotation reference for you in the Transcriptomics modules, but not in this Project module. You will need to download this one by yourself using WormBaseParaSite BioMart. :evil:
b.	See appendix XXX on how to do this
c.	You are also provided with the R script to reformat the downloaded GO reference into the format the topGO require
d.	However, this is not a complete and ready to use script, and will need some editing (real life scenario when you borrow script from someone else) 
e.	Have a look through the script and try to figure out which part need to be edited. A simplest was (not necessarily the most effective way, the most effective way would be to turn this script into an R function which can be run from a command line – like the run_topGO.R script) to do this would be to copy the whole code from GO_formatting.R and place it into your current R script and test your editing there. 


