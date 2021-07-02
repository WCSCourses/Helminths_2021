# Project Introduction and Planning

In just a few days, you will have reached the final module of the course, where you will apply what you learn in this course on a new problem. 

Instead of telling you exactly what to work on, we have provided 4 datasets for 4 different projects(!!). Each group (your Zoom breakout room) will discuss which project they collectively want to work on. We will do this project selection on the second day so that we are aware of the task required at the end of the course, and to have more time for absorbing and planning whatever are needed to complete the task.

In each project, on Friday, we recommend that you work individually on your laptop (so that you get your own practice), and share with the group when you hit an issue, or when you find something interesting. 

In the end, the group will work together to produce one group presentation. Each group will then have 8-10 minutes to present their analysis and their answers to the project questions. 

It could be helpful to have a group workspace on a shared document (such as Google Doc, Google Jamboard) where your group can share resources and ideas during the course, or even draw a workflow. Google Slide can be useful for co-creating group slides for the presentation.

## Here are the 4 projects and their questions:

### 1)	_S. mansoni_ in vivo timecourse <a name="proj1"></a>
This dataset come from the same RNA-seq data as your transcriptomic module (Module 7). These are different stages of _S. mansoni_ from infected mice, from early day-6 worms through to adult worms from [Wangwiwatsin et al. (2020)](https://journals.plos.org/plosntds/article/authors?id=10.1371/journal.pntd.0007743). However, this same dataset have been mapped to a different version of the reference genome. The data you analysed on Module 7 were mapped to an older version. For this exercise, we can find out if our interesting results from the previous genome version still remain when we use a newer version. In addition, have a go at different pairwise comparisons that were not explored in Module 7, run GO term enrichment on the list of differentially expressed genes, and explore any genes that look particularly interesting.

**Key questions**
-	What happened as the parasites developed to the next stage?
-	Are there any distinct differences between V5 and V7 results?
-	Challenge question: identify house-keeping genes based on this dataset

### 2)	The genetic variation of _S. mansoni_ genome after extended time in lab life cycle. <a name="proj2"></a>
Both the reference genome and the worms used for RNA-seq data of [this paper](https://journals.plos.org/plosntds/article/authors?id=10.1371/journal.pntd.0007743) are of the same strain (Puerto Rico), but they have been maintained in labs for more than 10 years apart. Does this lead to any changes in the genome? Can we see anything genomic changes from the RNA-seq data? 

**Key questions**
-	Are there any genes that contain SNPs (a form of genetic changes) and how might these SNPs affect the parasites?
-	What are the pros and cons of using RNA-seq data to study genetic variations?
-	Challenge question: Discuss any points of concern if you are using a reference genome that are from a very different time and place from your specimens. 

### 3)	Changes in gene expression of male and female from juvenile to adult worms <a name="proj3"></a>
Here we have another set of _S. mansoni_ RNA-seq data for schistosomules (juvenile worms) from [Picard et al. (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5038963/) and adult male and female from [Lu et al. (2016)](https://www.nature.com/articles/srep31150). The dataset have been mapped to the current version of _S. mansoni_ genome. We can perform various pairwise comparisons to look at how female and male worms change during their development, or how they behave differently at different stages. 

**Key questions**
-	What are the key differences between female and male adult worms?
-	What are the key differences between female and male schistosomules?
-	Challenge question 1: Are there any changes that unique to either male or female?
-	Challenge question 2: These RNA-seq data came from different studies, is it ok to analyse them together? 

### 4) _Brugia malayi_ male and female of different stages <a name="proj4"></a>
Here we have RNA-seq of _B. malayi_ from larval stages to various ages of male and female from [Grote et al. (2017)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005357). In their paper, Grote et al. investigated the relationship between _B.malayi_ and their symbiont bacteria _Wolbachia_, but they did not explain much about the gene expression changes of the worms in the different life cycle stages of the parasite. There are various pair-wise comparisons we could do here. Also, if you are intested in dual-RNA-seq, we also provide the FASTQ files that will have sequencing reads from the _Wolbachia_. You could use these FASTQ files to do the mapping to a relevant _Wolbachia_ genome.

**Key questions**
-	What happened to male and female worms as they develop into adults?
-	At which point do gene expression in male and female worms become very different?
-	Challenge question: If you want to do dual-RNA-seq analysis, what are the steps in the workflow that need to be done differently (e.g. different tools, different genome resources, etc)? 

---

<img src="https://www.cdc.gov/dpdx/schistosomiasis/modules/Schistomes_LifeCycle_lg.jpg" width="700">

**Figure 8.1:** Schistosoma spp. life cycle


<img src="https://www.cdc.gov/parasites/images/lymphaticfilariasis/B_malayi_LifeCycle.gif" width="700">

**Figure 8.2:** Brugia malayi life cycle

---

PS.  The questions, especially the **Challenge questions**, may not have exactly "the correct answer" but they are to allow discussion within groups, and with instructors. Plus, perhaps chances to reflex and discuss their current/future genomics projects

---
