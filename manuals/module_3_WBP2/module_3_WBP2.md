# WormBase ParaSite part 2

## Table of contents
1. [Tools](#tools)
    * [BLAST](#blast)
    * [The genome browser](#genome_browser)
    * [Expression data](#expression_data)
    * [VEP](#vep)
      * [EXERCISE](#vep_exercise)
2. [Accessing WormBase ParaSite data programmatically](#programmatic_access)
    * [Working with sequence and annotation files](#files)
    * [The REST API](#api)


## Tools <a name="tools"></a>

### BLAST <a name="blast"></a>

### The genome browser <a name="genome_browser"></a>

A genome browser is a tool that allows you to visualise a genome assembly and its features, together with experimental data aligned to the genome. Genome browsers are useful in many respects. For example, they are used for visualising gene models in their genomic context, and also to assess the correctness of the model. They are also used for visualising functional genomics data, such as the results of ChIP-Seq and RNA-Seq experiments.

There are several commonly used genome browsers in bioinformatics, each with different features. Examples include:

* The Integrative Genomics Viewer (IGV) - this software can be downloaded and installed locally to work with local data without internet access.
* The UCSC Genome Browser - this browser hosts selected vertebrate species and some model organism (and closely related) species.
* Ensembl - this can be used to browse a large catalog of genomes across the tree of life. WormBase ParaSite has an instance of the Ensembl browser built in.
* Artemis/ Artemis Comparison Tool (ACT) - similar to IGV, this can also run locally without the internet. It can also be used to edit (i.e. change) gene models (to correct errors). ACT is a tool for visualising synteny (regions of similarity) between sequences.
* JBrowse - this is the genome browser that we’ll be using today. WormBase ParaSite has an instance of JBrowse for every genome that it hosts. The Apollo project is a well known extension of JBrowse, which, like Artemis, can be used to edit gene models.

#### Using JBrowse: basic functionality

In this example we’ll introduce the basic functionality of JBrowse, and demonstrate how to use the various tracks.

* Navigate to the S. mansoni genome page and select the “Genome Browser (JBrowse)” icon.

![](figures/jbrowse_1.png)

In the genome browser, each scaffold is represented from its 5-prime end to its 3-prime end (relative to the forward strand). You can navigate to different scaffolds using the drop down menu in the middle of the screen, or by typing coordinates into the text box. Different types of data aligned to the genome are represented as tracks. When you first open JBrowse, one track will be on by default: the reference gene set.

For this example, we’ll consider that you’re interested in the gene Smp_312440.

* Start by typing the identifier into the search box and clicking “Go” to navigate to the gene.
* Zoom in by clicking the large magnifying glass with a “+” symbol until the reference sequence resolves.

![](figures/jbrowse_2.png)

Here, you can see the forward and reverse DNA strands, together with the six possible translational reading frames (3 forward and 3 reverse).

* Zoom out again so that you have the whole gene model in your field of view.
* To extract sequence information about the gene, click the gene model such that a dialogue box pops up.

![](figures/jbrowse_3.png)

Scrolling down the content of the box, you can extract genomic or cDNA sequence, or the sequence of specific subfeatures (specific exons or UTRs, for example).

Alternatively, you may wish to extract the genomic sequence of a whole region:

* Click the arrow next to the “Reference sequence” track label in the top left of the screen, select “Save track data”, then download the sequence as a FASTA
file.

![](figures/jbrowse_4.png)

#### Tracks

We can also use JBrowse to view other types of data aligned to the genome. 

* Click the “select tracks” button in the top left of the screen.

![](figures/jbrowse_5.png)

For most species, in addition to the gene model (“Genome Annotation”) track, there are two additional main types of track:

1. Repeat regions tracks - repetitive regions of the genome are annotated as part of WormBase ParaSite’s production process.

2. RNASeq tracks - WormBase ParaSite has a process of finding and aligning RNASeq data in the sequencing archives for our species of interest. These can be useful, for example, for checking that a gene model is well supported by expression data, or seeing in which life stages, or under which conditions, a gene of interest is transcribed.

For species with a lot of publicly available data, such as S. mansoni, the easiest way to explore the samples that are available is by using the facets on the left hand side. When sequencing experiments are submitted to the archives, submitters are asked to provide metadata on the samples (that is, to describe the samples in a detailed and structured way). This data is used to classify the samples in JBrowse.

Let’s say you want to see in which life stages Smp_312440 is expressed.

* Click the “” facet and select ..
* Select a few of the available libraries and click “back to browser”.

![](figures/jbrowse_6.png)

Each track represents a different sequencing library, and shows the number of reads that have been aligned at each position. By mousing over the histogram, you can see the exact number of aligned reads at each base. We can see that a lot of the tracks show biological replicates of the same condition (single sex infections or mixed sex infections). We can use combination tracks to combine replicate tracks “on the fly”, so we use up less space on the screen.

* In the main menu at the top of the page, select “Track” and “Add combination track”.

A new empty track should appear. You can then drag and drop existing tracks in to combine them. When you add additional tracks, a dialogue box should appear for you to select the type of operation to use to combine them. For this example, we’ll choose “addition”: you’ll now see the total number of reads across both selected libraries that aligned at each region. Note that different set operations can be performed, including subtraction, multiplication and division; these might make sense depending on the tracks that are being combined and the information that you’re interested in.

![](figures/jbrowse_7.png)

As well as seeing that Smp_312440 is expressed in these conditions, we can use the coverage histograms to assess the quality of the gene model. Most parasitic worm genomes are annotated with automated pipelines. Whilst annotation algorithms can often be very accurate, they are not infallible. Most of the gene models that you look at in WormBase ParaSite will not have been checked by a human curator, so it is important not to take them as “truth” unless you verify that they agree with any evidence that is available.

![](figures/jbrowse_8.png)

In this case we can see that each of the exons in the gene model have got good RNASeq coverage, with no additional exons suggested by the RNASeq data.

#### Motif searching

It might be useful to have a quick, visual way of showing where certain motifs (short, defined DNA sequences) are found in the reference sequence. JBrowse offers a quick and flexible way to do this. We’ll demonstrate this by generating a track for the TATA box sequence (a sequence found in the promoter region of many eukaryotic genes). The consensus TATA sequence is TATA[A/T]A[A/T] (where [A/T] indicates that either A or T could be present at that position).

* In JBrowse, select “Track” from the main menu bar, followed by “Add sequence search track”.
* Type the motif that we’re searching for in the dialogue box, in this format: TATA[AT]A[AT], and tick “Treat as regular expression”. This means that the [AT] section of the motif will be interpreted as a regular expression (ie, the base in this position can be either A or T). We will limit our search to the forward strand, because that’s the strand that our gene of interest is on. Click “Search”.

![](figures/jbrowse_9.png)

Going back to the main JBrowse window, a new track has appeared with all instances of the motif marked. Zooming in to the 5-prime end of Smp_312440, we can see that one of these is well positioned to be our TATA box.

#### Visualising your own data


As well as looking at publicly available data, you can use WormBase ParaSite JBrowse to visualise your own data. We’ll demonstrate how to do this using a BAM file that we have provided for you. BAM is a type of sequence file, in this case of RNA sequencing data. BAM files are binary (i.e., compressed and not human readable) versions of SAM files. SAM files are tab-delimited text files; each line in a SAM file represents a sequencing read, and (optionally) a description of how that read is aligned to a reference sequence.

In the module 3 data directory you should find a file named somules_isoseq_sorted.bam. This is a binary file, so trying to read it as it is won’t be very informative. samtools is a useful software package for manipulating SAM and BAM files. We’ll use a samtools command to convert the BAM file to a SAM file so we can have a look at how it’s structured. Move to the module 3 data directory and type the following into your terminal:


    samtools view -h somules_isoseq_sorted.bam | less

The SAM file starts with a header section. All header lines begin with a ‘@’ character.

![](figures/jbrowse_10.png)

Move down through the file (by pressing the space bar) until you come to the alignment section. Here, each line represents a sequencing read (though be aware that the lines are long, so a single line will probably wrap around your terminal window a few times). Some of the key fields are labelled below:

![](figures/jbrowse_11.png)

The full SAM specification is available here: http://samtools.github.io/hts-specs/

If you’ve looked at RNA sequencing data before, you may notice something unusual about the reads in this file: they’re very long! Until recently, next generation sequencing reads were typically ~100bp in length, so transcripts had to be sequenced in short sections at high coverage and reconstructed computationally. This BAM file contains “IsoSeq” data, from the Pacific Biosciences platform, whereby full length transcripts have been sequences in their entirety.

Before we can visualise the file in JBrowse, we need to create an index. An index is another file that often accompanies a BAM file, and acts like a table of contents. Software such as JBrowse can look inside the index file and find where exactly in the corresponding BAM file it needs to look, without having to go through all of the reads (which would be computationally very expensive).

BAM index files should have exactly the same name as their corresponding BAM file, with the addition of a .bai suffix. We can index our BAM file using samtools. Type:

    samtools index somules_isoseq_sorted.bam
    
You should now see a file called somules_isoseq_sorted.bam.bai in your working directory. We can now load the file into WormBase ParaSite JBrowse.

![](figures/jbrowse_12.png)

To add the BAM track, select the “Track” menu option in the top left of the screen. Selecting “Open track file or URL” will open a dialogue box giving you an option to view a file that is either on your file system, or accessible via a URL. Select both the BAM file and the index file. JBrowse guesses the file type from the name, but we have an option to correct it if it gets it wrong. We can see that it’s right this time. Click “Open”.

Now we can see the IsoSeq reads aligned to the genome. Notice that IsoSeq data is stranded- this means that the library preparation protocol preserved information on which end of the RNA molecule was 5-prime and which end was 3-prime, so we can infer which strand of DNA it was transcribed from. This information is encoded in the BAM file, and JBrowse colours the reads accordingly: reads aligning to the forward strand are pink, and reads aligning to the reverse strand are purple.


### Expression data <a name="expression_data"></a>

### VEP <a name="vep"></a>

The final WormBase ParaSite tool that we’ll look at today is the Variant Effect Predictor, or VEP. A common approach to understanding the genetic basis of phenotypic differences is to identify genetic variants that are overrepresented in some populations of individuals. For example, you might sequence two populations of worm: one that is susceptible to a drug and one that is resistant to the drug. You could then identify genomic positions where each of these populations differs from the reference genome. VEP is a tool that allows you to predict what the consequences of these variants are: whether they fall within or near genes, and whether they result in a change to the amino acid sequence of a protein.

The standard file format for storing variation data is VCF (variant call format); this is another tab-delimited text format. In the module 3 data directory, we have provided you with a Hymenolepis microstoma VCF file to demonstrate how to use VEP. Have a look at the file first to see how it’s structured (you'll have to scroll down beyond the headers to see the data lines):

    less h_microstoma.vcf
    
![](figures/vep_1.png)

* From the WormBase ParaSite homepage, select “Tools” from the toolbar.
* From the “Tools” page, select Variant Effect Predictor

![](figures/vep_2.png)

* To submit a VEP job, just select the correct species (Hymenolepis microstoma), upload your VCF file and click “Run”.

![](figures/vep_3.png)

The pie charts give a summary of the consequences of the variants found in the file. Variants with coding consequences are found in the protein-coding sequence of genes, whilst variants with non-coding consequences are in intergenic regions or non-coding regions of genes. These variants could still be functionally important; for example, variants in non-coding regions near genes can have effects on expression levels.
You can explore the results interactively on the webpage, or download them to file.

#### VEP exercise <a name="vep_exercise"></a>

Download the VEP results from the example above as a “VEP file”. Use this file and the original VCF file to answer the following questions:

1. How many variants were there in the original dataset?

2. What is their distribution across the scaffolds of the H. microstoma genome (hint: count how many times each scaffold appears in the VCF file)?

3. What are the different types of consequence that are found in the file, and how often does each occur?

4. List the genes where a ‘stop gained’ variant is found.

5. You’re interested in one particular gene, HmN_002063100. Does it have any variants in the file, and what are the reported consequences? Now view the VCF file in JBrowse and visualise where the variants are in the gene model.

Hint: to view the VCF in JBrowse you first need to compress and index it. Do:

    bgzip h_microstoma.vcf && tabix -p vcf h_microstoma.vcf.gz

## Accessing WormBase ParaSite data programmatically <a name="programmatic_access"></a>

### Working with sequence and annotation files <a name="files"></a>

### The REST API <a name="api"></a>

The other way to query WormBase ParaSite data is programmatically, via the REST API (Application Programming Interface). An API is just another way to retrieve data from a server, but this time via scripts or commands. You make a request to the server, but rather than returning a webpage, it returns the data in a structured format. We offer data in JSON (JavaScript Object Notation) and XML (Extensible Markup Language), which are both commonly used formats for data exchange. They are structured, so good for writing programs to interpret and manipulate them, but also human readable.

There are a few situations where accessing WormBase ParaSite data via the API might be the better choice over BioMart or the website:

1. For queries that you’re likely to have to run multiple times (for example, with different datasets, or against different genomes)

2. For queries that plug into a larger pipeline, it might be convenient to retrieve the data in an easily computer-processable format

3. Some types of data are not available in BioMart (such as CEGMA and BUSCO scores), and can only be accessed via the website or the API

In an earlier exercise, you used the assembly statistics widget on the genome page to compare Brugia sp. genome assemblies. In this example, we’ll do the same for the Meloidogyne sp. assemblies, using the API.

* From the WormBase ParaSite home page, select “REST API” from the toolbar.


