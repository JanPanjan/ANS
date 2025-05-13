*Metabarcoding* is a technique that combines [[High-throughput technologies|high throughput DNA sequencing]] with **computational analysis** to identify and characterize the genetic material of multiple organisms within a single environment or bulk sample. 

Unlike traditional barcoding that focuses on identifying a single organism with a specific genetic marker, this approach focuses on the entire organism comunity in the sample.

This revolutionized ecology, environmental monitoring and biodiversity research by enabling quick and good quality assessments of the *species composition*.

## More formal

It's defined as *automated identification of multiple species from a single bulk sample containing entire organisms or from a single environmental sample containing degraded DNA*.

==Microbiome== refers to *the study of the entirety of the microbial genetic material recovered directly from the environment* (environmental DNA - eDNA), also known as shotgun metagenomics. Provids information about composition and function of the microbial comunity.

==Microbiota== refers to the *taxonomic composition of the microbial comunity as determined by metabarcoding analysis*. Provides an answer to "who is there" when studying a microbial comunity.

## Properties of the DNA fragment used for the barcode

Metabarcoding relies on the amplification and sequencing of a short region of DNA (a barcode marker - a gene or part of a gene). It has to have some properties:

- variable between species
- conserved flanking sites for annealing universal primers (across species)
- a sufficiently short region that can be sequenced using newer [[High-throughput technologies]]

Most commonly used genes for metabarcoding:

- prokaryotes with 16S [[rRNA]] ^e8946e
- eukaryotes with 18S rRNA, [ITS1](Internal%20transcribed%20spacer%20(ITS)) and [ITS2](Internal%20transcribed%20spacer%20(ITS))

---

**Primers** bind somewhere around the barcode region so the **PCR** amplifies these fragments. 

**Adapters** (molecular barcodes/index sequences) are then bound to these fragments so the machine can detect them. They are **distinct** for each processed sample in the same run. ^efcd99

Once the sequencing is done, the **bioinformatics data analysis** part is next.

## Typical metabarcoding pipeline

[[Demultiplexing]] : separating sequences

General [[practical-notebook#Filtering and trimming based on quality parameters|quality control]] steps to trim and discard low quality reads.

[[Denoising sequences in metabarcoding|Denoising]] / [[Clustering of sequences in metabarcoding|clustering]]  : processing of sequences to identify unique biological sequences present in the sample.

- Taxonomy classification 
	- Alignment based methods 
	- Machine-learning based classification methods (the multinomial Naive Bayes machine learning classifier in q2-feature classifier) 
- Diversity analysis 
	- Alpha diversity 
	- Beta diversity