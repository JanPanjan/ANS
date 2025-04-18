The workflow assumes a reference genome is available.

# 1. quality control of reads

Analysis starts with raw sequence reads, typically in FASTQ format. This analysis checks the overall quality of the millions of reads. 
They are scaned for low-confidence bases, biased nucleotide composition, adapters, duplicates... 

Output are basic statistics, which guide the preproccessing decisions in the next step.

# 2. preproccessing of reads

Goal is to remote low-quality bases and artifacts (like adapter or library construction sequences and also poly A tails) from each read. They may also be trimmed because of their size.

# 3. aligning reads to a reference genome

Goal is to find a point of origin for every read. When a read is mapped, a sequence alignment is created. It's neccessary to have a ref-seq as one of the input files here in addition to preproccessed reads.
Mapping is intesive, so the genome sequence is often transformed and compressed into an index to speed up the proccess. Most common one in use is the Burrows-Wheeler transform.

Output is an alignment file, which lists the mapped reads and their mapping positions in the reference.

# 4. genome-guided transcriptome assembly

![[Pasted image 20250417125755.png]]