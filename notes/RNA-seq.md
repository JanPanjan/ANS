Possible path in RNA-seq data analysis depending on whether a reference genome or transcriptome is available or not.

1. **quality of reads** is checked and reads are preproccessed to remove low-quality data and artifacts if neccessary.
2. read's origin is identified by **alignment to a reference genome** (if available), otherwise to a **reference transcriptome**. And if this does not exist, it can be produced using *de novo* **transcriptome assembly**.
3. novel genes and transcripts are detected using **genome-guided transcriptome assembly** (can be skipped for known genes and transcripts).
4. gene and transcript expression is **quantified**.

After obtaining abundance estimates, expression differences between sample groups can be analyzed using statistical testing. 

This whole proccess is called [[differential expression]] (a ni RNA-seq...)

![[Pasted image 20250417124800.png]]

Different steps are performed by different programs (UNIX philosophy?), which may need specific data formats and external files.
