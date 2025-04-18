OUR ACCESSION NUMBER: **SRR30833097**

---

Look for the NCBI paper with the given title.

NCBI paper: *Wolbachia pipientis modulates germline stem cells and gene expression associated with ubiquitination and histone lysine trimethylation to rescue fertility defects in Drosophila* - https://academic.oup.com/genetics/article/229/3/iyae220/7934994#510949052

In the *data availability* section look for BioProject number. Copy it and paste it into NCBI to find the corresponding bioproject.

Sequence read archive bioproject with id PRJNA1166928: *The Drosophila melanogaster endosymbiont Wolbachia pipientis manipulates expression of genes associated with ubiquitin and histone lysine trimethylation when rescuing the fertility defect of a bag-of-marbles mutation.* - https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1166928/

To get our dataset, click on the number beside *SRA Experiments*.

![[Vaje-1 2025-04-17 13.57.02.excalidraw]]

files - https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR30833097&display=metadata

# 1. downloading through NCBI

Here we can see information about our data. Click on the link under *Run*, where the accession number is located.

![[Vaje-1 2025-04-17 13.58.26.excalidraw]]

> When we have paired-end files (here it's single), we have to use the `--split-3` flag if we are downloading data using `fastq-dump` and `prefetch`. When downloading from SRA-explorer this is done automatically.

Here we can directly download the FASTQ file.

![[Vaje-1 2025-04-17 13.59.13.excalidraw]]

![[Vaje-1 2025-04-17 13.59.39.excalidraw]]

# 2. downloading through SRA Explorer

If our files are too big, NCBI won't allow us to use the website interface. Instead we can use this tool: https://sra-explorer.info

![[Vaje-1 2025-04-17 14.01.41.excalidraw]]

Select the files and add them to collection. Then we can use `wget` from command line to download the dataset using the link it gives us.

![[Vaje-1 2025-04-17 14.02.55.excalidraw]]

# 3. downloading through the command line

Or if we're real programmers, we can use the command line. We can use `fastq-dump` or it's faster alternative `fasterq-dump`.

```
$ prefetch SRR30833097
2025-04-17T12:07:19 prefetch.3.0.3: Current preference is set to retrieve SRA Normalized Format files with full base quality scores.
2025-04-17T12:07:19 prefetch.3.0.3: 1) Downloading 'SRR30833097'...
2025-04-17T12:07:19 prefetch.3.0.3: SRA Normalized Format file is being retrieved, if this is different from your preference, it may be due to current file availability.
2025-04-17T12:07:19 prefetch.3.0.3:  Downloading via HTTPS...
2025-04-17T12:07:55 prefetch.3.0.3:  HTTPS download succeed
2025-04-17T12:07:56 prefetch.3.0.3:  'SRR30833097' is valid
2025-04-17T12:07:56 prefetch.3.0.3: 1) 'SRR30833097' was downloaded successfully
2025-04-17T12:07:56 prefetch.3.0.3: 'SRR30833097' has 0 unresolved dependencies

$ fasterq-dump SRR30833097
spots read      : 4,537,365
reads read      : 4,537,365
reads written   : 4,537,365
```

[[practical-notebook]]

# Homework questions (not graded)

## What is the aim of the study, described in the scientific article?

## Provide some info about the dataset you will work with and which sequencing technology was used.

## Which kit was used for sequencing library preparation? Does this kit preserve strand information (stranded library) or not?

## What is the advantage of stranded mRNA library preparation compared to non-stranded library?