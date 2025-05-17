We'll use **FeatureCounts** tool. It's a submodule of the [Subread](https://subread.sourceforge.net/) package. 

## Stranded or unstranded library

Before checking feature counts, we should check if our library is **stranded**, since we couldn't find this out before. 

We can do this with `infer_experiment.py`, which is a part of [RSeQC](https://rseqc.sourceforge.net/) package. 

> [!tip] RSeQC: An RNA-seq Quality Control Package
>
> RSeQC package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data. 
> 
> Some basic modules quickly inspect sequence quality, nucleotide composition bias, ..., while **RNA-seq specific modules** evaluate sequencing saturation, ..., ==strand specificity==, ..., etc.

Create a new environment and install `rseqc`:

```bash
conda create -n rseq
conda activate rseq
conda install -c bioconda rseqc
```

You can read about the script [here](https://rseqc.sourceforge.net/#infer-experiment-py). It requires a `bam` and a `bed` file. The `bed` file can be produced with [BEDOPS](https://bedops.readthedocs.io/en/latest/index.html) with the script  [gtf2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gtf2bed.html).

It transforms a file from a [Gene Transfer Format (GTF)](http://mblab.wustl.edu/GTF22.html) to a [Browser Extensible Data(BED)](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) one.

```bash
gtf2bed 
```