## Introduction

## What is HISAT2?

HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome, and exome sequencing data) against the general human population (as well as against a single reference genome). Based on [GCSA](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6698337&tag=1) (an extension of [BWT](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform) for a graph), we designed and implemented a graph FM index (GFM), an original approach and its first implementation to the best of our knowledge. In addition to using one global GFM index that represents general population, HISAT2 uses a large set of small GFM indexes that collectively cover the whole genome (each index representing a genomic region of 56 Kbp, with 55,000 indexes needed to cover human population). These small indexes (called local indexes) combined with several alignment strategies enable effective alignment of sequencing reads. This new indexing scheme is called Hierarchical Graph FM index (HGFM). We have developed HISAT 2 based on the [HISAT](http://ccb.jhu.edu/software/hisat) and [Bowtie2](http://bowtie-bio.sf.net/bowtie2) implementations. HISAT2 outputs alignments in [SAM](http://samtools.sourceforge.net/SAM1.pdf) format, enabling interoperation with a large number of other tools (e.g. [SAMtools](http://samtools.sourceforge.net/), [GATK](http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit)) that use SAM. HISAT2 is distributed under the [GPLv3 license](http://www.gnu.org/licenses/gpl-3.0.html), and it runs on the command line under Linux, Mac OS X and Windows.

## Obtaining HISAT2

Download HISAT2 sources and binaries from the Releases sections on the right side. Binaries are available for Intel architectures (`x86_64`) running Linux, and Mac OS X.

## Building from source

Building HISAT2 from source requires a GNU-like environment with GCC, GNU Make and other basics. It should be possible to build HISAT2 on most vanilla Linux installations or on a Mac installation with [Xcode](http://developer.apple.com/xcode/) installed. HISAT2 can also be built on Windows using [Cygwin](http://www.cygwin.com/) or [MinGW](http://www.mingw.org/) (MinGW recommended). For a MinGW build the choice of what compiler is to be used is important since this will determine if a 32 or 64 bit code can be successfully compiled using it. If there is a need to generate both 32 and 64 bit on the same machine then a multilib MinGW has to be properly installed. [MSYS](http://www.mingw.org/wiki/msys), the [zlib](http://cygwin.com/packages/mingw-zlib/) library, and depending on architecture [pthreads](http://sourceware.org/pthreads-win32/) library are also required. We are recommending a 64 bit build since it has some clear advantages in real life research problems. In order to simplify the MinGW setup it might be worth investigating popular MinGW personal builds since these are coming already prepared with most of the toolchains needed.

First, download the [source package](https://daehwankimlab.github.io/hisat2/download/) from the Download section on the right side. Unzip the file, change to the unzipped directory, and build the HISAT2 tools by running GNU `make` (usually with the command `make`, but sometimes with `gmake`) with no arguments. If building with MinGW, run `make` from the MSYS environment.

HISAT2 is using the multithreading software model in order to speed up execution times on SMP architectures where this is possible. On POSIX platforms (like linux, Mac OS, etc) it needs the pthread library. Although it is possible to use pthread library on non-POSIX platform like Windows, due to performance reasons HISAT2 will try to use Windows native multithreading if possible.

For the support of SRA data access in HISAT2, please download and install the [NCBI-NGS](https://github.com/ncbi/ngs/wiki/Downloads) toolkit. When running `make`, specify additional variables as follow. `make USE_SRA=1 NCBI_NGS_DIR=/path/to/NCBI-NGS-directory NCBI_VDB_DIR=/path/to/NCBI-NGS-directory`, where `NCBI_NGS_DIR` and `NCBI_VDB_DIR` will be used in Makefile for `-I` and `-L` compilation options. For example, `$(NCBI_NGS_DIR)/include` and `$(NCBI_NGS_DIR)/lib64` will be used.

## Running HISAT2

## Adding to PATH

By adding your new HISAT2 directory to your [PATH environment variable](http://en.wikipedia.org/wiki/PATH_\(variable\)), you ensure that whenever you run `hisat2`, `hisat2-build` or `hisat2-inspect` from the command line, you will get the version you just installed without having to specify the entire path. This is recommended for most users. To do this, follow your operating system’s instructions for adding the directory to your [PATH](http://en.wikipedia.org/wiki/PATH_\(variable\)).

If you would like to install HISAT2 by copying the HISAT2 executable files to an existing directory in your [PATH](http://en.wikipedia.org/wiki/PATH_\(variable\)), make sure that you copy all the executables, including `hisat2`, `hisat2-align-s`, `hisat2-align-l`, `hisat2-build`, `hisat2-build-s`, `hisat2-build-l`, `hisat2-inspect`, `hisat2-inspect-s` and `hisat2-inspect-l`.

## Reporting

The reporting mode governs how many alignments HISAT2 looks for, and how to report them.

In general, when we say that a read has an alignment, we mean that it has a [valid alignment](#valid-alignments-meet-or-exceed-the-minimum-score-threshold). When we say that a read has multiple alignments, we mean that it has multiple alignments that are valid and distinct from one another.

By default, HISAT2 may soft-clip reads near their 5’ and 3’ ends. Users can control this behavior by setting different penalties for soft-clipping ([`--sp`](#hisat2-options-no-softclip)) or by disallowing soft-clipping ([`--no-softclip`]).

### Distinct alignments map a read to different places

Two alignments for the same individual read are “distinct” if they map the same read to different places. Specifically, we say that two alignments are distinct if there are no alignment positions where a particular read offset is aligned opposite a particular reference offset in both alignments with the same orientation. E.g. if the first alignment is in the forward orientation and aligns the read character at read offset 10 to the reference character at chromosome 3, offset 3,445,245, and the second alignment is also in the forward orientation and also aligns the read character at read offset 10 to the reference character at chromosome 3, offset 3,445,245, they are not distinct alignments.

Two alignments for the same pair are distinct if either the mate 1s in the two paired-end alignments are distinct or the mate 2s in the two alignments are distinct or both.

### Default mode: search for one or more alignments, report each

HISAT2 searches for up to N distinct, primary alignments for each read, where N equals the integer specified with the `-k` parameter. Primary alignments mean alignments whose alignment score is equal or higher than any other alignments. It is possible that multiple distinct alignments have the same score. That is, if `-k 2` is specified, HISAT2 will search for at most 2 distinct alignments. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM ‘secondary’ bit (which equals 256) set in its FLAGS field. See the [SAM specification](http://samtools.sourceforge.net/SAM1.pdf) for details.

HISAT2 does not “find” alignments in any specific order, so for reads that have more than N distinct, valid alignments, HISAT2 does not guarantee that the N alignments reported are the best possible in terms of alignment score. Still, this mode can be effective and fast in situations where the user cares more about whether a read aligns (or aligns a certain number of times) than where exactly it originated.

## Alignment summary

When HISAT2 finishes running, it prints messages summarizing what happened. These messages are printed to the “standard error” (“stderr”) filehandle. For datasets consisting of unpaired reads, the summary might look like this:

```
20000 reads; of these:
  20000 (100.00%) were unpaired; of these:
    1247 (6.24%) aligned 0 times
    18739 (93.69%) aligned exactly 1 time
    14 (0.07%) aligned >1 times
93.77% overall alignment rate
```

For datasets consisting of pairs, the summary might look like this:

```
10000 reads; of these:
  10000 (100.00%) were paired; of these:
    650 (6.50%) aligned concordantly 0 times
    8823 (88.23%) aligned concordantly exactly 1 time
    527 (5.27%) aligned concordantly >1 times
    ----
    650 pairs aligned concordantly 0 times; of these:
      34 (5.23%) aligned discordantly 1 time
    ----
    616 pairs aligned 0 times concordantly or discordantly; of these:
      1232 mates make up the pairs; of these:
        660 (53.57%) aligned 0 times
        571 (46.35%) aligned exactly 1 time
        1 (0.08%) aligned >1 times
96.70% overall alignment rate
```

The indentation indicates how subtotals relate to totals.

## Wrapper

The `hisat2`, `hisat2-build` and `hisat2-inspect` executables are actually wrapper scripts that call binary programs as appropriate. The wrappers shield users from having to distinguish between “small” and “large” index formats, discussed briefly in the following section. Also, the `hisat2` wrapper provides some key functionality, like the ability to handle compressed inputs, and the functionality for [`--un`](#hisat2-options-un), [`--al`](#hisat2-options-al) and related options.

It is recommended that you always run the hisat2 wrappers and not run the binaries directly.

## Small and large indexes

`hisat2-build` can index reference genomes of any size. For genomes less than about 4 billion nucleotides in length, `hisat2-build` builds a “small” index using 32-bit numbers in various parts of the index. When the genome is longer, `hisat2-build` builds a “large” index using 64-bit numbers. Small indexes are stored in files with the `.ht2` extension, and large indexes are stored in files with the `.ht2l` extension. The user need not worry about whether a particular index is small or large; the wrapper scripts will automatically build and use the appropriate index.

## Performance tuning

1. If your computer has multiple processors/cores, use `-p`  
    The [`-p`](#hisat2-options-p) option causes HISAT2 to launch a specified number of parallel search threads. Each thread runs on a different processor/core and all threads find alignments in parallel, increasing alignment throughput by approximately a multiple of the number of threads (though in practice, speedup is somewhat worse than linear).

## Command Line

### Setting function options

Some HISAT2 options specify a function rather than an individual number or setting. In these cases the user specifies three parameters: (a) a function type `F`, (b) a constant term `B`, and (c) a coefficient `A`. The available function types are constant (`C`), linear (`L`), square-root (`S`), and natural log (`G`). The parameters are specified as `F,B,A` - that is, the function type, the constant term, and the coefficient are separated by commas with no whitespace. The constant term and coefficient may be negative and/or floating-point numbers.

For example, if the function specification is `L,-0.4,-0.6`, then the function defined is:

If the function specification is `G,1,5.4`, then the function defined is:

See the documentation for the option in question to learn what the parameter `x` is for. For example, in the case if the [`--score-min`](#hisat2-options-score-min) option, the function `f(x)` sets the minimum alignment score necessary for an alignment to be considered valid, and `x` is the read length.

### Usage

```
hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]
```

### Main arguments

- `-x <hisat2-idx>`  
    The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final `.1.ht2` / etc. `hisat2` looks for the specified index first in the current directory, then in the directory specified in the `HISAT2_INDEXES` environment variable.

- `-1 <m1>`  
    Comma-separated list of files containing mate 1s (filename usually includes `_1`), e.g. `-1 flyA_1.fq,flyB_1.fq`. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in `<m2>`. Reads may be a mix of different lengths. If `-` is specified, `hisat2` will read the mate 1s from the “standard in” or “stdin” filehandle.

- `-2 <m2>`  
    Comma-separated list of files containing mate 2s (filename usually includes `_2`), e.g. `-2 flyA_2.fq,flyB_2.fq`. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in `<m1>`. Reads may be a mix of different lengths. If `-` is specified, `hisat2` will read the mate 2s from the “standard in” or “stdin” filehandle.

- `-U <r>`  
    Comma-separated list of files containing unpaired reads to be aligned, e.g. `lane1.fq,lane2.fq,lane3.fq,lane4.fq`. Reads may be a mix of different lengths. If `-` is specified, `hisat2` gets the reads from the “standard in” or “stdin” filehandle.

- `--sra-acc <SRA accession number>`  
    Comma-separated list of SRA accession numbers, e.g. `--sra-acc SRR353653,SRR353654`. Information about read types is available at `http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?sp=runinfo&acc=**sra-acc**&retmode=xml`, where **sra-acc** is SRA accession number. If users run HISAT2 on a computer cluster, it is recommended to disable SRA-related caching (see the instruction at [SRA-MANUAL](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration)).

- `-S <hit>`  
    File to write SAM alignments to. By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console).

### Options

#### Input options

- `-q`  
    Reads (specified with `<m1>`, `<m2>`, `<s>`) are FASTQ files. FASTQ files usually have extension `.fq` or `.fastq`. FASTQ is the default format. See also: [`--solexa-quals`](#hisat2-options-solexa-quals) and [`--int-quals`](#hisat2-options-int-quals).

- `--qseq`  
    Reads (specified with `<m1>`, `<m2>`, `<s>`) are QSEQ files. QSEQ files usually end in `_qseq.txt`. See also: [`--solexa-quals`](#hisat2-options-solexa-quals) and [`--int-quals`](#hisat2-options-int-quals).

- `-f`  
    Reads (specified with `<m1>`, `<m2>`, `<s>`) are FASTA files. FASTA files usually have extension `.fa`, `.fasta`, `.mfa`, `.fna` or similar. FASTA files do not have a way of specifying quality values, so when `-f` is set, the result is as if `--ignore-quals` is also set.

- `-r`  
    Reads (specified with `<m1>`, `<m2>`, `<s>`) are files with one input sequence per line, without any other information (no read names, no qualities). When `-r` is set, the result is as if `--ignore-quals` is also set.

- `-c`  
    The read sequences are given on command line. I.e. `<m1>`, `<m2>` and `<singles>` are comma-separated lists of reads rather than lists of read files. There is no way to specify read names or qualities, so `-c` also implies `--ignore-quals`.

- `-s/--skip <int>`  
    Skip (i.e. do not align) the first `<int>` reads or pairs in the input.

- `-u/--upto <int>`  
    Align the first `<int>` reads or read pairs from the input (after the [`-s`/`--skip`](#hisat2-options-s) reads or pairs have been skipped), then stop. Default: no limit.

- `-5/--trim5 <int>`  
    Trim `<int>` bases from 5’ (left) end of each read before alignment (default: 0).

- `-3/--trim3 <int>`  
    Trim `<int>` bases from 3’ (right) end of each read before alignment (default: 0).

- `--phred33`  
    Input qualities are ASCII chars equal to the [Phred quality](http://en.wikipedia.org/wiki/Phred_quality_score) plus 33. This is also called the “Phred+33” encoding, which is used by the very latest Illumina pipelines.

- `--phred64`  
    Input qualities are ASCII chars equal to the [Phred quality](http://en.wikipedia.org/wiki/Phred_quality_score) plus 64. This is also called the “Phred+64” encoding.

- `--solexa-quals`  
    Convert input qualities from [Solexa](http://en.wikipedia.org/wiki/Phred_quality_score) (which can be negative) to [Phred](http://en.wikipedia.org/wiki/Phred_quality_score) (which can’t). This scheme was used in older Illumina GA Pipeline versions (prior to 1.3). Default: off.

- `--int-quals`  
    Quality values are represented in the read input file as space-separated ASCII integers, e.g., `40 40 30 40`…, rather than ASCII characters, e.g., `II?I`…. Integers are treated as being on the [Phred quality](http://en.wikipedia.org/wiki/Phred_quality_score) scale unless [`--solexa-quals`](#hisat2-options-solexa-quals) is also specified. Default: off.

#### Alignment options

- `--n-ceil <func>`  
    Sets a function governing the maximum number of ambiguous characters (usually `N`s and/or `.`s) allowed in a read as a function of read length. For instance, specifying `-L,0,0.15` sets the N-ceiling function `f` to `f(x) = 0 + 0.15 * x`, where x is the read length. See also: [setting function options]. Reads exceeding this ceiling are [filtered out](#filtering). Default: `L,0,0.15`.

- `--ignore-quals`  
    When calculating a mismatch penalty, always consider the quality value at the mismatched position to be the highest possible, regardless of the actual value. I.e. input is treated as though all quality values are high. This is also the default behavior when the input doesn’t specify quality values (e.g. in [`-f`](#hisat2-options-f), [`-r`](#hisat2-options-r), or [`-c`](#hisat2-options-c) modes).

- `--nofw/--norc`  
    If `--nofw` is specified, `hisat2` will not attempt to align unpaired reads to the forward (Watson) reference strand. If `--norc` is specified, `hisat2` will not attempt to align unpaired reads against the reverse-complement (Crick) reference strand. In paired-end mode, `--nofw` and `--norc` pertain to the fragments; i.e. specifying `--nofw` causes `hisat2` to explore only those paired-end configurations corresponding to fragments from the reverse-complement (Crick) strand. Default: both strands enabled.

#### Scoring options

- `--mp MX,MN`  
    Sets the maximum (`MX`) and minimum (`MN`) mismatch penalties, both integers. A number less than or equal to `MX` and greater than or equal to `MN` is subtracted from the alignment score for each position where a read character aligns to a reference character, the characters do not match, and neither is an `N`. If [`--ignore-quals`](#hisat2-options-ignore-quals) is specified, the number subtracted quals `MX`. Otherwise, the number subtracted is `MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) )` where Q is the Phred quality value. Default: `MX` = 6, `MN` = 2.

- `--sp MX,MN`  
    Sets the maximum (`MX`) and minimum (`MN`) penalties for soft-clipping per base, both integers. A number less than or equal to `MX` and greater than or equal to `MN` is subtracted from the alignment score for each position. The number subtracted is `MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) )` where Q is the Phred quality value. Default: `MX` = 2, `MN` = 1.

- `--no-softclip`  
    Disallow soft-clipping.

- `--np <int>`  
    Sets penalty for positions where the read, reference, or both, contain an ambiguous character such as `N`. Default: 1.

- `--rdg <int1>,<int2>`  
    Sets the read gap open (`<int1>`) and extend (`<int2>`) penalties. A read gap of length N gets a penalty of `<int1>` + N * `<int2>`. Default: 5, 3.

- `--rfg <int1>,<int2>`  
    Sets the reference gap open (`<int1>`) and extend (`<int2>`) penalties. A reference gap of length N gets a penalty of `<int1>` + N * `<int2>`. Default: 5, 3.

- `--score-min <func>`  
    Sets a function governing the minimum alignment score needed for an alignment to be considered “valid” (i.e. good enough to report). This is a function of read length. For instance, specifying `L,0,-0.6` sets the minimum-score function `f` to `f(x) = 0 + -0.6 * x`, where `x` is the read length. See also: [setting function options]. The default is `L,0,-0.2`.

#### Spliced alignment options

- `--pen-cansplice <int>`  
    Sets the penalty for each pair of canonical splice sites (e.g. GT/AG). Default: 0.

- `--pen-noncansplice <int>`  
    Sets the penalty for each pair of non-canonical splice sites (e.g. non-GT/AG). Default: 12.

- `--pen-canintronlen <func>`  
    Sets the penalty for long introns with canonical splice sites so that alignments with shorter introns are preferred to those with longer ones. Default: G,-8,1

- `--pen-noncanintronlen <func>`  
    Sets the penalty for long introns with noncanonical splice sites so that alignments with shorter introns are preferred to those with longer ones. Default: G,-8,1

- `--min-intronlen <int>`  
    Sets minimum intron length. Default: 20

- `--max-intronlen <int>`  
    Sets maximum intron length. Default: 500000

- `--known-splicesite-infile <path>`  
    With this mode, you can provide a list of known splice sites, which HISAT2 makes use of to align reads with small anchors.  
    You can create such a list using `python hisat2_extract_splice_sites.py genes.gtf > splicesites.txt`, where `hisat2_extract_splice_sites.py` is included in the HISAT2 package, `genes.gtf` is a gene annotation file, and `splicesites.txt` is a list of splice sites with which you provide HISAT2 in this mode. Note that it is better to use indexes built using annotated transcripts (such as _genome_tran_ or _genome_snp_tran_), which works better than using this option. It has no effect to provide splice sites that are already included in the indexes.

- `--novel-splicesite-outfile <path>`  
    In this mode, HISAT2 reports a list of splice sites in the file
    
    > chromosome name `<tab>` genomic position of the flanking base on the left side of an intron `<tab>` genomic position of the flanking base on the right `<tab>` strand (+, -, and .)
    
    ’.’ indicates an unknown strand for non-canonical splice sites.
    

- `--novel-splicesite-infile <path>`  
    With this mode, you can provide a list of novel splice sites that were generated from the above option “–novel-splicesite-outfile”.

- `--no-temp-splicesite`  
    HISAT2, by default, makes use of splice sites found by earlier reads to align later reads in the same run, in particular, reads with small anchors (<= 15 bp).  
    The option disables this default alignment strategy.

- `--no-spliced-alignment`  
    Disable spliced alignment.

- `--rna-strandness <string>`  
    Specify strand-specific information: the default is unstranded.  
    For single-end reads, use F or R.
    
    > ‘F’ means a read corresponds to a transcript.  
    > ‘R’ means a read corresponds to the reverse complemented counterpart of a transcript.
    
    For paired-end reads, use either FR or RF.  
    With this option being used, every read alignment will have an XS attribute tag:
    
    > ’+’ means a read belongs to a transcript on ‘+’ strand of genome.  
    > ‘-‘ means a read belongs to a transcript on ‘-‘ strand of genome.
    
    (TopHat has a similar option, –library-type option, where fr-firststrand corresponds to R and RF; fr-secondstrand corresponds to F and FR.)
    

- `--tmo/--transcriptome-mapping-only`  
    Report only those alignments within known transcripts.

- `--dta/--downstream-transcriptome-assembly`  
    Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

- `--dta-cufflinks`  
    Report alignments tailored specifically for Cufflinks. In addition to what HISAT2 does with the above option (–dta), With this option, HISAT2 looks for novel splice sites with three signals (GT/AG, GC/AG, AT/AC), but all user-provided splice sites are used irrespective of their signals. HISAT2 produces an optional field, XS:A:[+-], for every spliced alignment.

- `--avoid-pseudogene`  
    Try to avoid aligning reads to pseudogenes. Note this option is experimental and needs further investigation.

- `--no-templatelen-adjustment`  
    Disables template length adjustment for RNA-seq reads.

#### Reporting options

- `-k <int>`  
    It searches for at most `<int>` distinct, primary alignments for each read. Primary alignments mean alignments whose alignment score is equal to or higher than any other alignments. The search terminates when it cannot find more distinct valid alignments, or when it finds `<int>`, whichever happens first. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM ‘secondary’ bit (which equals 256) set in its FLAGS field. For reads that have more than `<int>` distinct, valid alignments, **hisat2** does not guarantee that the `<int>` alignments reported are the best possible in terms of alignment score. Default: 5 (linear index) or 10 (graph index).
    
    Note: HISAT2 is not designed with large values for `-k` in mind, and when aligning reads to long, repetitive genomes, large `-k` could make alignment much slower.
    

- `--max-seeds <int>`  
    HISAT2, like other aligners, uses seed-and-extend approaches. HISAT2 tries to extend seeds to full-length alignments. In HISAT2, `--max-seeds` is used to control the maximum number of seeds that will be extended. For DNA-read alignment ([`--no-spliced-alignment`](#hisat2-options-no-spliced-alignment)), HISAT2 extends up to these many seeds and skips the rest of the seeds. For RNA-read alignment, HISAT2 skips extending seeds and reports no alignments if the number of seeds is larger than the number specified with the option, to be compatible with previous versions of HISAT2. Large values for `--max-seeds` may improve alignment sensitivity, but HISAT2 is not designed with large values for `--max-seeds` in mind, and when aligning reads to long, repetitive genomes, large `--max-seeds` could make alignment much slower. The default value is the maximum of 5 and the value that comes with `-k` times 2.

- `-a/--all`  
    HISAT2 reports all alignments it can find. Using the option is equivalent to using both [`--max-seeds`](#hisat2-options-max-seeds) and [`-k`](#hisat2-options-k) with the maximum value that a 64-bit signed integer can represent (9,223,372,036,854,775,807).

- `--secondary`  
    Report secondary alignments.

#### Paired-end options

- `-I/--minins <int>`  
    The minimum fragment length for valid paired-end alignments.This option is valid only with [`--no-spliced-alignment`](#hisat2-options-no-spliced-alignment). E.g. if `-I 60` is specified and a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with a 20-bp gap between them, that alignment is considered valid (as long as [`-X`](#hisat2-options-X) is also satisfied). A 19-bp gap would not be valid in that case. If trimming options [`-3`](#hisat2-options-3) or [`-5`](#hisat2-options-5) are also used, the [`-I`](#hisat2-options-I) constraint is applied with respect to the untrimmed mates.
    
    The larger the difference between [`-I`](#hisat2-options-I) and [`-X`](#hisat2-options-X), the slower HISAT2 will run. This is because larger differences between [`-I`](#hisat2-options-I) and [`-X`](#hisat2-options-X) require that HISAT2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), HISAT2 is very efficient.
    
    Default: 0 (essentially imposing no minimum)
    

- `-X/--maxins <int>`  
    The maximum fragment length for valid paired-end alignments. This option is valid only with [`--no-spliced-alignment`](#hisat2-options-no-spliced-alignment). E.g. if `-X 100` is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as [`-I`](#hisat2-options-I) is also satisfied). A 61-bp gap would not be valid in that case. If trimming options [`-3`](#hisat2-options-3) or [`-5`](#hisat2-options-5) are also used, the `-X` constraint is applied with respect to the untrimmed mates, not the trimmed mates.
    
    The larger the difference between [`-I`](#hisat2-options-I) and [`-X`](#hisat2-options-X), the slower HISAT2 will run. This is because larger differences between [`-I`](#hisat2-options-I) and [`-X`](#hisat2-options-X) require that HISAT2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), HISAT2 is very efficient.
    
    Default: 500.
    

- `--fr/--rf/--ff`  
    The upstream/downstream mate orientations for a valid paired-end alignment against the forward reference strand. E.g., if `--fr` is specified and there is a candidate paired-end alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment length constraints ([`-I`](#hisat2-options-I) and [`-X`](#hisat2-options-X)) are met, that alignment is valid. Also, if mate 2 appears upstream of the reverse complement of mate 1 and all other constraints are met, that too is valid. `--rf` likewise requires that an upstream mate1 be reverse-complemented and a downstream mate2 be forward-oriented. `--ff` requires both an upstream mate 1 and a downstream mate 2 to be forward-oriented. Default: `--fr` (appropriate for Illumina’s Paired-end Sequencing Assay).

- `--no-mixed`  
    By default, when `hisat2` cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.

- `--no-discordant`  
    By default, `hisat2` looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints ([`--fr/--rf/--ff`](#hisat2-options-fr), [`-I`](#hisat2-options-I), [`-X`](#hisat2-options-X)). This option disables that behavior.

#### Output options

- `-t/--time`  
    Print the wall-clock time required to load the index files and align the reads. This is printed to the “standard error” (“stderr”) filehandle. Default: off.

- `--un <path>`, `--un-gz <path>`, `--un-bz2 <path>`  
    Write unpaired reads that fail to align to file at `<path>`. These reads correspond to the SAM records with the FLAGS `0x4` bit set and neither the `0x40` nor `0x80` bits set. If `--un-gz` is specified, output will be gzip compressed. If `--un-bz2` is specified, output will be bzip2 compressed. Reads written in this way will appear exactly as they did in the input file, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the input.

- `--al <path>`, `--al-gz <path>`, `--al-bz2 <path>`  
    Write unpaired reads that align at least once to file at `<path>`. These reads correspond to the SAM records with the FLAGS `0x4`, `0x40`, and `0x80` bits unset. If `--al-gz` is specified, output will be gzip compressed. If `--al-bz2` is specified, output will be bzip2 compressed. Reads written in this way will appear exactly as they did in the input file, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the input.

- `--un-conc <path>`, `--un-conc-gz <path>`, `--un-conc-bz2 <path>`  
    Write paired-end reads that fail to align concordantly to file(s) at `<path>`. These reads correspond to the SAM records with the FLAGS `0x4` bit set and either the `0x40` or `0x80` bit set (depending on whether it’s mate #1 or #2). `.1` and `.2` strings are added to the filename to distinguish which file contains mate #1 and mate #2. If a percent symbol, `%`, is used in `<path>`, the percent symbol is replaced with `1` or `2` to make the per-mate filenames. Otherwise, `.1` or `.2` are added before the final dot in `<path>` to make the per-mate filenames. Reads written in this way will appear exactly as they did in the input files, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the inputs.

- `--al-conc <path>`, `--al-conc-gz <path>`, `--al-conc-bz2 <path>`  
    Write paired-end reads that align concordantly at least once to file(s) at `<path>`. These reads correspond to the SAM records with the FLAGS `0x4` bit unset and either the `0x40` or `0x80` bit set (depending on whether it’s mate #1 or #2). `.1` and `.2` strings are added to the filename to distinguish which file contains mate #1 and mate #2. If a percent symbol, `%`, is used in `<path>`, the percent symbol is replaced with `1` or `2` to make the per-mate filenames. Otherwise, `.1` or `.2` are added before the final dot in `<path>` to make the per-mate filenames. Reads written in this way will appear exactly as they did in the input files, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the inputs.

- `--quiet`  
    Print nothing besides alignments and serious errors.

- `--summary-file`  
    Print alignment summary to this file.

- `--new-summary`  
    Print alignment summary in a new style, which is more machine-friendly.

- `--met-file <path>`  
    Write `hisat2` metrics to file `<path>`. Having alignment metric can be useful for debugging certain problems, especially performance issues. See also: [`--met`](#hisat2-options-met). Default: metrics disabled.

- `--met-stderr`  
    Write `hisat2` metrics to the “standard error” (“stderr”) filehandle. This is not mutually exclusive with [`--met-file`](#hisat2-options-met-file). Having alignment metric can be useful for debugging certain problems, especially performance issues. See also: [`--met`](#hisat2-options-met). Default: metrics disabled.

- `--met <int>`  
    Write a new `hisat2` metrics record every `<int>` seconds. Only matters if either [`--met-stderr`](#hisat2-options-met-stderr) or [`--met-file`](#hisat2-options-met-file) are specified. Default: 1.

#### SAM options

- `--no-unal`  
    Suppress SAM records for reads that failed to align.

- `--no-hd`  
    Suppress SAM header lines (starting with `@`).

- `--no-sq`  
    Suppress `@SQ` SAM header lines.

- `--rg-id <text>`  
    Set the read group ID to `<text>`. This causes the SAM `@RG` header line to be printed, with `<text>` as the value associated with the `ID:` tag. It also causes the `RG:Z:` extra field to be attached to each SAM output record, with value set to `<text>`.

- `--rg <text>`  
    Add `<text>` (usually of the form `TAG:VAL`, e.g. `SM:Pool1`) as a field on the `@RG` header line. Note: in order for the `@RG` line to appear, [`--rg-id`](#hisat2-options-rg-id) must also be specified. This is because the `ID` tag is required by the [SAM Spec](http://samtools.sourceforge.net/SAM1.pdf). Specify `--rg` multiple times to set multiple fields. See the [SAM Spec](http://samtools.sourceforge.net/SAM1.pdf) for details about what fields are legal.

- `--remove-chrname`  
    Remove ‘chr’ from reference names in alignment (e.g., chr18 to 18)

- `--add-chrname`  
    Add ‘chr’ to reference names in alignment (e.g., 18 to chr18)

- `--omit-sec-seq`  
    When printing secondary alignments, HISAT2 by default will write out the `SEQ` and `QUAL` strings. Specifying this option causes HISAT2 to print an asterisk in those fields instead.

#### Performance options

- `-o/--offrate <int>`  
    Override the offrate of the index with `<int>`. If `<int>` is greater than the offrate used to build the index, then some row markings are discarded when the index is read into memory. This reduces the memory footprint of the aligner but requires more time to calculate text offsets. `<int>` must be greater than the value used to build the index.

- `-p/--threads NTHREADS`  
    Launch `NTHREADS` parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing `-p` increases HISAT2’s memory footprint. E.g. when aligning to a human genome index, increasing `-p` from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if HISAT2 is linked with the `pthreads` library.

- `--reorder`  
    Guarantees that output SAM records are printed in an order corresponding to the order of the reads in the original input file, even when [`-p`](#hisat2-options-p) is set greater than 1. Specifying `--reorder` and setting [`-p`](#hisat2-options-p) greater than 1 causes HISAT2 to run somewhat slower and use somewhat more memory then if `--reorder` were not specified. Has no effect if [`-p`](#hisat2-options-p) is set to 1, since output order will naturally correspond to input order in that case.

- `--mm`  
    Use memory-mapped I/O to load the index, rather than typical file I/O. Memory-mapping allows many concurrent `bowtie` processes on the same computer to share the same memory image of the index (i.e. you pay the memory overhead just once). This facilitates memory-efficient parallelization of `bowtie` in situations where using [`-p`](#hisat2-options-p) is not possible or not preferable.

#### Other options

- `--qc-filter`  
    Filter out reads for which the QSEQ filter field is non-zero. Only has an effect when read format is [`--qseq`](#hisat2-options-qseq). Default: off.

- `--seed <int>`  
    Use `<int>` as the seed for pseudo-random number generator. Default: 0.

- `--non-deterministic`  
    Normally, HISAT2 re-initializes its pseudo-random generator for each read. It seeds the generator with a number derived from (a) the read name, (b) the nucleotide sequence, (c) the quality sequence, (d) the value of the [`--seed`](#hisat2-options-seed) option. This means that if two reads are identical (same name, same nucleotides, same qualities) HISAT2 will find and report the same alignment(s) for both, even if there was ambiguity. When `--non-deterministic` is specified, HISAT2 re-initializes its pseudo-random generator for each read using the current time. This means that HISAT2 will not necessarily report the same alignment for two identical reads. This is counter-intuitive for some users, but might be more appropriate in situations where the input consists of many identical reads.

- `--version`  
    Print version information and quit.

- `-h/--help`  
    Print usage information and quit.

## SAM output

Following is a brief description of the [SAM](http://samtools.sourceforge.net/SAM1.pdf) format as output by `hisat2`. For more details, see the [SAM format specification](http://samtools.sourceforge.net/SAM1.pdf).

By default, `hisat2` prints a SAM header with `@HD`, `@SQ` and `@PG` lines. When one or more [`--rg`](#hisat2-options-rg) arguments are specified, `hisat2` will also print an `@RG` line that includes all user-specified [`--rg`](#hisat2-options-rg) tokens separated by tabs.

Each subsequent line describes an alignment or, if the read failed to align, a read. Each line is a collection of at least 12 fields separated by tabs; from left to right, the fields are:

1. Name of read that aligned. Note that the [SAM specification](http://samtools.sourceforge.net/SAM1.pdf) disallows whitespace in the read name. If the read name contains any whitespace characters, HISAT2 will truncate the name at the first whitespace character. This is similar to the behavior of other tools.
2. Sum of all applicable flags.  
    Flags relevant to HISAT2 are
    
    - 1: The read is one of a pair
    - 2: The alignment is one end of a proper paired-end alignment
    - 4: The read has no reported alignments
    - 8: The read is one of a pair and has no reported alignments
    - 16: The alignment is to the reverse reference strand
    - 32: The other mate in the paired-end alignment is aligned to the reverse reference strand
    - 64: The read is mate 1 in a pair
    - 128: The read is mate 2 in a pair
    
    Thus, an unpaired read that aligns to the reverse reference strand will have flag 16. A paired-end read that aligns and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1).
    
3. Name of reference sequence where alignment occurs
4. 1-based offset into the forward reference strand where leftmost character of the alignment occurs
5. Mapping quality. Mapping quality of HISAT2
6. CIGAR string representation of alignment
7. Name of reference sequence where mate’s alignment occurs. Set to `=` if the mate’s reference sequence is the same as this alignment’s, or `*` if there is no mate.
8. 1-based offset into the forward reference strand where leftmost character of the mate’s alignment occurs. Offset is 0 if there is no mate.
9. Inferred fragment length. Size is negative if the mate’s alignment occurs upstream of this alignment. Size is 0 if the mates did not align concordantly. However, size is non-0 if the mates aligned discordantly to the same chromosome.
10. Read sequence (reverse-complemented if aligned to the reverse strand)
11. ASCII-encoded read qualities (reverse-complemented if the read aligned to the reverse strand). The encoded quality values are on the [Phred quality](http://en.wikipedia.org/wiki/Phred_quality_score) scale and the encoding is ASCII-offset by 33 (ASCII char `!`), similarly to a [FASTQ](http://en.wikipedia.org/wiki/FASTQ_format) file.
12. Optional fields. Fields are tab-separated. `hisat2` outputs zero or more of these optional fields for each alignment, depending on the type of the alignment:
    - `AS:i:<N>` : Alignment score. Can be negative. Only present if SAM record is for an aligned read.
    - `ZS:i:<N>` : Alignment score for the best-scoring alignment found other than the alignment reported. Can be negative. Only present if the SAM record is for an aligned read and more than one alignment was found for the read. Note that, when the read is part of a concordantly-aligned pair, this score could be greater than [`AS:i`](#hisat2-opt-fields-as).
    - `YS:i:<N>` : Alignment score for opposite mate in the paired-end alignment. Only present if the SAM record is for a read that aligned as part of a paired-end alignment.
    - `XN:i:<N>` : The number of ambiguous bases in the reference covering this alignment. Only present if SAM record is for an aligned read.
    - `XM:i:<N>` : The number of mismatches in the alignment. Only present if SAM record is for an aligned read.
    - `XO:i:<N>` : The number of gap opens, for both read and reference gaps, in the alignment. Only present if SAM record is for an aligned read.
    - `XG:i:<N>` : The number of gap extensions, for both read and reference gaps, in the alignment. Only present if SAM record is for an aligned read.
    - `NM:i:<N>` : The edit distance; that is, the minimal number of one-nucleotide edits (substitutions, insertions and deletions) needed to transform the read string into the reference string. Only present if SAM record is for an aligned read.
    - `YF:Z:<S>` : String indicating reason why the read was filtered out. See also: [Filtering]. Only appears for reads that were filtered out.
    - `YT:Z:<S>` : Value of `UU` indicates the read was not part of a pair. Value of `CP` indicates the read was part of a pair and the pair aligned concordantly. Value of `DP` indicates the read was part of a pair and the pair aligned discordantly. Value of `UP` indicates the read was part of a pair but the pair failed to aligned either concordantly or discordantly.
    - `MD:Z:<S>` : A string representation of the mismatched reference bases in the alignment. See [SAM](http://samtools.sourceforge.net/SAM1.pdf) format specification for details. Only present if SAM record is for an aligned read.
    - `XS:A:<A>` : Values of `+` and `-` indicate the read is mapped to transcripts on sense and anti-sense strands, respectively. Spliced alignments need to have this field, which is required in Cufflinks and StringTie.  
        We can report this field for the canonical-splice site (GT/AG), but not for non-canonical splice sites. You can direct HISAT2 not to output such alignments (involving non-canonical splice sites) using “–pen-noncansplice 1000000”.
    - `NH:i:<N>` : The number of mapped locations for the read or the pair.
    - `Zs:Z:<S>` : When the alignment of a read involves SNPs that are in the index, this option is used to indicate where exactly the read involves the SNPs. This optional field is similar to the above MD:Z field. For example, `Zs:Z:1|S|rs3747203,97|S|rs16990981` indicates the second base of the read corresponds to a known SNP (ID: rs3747203). 97 bases after the third base (the base after the second one), the read at 100th base involves another known SNP (ID: rs16990981). ‘S’ indicates a single nucleotide polymorphism. ‘D’ and ‘I’ indicate a deletion and an insertion, respectively.

## The `hisat2-build` indexer

`hisat2-build` builds a HISAT2 index from a set of DNA sequences. `hisat2-build` outputs a set of 6 files with suffixes `.1.ht2`, `.2.ht2`, `.3.ht2`, `.4.ht2`, `.5.ht2`, `.6.ht2`, `.7.ht2`, and `.8.ht2`. In the case of a large index these suffixes will have a `ht2l` termination. These files together constitute the index: they are all that is needed to align reads to that reference. The original sequence FASTA files are no longer used by HISAT2 once the index is built.

Use of Karkkainen’s [blockwise algorithm](http://portal.acm.org/citation.cfm?id=1314852) allows `hisat2-build` to trade off between running time and memory usage. `hisat2-build` has three options governing how it makes this trade: [`-p`/`--packed`], [`--bmax`](#hisat2-build-options-bmax)/[`--bmaxdivn`](#hisat2-build-options-bmaxdivn), and [`--dcv`](#hisat2-build-options-dcv). By default, `hisat2-build` will automatically search for the settings that yield the best running time without exhausting memory. This behavior can be disabled using the [`-a`/`--noauto`](#hisat2-build-options-a) option.

The indexer provides options pertaining to the “shape” of the index, e.g. [`--offrate`](#hisat2-build-options-o) governs the fraction of [Burrows-Wheeler](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform) rows that are “marked” (i.e., the density of the suffix-array sample; see the original [FM Index](http://en.wikipedia.org/wiki/FM-index) paper for details). All of these options are potentially profitable trade-offs depending on the application. They have been set to defaults that are reasonable for most cases according to our experiments. See [Performance tuning](#performance-tuning) for details.

`hisat2-build` can generate either [small or large indexes](#small-and-large-indexes). The wrapper will decide which based on the length of the input genome. If the reference does not exceed 4 billion characters but a large index is preferred, the user can specify [`--large-index`](#hisat2-build-options-large-index) to force `hisat2-build` to build a large index instead.

The HISAT2 index is based on the [FM Index](http://en.wikipedia.org/wiki/FM-index) of Ferragina and Manzini, which in turn is based on the [Burrows-Wheeler](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform) transform. The algorithm used to build the index is based on the [blockwise algorithm](http://portal.acm.org/citation.cfm?id=1314852) of Karkkainen.

## Command Line

Usage:

```
hisat2-build [options]* <reference_in> <ht2_base>
```

### Notes

```
If you use --snp, --ss, and/or --exon, hisat2-build will need about 200GB RAM for the human genome size as index building involves a graph construction. 
Otherwise, you will be able to build an index on your desktop with 8GB RAM.
```

### Main arguments

- `<reference_in>`  
    A comma-separated list of FASTA files containing the reference sequences to be aligned to, or, if [`-c`](#hisat2-build-options-c) is specified, the sequences themselves. E.g., `<reference_in>` might be `chr1.fa,chr2.fa,chrX.fa,chrY.fa`, or, if [`-c`](#hisat2-build-options-c) is specified, this might be `GGTCATCCT,ACGGGTCGT,CCGTTCTATGCGGCTTA`.

- `<ht2_base>`  
    The basename of the index files to write. By default, `hisat2-build` writes files named `NAME.1.ht2`, `NAME.2.ht2`, `NAME.3.ht2`, `NAME.4.ht2`, `NAME.5.ht2`, `NAME.6.ht2`, `NAME.7.ht2`, and `NAME.8.ht2` where `NAME` is `<ht2_base>`.

### Options

- `-f`  
    The reference input files (specified as `<reference_in>`) are FASTA files (usually having extension `.fa`, `.mfa`, `.fna` or similar).

- `-c`  
    The reference sequences are given on the command line. I.e. `<reference_in>` is a comma-separated list of sequences rather than a list of FASTA files.

- `--large-index`  
    Force `hisat2-build` to build a [large index](#small-and-large-indexes), even if the reference is less than ~ 4 billion nucleotides long.

- `-a/--noauto`  
    Disable the default behavior whereby `hisat2-build` automatically selects values for the [`--bmax`](#hisat2-build-options-bmax), [`--dcv`](#hisat2-build-options-dcv) and [`--packed`] parameters according to available memory. Instead, user may specify values for those parameters. If memory is exhausted during indexing, an error message will be printed; it is up to the user to try new parameters.

- `--bmax <int>`  
    The maximum number of suffixes allowed in a block. Allowing more suffixes per block makes indexing faster, but increases peak memory usage. Setting this option overrides any previous setting for [`--bmax`](#hisat2-build-options-bmax), or [`--bmaxdivn`](#hisat2-build-options-bmaxdivn). Default (in terms of the [`--bmaxdivn`](#hisat2-build-options-bmaxdivn) parameter) is [`--bmaxdivn`](#hisat2-build-options-bmaxdivn) 4. This is configured automatically by default; use [`-a`/`--noauto`](#hisat2-build-options-a) to configure manually.

- `--bmaxdivn <int>`  
    The maximum number of suffixes allowed in a block, expressed as a fraction of the length of the reference. Setting this option overrides any previous setting for [`--bmax`](#hisat2-build-options-bmax), or [`--bmaxdivn`](#hisat2-build-options-bmaxdivn). Default: [`--bmaxdivn`](#hisat2-build-options-bmaxdivn) 4. This is configured automatically by default; use [`-a`/`--noauto`](#hisat2-build-options-a) to configure manually.

- `--dcv <int>`  
    Use `<int>` as the period for the difference-cover sample. A larger period yields less memory overhead, but may make suffix sorting slower, especially if repeats are present. Must be a power of 2 no greater than 4096. Default: 1024. This is configured automatically by default; use [`-a`/`--noauto`](#hisat2-build-options-a) to configure manually.

- `--nodc`  
    Disable use of the difference-cover sample. Suffix sorting becomes quadratic-time in the worst case (where the worst case is an extremely repetitive reference). Default: off.

- `-r/--noref`  
    Do not build the `NAME.3.ht2` and `NAME.4.ht2` portions of the index, which contain a bitpacked version of the reference sequences and are used for paired-end alignment.

- `-3/--justref`  
    Build only the `NAME.3.ht2` and `NAME.4.ht2` portions of the index, which contain a bitpacked version of the reference sequences and are used for paired-end alignment.

- `-o/--offrate <int>`  
    To map alignments back to positions on the reference sequences, it’s necessary to annotate (“mark”) some or all of the [Burrows-Wheeler](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform) rows with their corresponding location on the genome. [`-o`/`--offrate`](#hisat2-build-options-o) governs how many rows get marked: the indexer will mark every 2^`<int>` rows. Marking more rows makes reference-position lookups faster, but requires more memory to hold the annotations at runtime. The default is 4 (every 16th row is marked; for human genome, annotations occupy about 680 megabytes).

- `-t/--ftabchars <int>`  
    The ftab is the lookup table used to calculate an initial [Burrows-Wheeler](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform) range with respect to the first `<int>` characters of the query. A larger `<int>` yields a larger lookup table but faster query times. The ftab has size 4^(`<int>`+1) bytes. The default setting is 10 (ftab is 4MB).

- `--localoffrate <int>`  
    This option governs how many rows get marked in a local index: the indexer will mark every 2^`<int>` rows. Marking more rows makes reference-position lookups faster, but requires more memory to hold the annotations at runtime. The default is 3 (every 8th row is marked, this occupies about 16KB per local index).

- `--localftabchars <int>`  
    The local ftab is the lookup table in a local index. The default setting is 6 (ftab is 8KB per local index).

- `-p <int>`  
    Launch `NTHREADS` parallel build threads (default: 1).

- `--snp <path>`  
    Provide a list of SNPs (in the HISAT2’s own format) as follows (five columns).
    
    ```
    SNP ID <tab> snp type (single, deletion, or insertion) <tab> chromosome name <tab> zero-offset based genomic position of a SNP <tab> alternative base (single), the length of SNP (deletion), or insertion sequence (insertion)
    ```
    
    For example,
    
    ```
    rs58784443      single  13      18447947        T  
    ```
    
    Use `hisat2_extract_snps_haplotypes_UCSC.py` (in the HISAT2 package) to extract SNPs and haplotypes from a dbSNP file (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz). or `hisat2_extract_snps_haplotypes_VCF.py` to extract SNPs and haplotypes from a VCF file (e.g. ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz).
    

- `--haplotype <path>`  
    Provide a list of haplotypes (in the HISAT2’s own format) as follows (five columns).
    
    ```
    Haplotype ID <tab> chromosome name <tab> zero-offset based left coordinate of haplotype <tab> zero-offset based right coordinate of haplotype <tab> a comma separated list of SNP ids in the haplotype
    ```
    
    For example,
    
    ```
    ht35    13      18446877        18446945        rs12381094,rs12381056,rs192016659,rs538569910  
    ```
    
    See the above option, –snp, about how to extract haplotypes. This option is not required, but haplotype information can keep the index construction from exploding and reduce the index size substantially.
    

- `--ss <path>`  
    Note this option should be used with the following [`--exon`](#hisat2-build-options-exon) option. Provide a list of splice sites (in the HISAT2’s own format) as follows (four columns).
    
    ```
    chromosome name <tab> zero-offset based genomic position of the flanking base on the left side of an intron <tab> zero-offset based genomic position of the flanking base on the right <tab> strand
    ```
    
    Use `hisat2_extract_splice_sites.py` (in the HISAT2 package) to extract splice sites from a GTF file.
    

- `--exon <path>`  
    Note this option should be used with the above [`--ss`](#hisat2-build-options-ss) option. Provide a list of exons (in the HISAT2’s own format) as follows (three columns).
    
    ```
    chromosome name <tab> zero-offset based left genomic position of an exon <tab> zero-offset based right genomic position of an exon
    ```
    
    Use `hisat2_extract_exons.py` (in the HISAT2 package) to extract exons from a GTF file.
    

- `--seed <int>`  
    Use `<int>` as the seed for pseudo-random number generator.

- `--cutoff <int>`  
    Index only the first `<int>` bases of the reference sequences (cumulative across sequences) and ignore the rest.

- `-q/--quiet`  
    `hisat2-build` is verbose by default. With this option `hisat2-build` will print only error messages.

- `-h/--help`  
    Print usage information and quit.

- `--version`  
    Print version information and quit.

## The `hisat2-inspect` index inspector

`hisat2-inspect` extracts information from a HISAT2 index about what kind of index it is and what reference sequences were used to build it. When run without any options, the tool will output a FASTA file containing the sequences of the original references (with all non-`A`/`C`/`G`/`T` characters converted to `N`s). It can also be used to extract just the reference sequence names using the [`-n`/`--names`](#hisat2-inspect-options-n) option or a more verbose summary using the [`-s`/`--summary`](#hisat2-inspect-options-s) option.

## Command Line

Usage:

```
hisat2-inspect [options]* <ht2_base>
```

### Main arguments

- `<ht2_base>`  
    The basename of the index to be inspected. The basename is name of any of the index files but with the `.X.ht2` suffix omitted. `hisat2-inspect` first looks in the current directory for the index files, then in the directory specified in the `HISAT2_INDEXES` environment variable.

### Options

- `-a/--across <int>`  
    When printing FASTA output, output a newline character every `<int>` bases (default: 60).

- `-n/--names`  
    Print reference sequence names, one per line, and quit.

- `-s/--summary`  
    Print a summary that includes information about index settings, as well as the names and lengths of the input sequences. The summary has this format:
    
    ```
    Colorspace	<0 or 1>
    SA-Sample	1 in <sample>
    FTab-Chars	<chars>
    Sequence-1	<name>	<len>
    Sequence-2	<name>	<len>
    ...
    Sequence-N	<name>	<len>
    ```
    
    Fields are separated by tabs. Colorspace is always set to 0 for HISAT2.
    

- `--snp`  
    Print SNPs, and quit.

- `--ss`  
    Print splice sites, and quit.

- `--ss-all`  
    Print splice sites including those not in the global index, and quit.

- `--exon`  
    Print exons, and quit.

- `-v/--verbose`  
    Print verbose output (for debugging).
    
- `--version`  
    Print version information and quit.
    
- `-h/--help`  
    Print usage information and quit.
    

## Getting started with HISAT2

HISAT2 comes with some example files to get you started. The example files are not scientifically significant; these files will simply let you start running HISAT2 and downstream tools right away.

First follow the manual instructions to [obtain HISAT2](#obtaining-hisat2). Set the `HISAT2_HOME` environment variable to point to the new HISAT2 directory containing the `hisat2`, `hisat2-build` and `hisat2-inspect` binaries. This is important, as the `HISAT2_HOME` variable is used in the commands below to refer to that directory.

## Indexing a reference genome

To create an index for the genomic region (1 million bps from the human chromosome 22 between 20,000,000 and 20,999,999) included with HISAT2, create a new temporary directory (it doesn’t matter where), change into that directory, and run:

```
$HISAT2_HOME/hisat2-build $HISAT2_HOME/example/reference/22_20-21M.fa --snp $HISAT2_HOME/example/reference/22_20-21M.snp 22_20-21M_snp
```

The command should print many lines of output then quit. When the command completes, the current directory will contain ten new files that all start with `22_20-21M_snp` and end with `.1.ht2`, `.2.ht2`, `.3.ht2`, `.4.ht2`, `.5.ht2`, `.6.ht2`, `.7.ht2`, and `.8.ht2`. These files constitute the index - you’re done!

You can use `hisat2-build` to create an index for a set of FASTA files obtained from any source, including sites such as [UCSC](http://genome.ucsc.edu/cgi-bin/hgGateway), [NCBI](http://www.ncbi.nlm.nih.gov/sites/genome), and [Ensembl](http://www.ensembl.org/). When indexing multiple FASTA files, specify all the files using commas to separate file names. For more details on how to create an index with `hisat2-build`, see the [manual section on index building](#the-hisat2-build-indexer). You may also want to bypass this process by obtaining a pre-built index.

## Aligning example reads

Stay in the directory created in the previous step, which now contains the `22_20-21M` index files. Next, run:

```
$HISAT2_HOME/hisat2 -f -x $HISAT2_HOME/example/index/22_20-21M_snp -U $HISAT2_HOME/example/reads/reads_1.fa -S eg1.sam
```

This runs the HISAT2 aligner, which aligns a set of unpaired reads to the genome region using the index generated in the previous step. The alignment results in SAM format are written to the file `eg1.sam`, and a short alignment summary is written to the console. (Actually, the summary is written to the `"standard error"` or `"stderr"` filehandle, which is typically printed to the console.)

To see the first few lines of the SAM output, run:

You will see something like this:

```
@HD     VN:1.0  SO:unsorted
@SQ     SN:22:20000001-21000000 LN:1000000
@PG     ID:hisat2       PN:hisat2       VN:2.0.0-beta
1       0       22:20000001-21000000    397984  255     100M    *       0       0       GCCTGTGAGGGAGCCCCGGACCCGGTCAGAGCAGGAGCCTGGCCTGGGGCCAAGTTCACCTTATGGACTCTCTTCCCTGCCCTTCCAGGAGCAGCTCACT    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YT:Z:UU NH:i:1
2       16      22:20000001-21000000    398131  255     100M    *       0       0       ATGACACACTGTACACACCAGGGGCCCTGTGCTCCCCAGGAAGAGGGCCCTCACTTGAAGCGGGGCCCGATGGCCGCCACGTGCCGGTTCATGCTCCCCT    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:80A19      YT:Z:UU NH:i:1  Zs:Z:80|S|rs576159895
3       16      22:20000001-21000000    398222  255     100M    *       0       0       TGCTCCCCTTGGCCCCGCCGATGTTCAGGGACATGGAGCGCTGCAGCAGGCTGGAGAAGATCTCCACTTGGTCAGAGCTGCAGTACTTGGCGATCTCAAA    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:16A83      YT:Z:UU NH:i:1  Zs:Z:16|S|rs2629364
4       16      22:20000001-21000000    398247  255     90M200N10M      *       0       0       CAGGGACATGGAGCGCTGCAGCAGGCTGGAGAAGATCTCCACTTGGTCAGAGCTGCAGTACTTGGCGATCTCAAACCGCTGCACCAGGAAGTCGATCCAG    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YT:Z:UU XS:A:-  NH:i:1
5       16      22:20000001-21000000    398194  255     100M    *       0       0       GGCCCGATGGCCGCCACGTGCCGGTTCATGCTCCCCTTGGCCCCGCCGATGTTCAGGGACATGGAGCGCTGCAGCAGGCTGGAGAAGATCTCCACTTGGT    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:17A26A55   YT:Z:UU NH:i:1  Zs:Z:17|S|rs576159895,26|S|rs2629364
6       0       22:20000001-21000000    398069  255     100M    *       0       0       CAGGAGCAGCTCACTGAAATGTGTTCCCCGTCTACAGAAGTACCGTGATACACAGACGCCCCATGACACACTGTACACACCAGGGGCCCTGTGCTCCCCA    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YT:Z:UU NH:i:1
7       0       22:20000001-21000000    397896  255     100M    *       0       0       GTGGAGTAGATCTTCTCGCGAAGCACATTGCAGATGGTTGCATTTGGAACCACATCGGCATGCAGGAGGGACAGCCCCAGGGTCAGCAGCCTGTGAGGGA    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:31G68      YT:Z:UU NH:i:1  Zs:Z:31|S|rs562662261
8       0       22:20000001-21000000    398150  255     100M    *       0       0       AGGGGCCCTGTGCTCCCCAGGAAGAGGGCCCTCACTTGAAGCGGGGCCCGATGGCCGCCACGTGCCGGTTCATGCTCCCCTTGGCCCCGCCGATGTTCAG    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:61A26A11   YT:Z:UU NH:i:1  Zs:Z:61|S|rs576159895,26|S|rs2629364
9       16      22:20000001-21000000    398329  255     8M200N92M       *       0       0       ACCAGGAAGTCGATCCAGATGTAGTGGGGGGTCACTTCGGGGGGACAGGGTTTGGGTTGACTTGCTTCCGAGGCAGCCAGGGGGTCTGCTTCCTTTATCT    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YT:Z:UU XS:A:-  NH:i:1
10      16      22:20000001-21000000    398184  255     100M    *       0       0       CTTGAAGCGGGGCCCGATGGCCGCCACGTGCCGGTTCATGCTCCCCTTGGCCCCGCCGATGTTCAGGGACATGGAGCGCTGCAGCAGGCTGGAGAAGATC    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:27A26A45   YT:Z:UU NH:i:1  Zs:Z:27|S|rs576159895,26|S|rs2629364
```

The first few lines (beginning with `@`) are SAM header lines, and the rest of the lines are SAM alignments, one line per read or mate. See the [HISAT2 manual section on SAM output](#sam-output) and the [SAM specification](http://samtools.sourceforge.net/SAM1.pdf) for details about how to interpret the SAM file format.

## Paired-end example

To align paired-end reads included with HISAT2, stay in the same directory and run:

```
$HISAT2_HOME/hisat2 -f -x $HISAT2_HOME/example/index/22_20-21M_snp -1 $HISAT2_HOME/example/reads/reads_1.fa -2 $HISAT2_HOME/example/reads/reads_2.fa -S eg2.sam
```

This aligns a set of paired-end reads to the reference genome, with results written to the file `eg2.sam`.

[SAMtools](http://samtools.sourceforge.net/) is a collection of tools for manipulating and analyzing SAM and BAM alignment files. [BCFtools](http://samtools.sourceforge.net/mpileup.shtml) is a collection of tools for calling variants and manipulating VCF and BCF files, and it is typically distributed with [SAMtools](http://samtools.sourceforge.net/). Using these tools together allows you to get from alignments in SAM format to variant calls in VCF format. This example assumes that `samtools` and `bcftools` are installed and that the directories containing these binaries are in your [PATH environment variable](http://en.wikipedia.org/wiki/PATH_\(variable\)).

Run the paired-end example:

```
$HISAT2_HOME/hisat -f -x $HISAT2_HOME/example/index/22_20-21M_snp -1 $HISAT2_HOME/example/reads/reads_1.fa -2 $HISAT2_HOME/example/reads/reads_2.fa -S eg2.sam
```

Use `samtools view` to convert the SAM file into a BAM file. BAM is a the binary format corresponding to the SAM text format. Run:

```
samtools view -bS eg2.sam > eg2.bam
```

Use `samtools sort` to convert the BAM file to a sorted BAM file. The following command requires samtools version 1.2 or higher.

```
samtools sort eg2.bam -o eg2.sorted.bam
```

We now have a sorted BAM file called `eg2.sorted.bam`. Sorted BAM is a useful format because the alignments are (a) compressed, which is convenient for long-term storage, and (b) sorted, which is convenient for variant discovery. To generate variant calls in VCF format, run:

```
samtools mpileup -uf $HISAT2_HOME/example/reference/22_20-21M.fa eg2.sorted.bam | bcftools view -bvcg - > eg2.raw.bcf
```

Then to view the variants, run:

```
bcftools view eg2.raw.bcf
```

See the official SAMtools guide to [Calling SNPs/INDELs with SAMtools/BCFtools](http://samtools.sourceforge.net/mpileup.shtml) for more details and variations on this process.