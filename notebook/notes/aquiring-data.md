---
id: aquiring-data
aliases: []
tags: []
---

# Setting up your environment

**Virtual environments** are great, because they let you have separate environments for separate projects. This is advantageous, since one project could rely on a certain library version 3, while some other may require version 4 (or even on a Python interpreter).

## Conda

An advantage that **Conda** provides is not only for managing Python libraries, but also command line tools. This can make the tool's instalation process uniform and more generalized for users that don't work on the same systems.

I recommend installing Conda with these instructions: [docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
Essentially the difference between Miniconda and Anaconda is that with Miniconda you have to install many tools manually. Install whichever you like.

To use certain bioinformatics tools, use the **Bioconda** channel. No installation is needed, only this 3 commands that alter your `.condarc` configuration:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Finally, create an environment for this project using:

```bash
conda create --name ans
```

Answer the prompt with yes to create an environment and then **activate** the environment with `conda activate ans`. Activation of the environment is local to a terminal session, so you will have to use this command again if you open another terminal instance. To deactivate it simply run `conda deactivate`.

# Finding our data

During this course we were tasked to work on the RNA-seq dataset linked to the article titled *"Wolbachia pipientis modulates germline stem cells and gene expression associated with ubiquitination and histone lysine trimethylation to rescue fertility defects in Drosophila"*.

First thing to do is to find this article on **NCBI** or **PubMed**. This is done with a quick Google search: https://academic.oup.com/genetics/article/229/3/iyae220/7934994#510949052.

In the article, find the section *Data availability*. There will be an accession number for a **BioProject** (a BioProject is a collection of biological data related to a large-scale research effort). Create a directory where we'll store all accession numbers and save it into a file:

```bash
mkdir accession_numbers
echo "PRJNA1166928" > accession_numbers/BioProjAcc.txt
```

Search for sequencing data related to this BioProject through the command line using **NCBI Entrez Direct tool**. First install it in the environment with. Then make a query and save the output in `csv` format.

```bash
conda install bioconda::entrez-direct
```

```bash
esearch -db sra -query $(cat accession_numbers/BioProjAcc.txt) |\
efetch -format runinfo > runinfo.csv
```

- `esearch` queries the database and returns all accession numbers that match that query
- `efetch` fetches the data linked to those accession numbers
- `-db` specifies an NCBI database to search
- `-query` specifies the search query
- `|` pipes the output of the previous command as input for this command
- `-format` specifies the output format

Taking a quick look at this file, it has many different columns, while we are only interested in the SRA (sequence read archive) accession numbers. Extract the first column with `cut` on the output of `tail` that will skip the first line (the header):

```bash
tail -n +2 runinfo.csv |\
cut -d ',' -f1 > accession_numbers/SraAccList.txt
```

- `-n +2` tells `tail` to start displaying from the second line
- The `-d` flag is used to specify a delimiter and the `-f1` tells `cut` to extract the first field.

Each of us students had to choose one accession number. I choose **SRR30833097**.  Save this into a file with

```bash
echo "SRR30833097" > accession_numbers/OurAcc.txt
```

## Questions regarding research and sequencing dataset

### a) What is the aim of the study, described in the scientific article?

To study how *Wolbachia pipientis* influences **gene expression** (particulary the gene expression of genes associated with ubiquitin and histone lysine trimethylation) when rescuing the fertility defect on female hypomorph **bam gene** (*bag of marbles* - a gene involved in germ cell differentiation).

### b) Provide some info about the dataset you will work with and which sequencing technology was used.

![[Pasted image 20250421105830.png]]

**Name:**
- This dataset contains reads from ovarie tissue of infected and uninfected *Drosophiila menogaster* organisms with wild-type genotypes and genotypes with hypomorphed (functional gene with reduced function) *bam* gene.
- Since this is a [3' RNA-seq](https://en.wikipedia.org/wiki/3%27_mRNA-seq), this means that the 3' untranslated regions of polyadenylated mRNA's were sequenced instead of the whole trancript.

**Sequencing technology:**
- ILLUMINA (NextSeq 500).

**Run statistics:**
- The data comes from a single sequencing run with the mentioned technology.
- 4.5 million reads were generated.
- Total number of nucleotides across all reads is 390.2 million.

**Library:**
- Source material are transcripts.
- Selection (filtering of RNA) by polyA tail was used.
- SINGLE layout means our reads were single-ended.

### c) Which kit was used for sequencing library preparation? Does this kit preserve strand information (stranded library) or not?

The kit is not listed, so I cannot answer this question.

### d) What is the advantage of stranded mRNA library preparation compared to non-stranded library?

**Stranded libraries** are prepared in a way to contain information about the strand of cDNA from which the transcript originates.

This kind of library:
- allows you to **distinguish between the** **sense** and **antisense** strands of cDNA
- is useful for
	- identifying **antisense trancription**
	- **gene annotation** and **novel gene discovery**
	- determining the **origin of RNA** molecules in **overlapping regions**
	- **accurately** quantifying gene expression

> source: https://youtu.be/yp9A5E-Y49Y?si=qzVTqlXUrEowuIKG, https://lcsciences.com/why-is-strand-specific-library-preparation-important/

# Downloading data

Sequencing data is obtained from the SRA database with the SRA Toolkit. Install it with

```bash
conda install -c bioconda sra-tools
```

**Prefetch** is used to obtain *runs* (sequence files in compressed SRA format). The `--output-file` flag is used to use a file with a list of accession numbers as input. The command creates a directory named after the given accession number, where the downloaded files reside.

```bash
prefetch --option-file accession_numbers/OurAcc.txt
```

> [!CAUTION]
>
> `perl: error while loading shared libraries: libcrypt.so.1: cannot open shared object file: No such file or directory`
>
> Modern linux distributions have moved away from `libcrypt.so.1`  to `libxcrypt`. If you get this error, the most straightforward fix is to install the missing package, e.g.
>
> ```bash
> sudo pacman -S libxcrypt-compat
> ```
>
> Try running the prefetch command again. If the error persists, try installing through conda-forge with `conda install sra-tools` since it has top priority.

The prefetched runs can be converted into FastQ format using `fasterq-dump` (included in sra-tools), that takes the created directory as input:

```bash
fasterq-dump --skip-technical SRR30833097/
```

`--skip-technical` returns only biological reads.

Since we have only single-end sequences, it should output a single `.fastq` file in the current directory. Check that the line count matches 18149460 with:

```bash
wc -l SRR30833097.fastq
```

Remove the prefetched runs and run info since we won't need them anymore.

```bash
rm -rf SRR30833097 runinfo.csv
```

# Generating a quality report

To get a **full report** on all sequences in our dataset use the **FastQC tool**.

```bash
conda install -c bioconda fastqc
```

Make a new directory for FastQC reports.

```bash
mkdir fastqc_reports
```

Run it on the generated FastQ file.

```bash
fastqc SRR30833097.fastq --outdir fastqc_reports
```

Rename the report to something more meaningful and remove generated `zip` file.

```bash
cd fastqc_reports
mv SRR30833097_fastqc.html pre-filtering_fastqc.html
rm SRR30833097_fastqc.zip
cd ..
```

The quality report is the generated `.html` file. To view the rendered file, you can open it with your browser (e.g. `opera *.html`).

Alternativelly, you can generate the report with a **graphical user interface** of FastQC by executing `fastqc` in the terminal. This will open a new window where you can open your `.fastq` file and generate and view the report.

![[Pasted image 20250418203153.png]]

Click on `File > Open` and select your file. When it's finished, you should see the report. Click on `File > Save report` and save the report in the working directory as `pre-filtertering_fastqc.html`.

## Explaining the results

FastQC shows a summary of modules that were run on our data. On the left there is a **green tick** if the module seems entirely normal, an **orange exclamation mark** if results are slightly abnormal and a **red cross** if they are very unsusual.

It's important to analyse these results in detail, because these three symbols might not show the whole context.

### Basic Statistics

![[Pasted image 20250421121853.png]]

Nothing that hasn't been mentioned before, except the `%GC` number, which tells us the overall **GC (guanosine + citosine) content** across all sequences.

It is the simplest known HI (homology independent) metric used as a genome signiature. It can be easily calculated from sequence data alone. It displays a huge variation across genomes and is reasonably constant within a given genome

Based on the table below, it seems GC content in our sequences is valid.

![[Pasted image 20250422113052.png|400]]
*source: www.researchgate.net/figure/GC-content-of-Human-Mouse-Drosophila-melanogaster-Caenorhabditis-elegans_tbl1_5485258?__cf_chl_tk=CQJZOaB7o7m7DdJQ09Z_pGqtyr92uh63l1kqDyk4nvI-1745314185-1.0.1.1-Js8Lw9JrJD2ALRLlK.cLmPitbjmIqgkt5bpJ_4_iZuc*

> It is also worth noting that the manuals says the Basic Statistics module never raises a warning or an error.

### Per base sequence quality

![[Pasted image 20250421121906.png]]

It shows a range of quality values across all bases at each position. **Green area** are good quality calls, **orange area** are calls of reasonable quality and **red area** are calls of poor quality. The quality of calls will degrade as the run progresses, so the drop in quality at the end, as seen above, is not uncommon.

This plot can alert us to whether there were any problems occuring during sequencing. Our mean quality (blue line) is consistently in the green area, which means our per base sequence quality of good.

### Per sequence quality scores

![[Pasted image 20250421121926.png]]

This plot shows us the average quality score per read on the x-axis and the number of sequences with that average on the y-axis. Most of our reads have the highest quality scores.

###  Per base sequence content

![[Pasted image 20250421121941.png]]

This module presents **the percentage of each base position for which each of the 4 DNA bases had been called**. For RNA-seq the beginning is usually all over the place, because of adapters added during library preparation, although in our case there isn't a lot of noise in the beginning, but rather at the end, where the `%A` starts growing rapidly.

This module will fail when the percentage difference between AT and GC is greater that 20% in any position.

### Per sequence GC content

![[Pasted image 20250421121950.png]]

This plot shows the GC distribution across all sequences (red) compared to a theoretical normal distribution (blue). As I mentioned, our organisms GC content is around 40%. A shifted theoretical distribution indicates this bias, but this is not the reason this module fails.

The module will show a failure when our GC distribution doesn't follow the theoretical one. This could indicate **contamination** with another organism within the library (broad peaks) or presence of **over-expressed sequences** (sharp peaks).

### Per base N content

![[Pasted image 20250421121959.png]]

If a sequencer is unable to make a base call with sufficient confidence then it will substitute an N rather than a conventional base. This plot show that no N's were substituted in our sequences.

### Sequence Length Distribution

![[Pasted image 20250421122009.png]]

This module generates a graph showing the distribution of read sizes. Generally all reads will be the same length (prior to trimming and preprocessing).

### Sequence Duplication Levels

![[Pasted image 20250421122022.png]]

In a diverse library most sequences will occur only once in the final set. A high level of duplication is likely to indicate some kind of enrichment bias (e.g. PCR over amplification).

There are quite a few duplicate sequences in our dataset, but I am not sure how to interpret this.

### Overrepresented sequences

![[Pasted image 20250421122035.png]]

**Overrepresented sequences** are sequences that are found more frequently than statistical models predict it should occur by chance in a random distribution of sequences of that length.

This model lists all of the sequences which make up more than 0.1 % of the total, but it tracks only the **first 200,000 sequences**, therefore some overrepresented sequences might be missed. Minimum read length is **25bp** and reads longer than **75bp** are **truncated to 50bp**.

The module looks for matches in a database of common contaminants and reports hits if they're found (*possible source* column). If a single sequence is overrepresented, this could mean that it is higly biologically significant or that the library is contaminated. A hit **may** indicate some form of contamination (e.g. adapters).

There are not hits detected, which could mean that they are less common primers or that some other artifact is present.

### Adapter content

![[Pasted image 20250421122129.png]]

This module shows the mean percentage of adapter content or polyA or polyG tail. There was a very small percentage of adapters detected in our reads, which is almost neglegible but will be dealt with in the next section.

A much higher percent of polyA tails was detected. This is not untypical, since our reads came from 3' ends of mRNA molecules and the polyadenilation is a known chemical modification of pre-mRNA molecules.

# Checking quality of sequences

## Questions regarding read quality

### a) How is a fastq file composed?

FastQ files are text files that combine **FASTA formatted sequences** with their **quality scores** and is the standard format for storing the output of high-throughput sequencing instruments.

Since FastQ files bundle quality scores of sequences, we can take a quick look at the first stored sequence with

```bash
head -n 4 SRR30833097.fastq
```

```
@SRR30833097.1 NB500947:1144:HM5GMBGXH:1:11101:4455:1091 length=86
GGCGGTCGAGTGCCTCACAGTGTATCAAGGGTNGGCCACGNTCCTNACTAATNGNGGCTNNTTGCGCCATCGTCTCANGCAATGTT
+SRR30833097.1 NB500947:1144:HM5GMBGXH:1:11101:4455:1091 length=86
AAAAAEEEEEEEEEEEEEEEEEEEE<E//EEE#EEEEEE<#EEEE#E/AEEE#/#EEEE##EEE/EEAEEEEAE6EE#/EEEAAEA
```

Each entry starts with `@` and is followed by a sequence identifier and some other information about the sequence. The line directly below it shows the raw sequence. The `+` line again can contain information about the sequence and below are the quality values for each nucleotide.

### b) How can I count the number of reads in a fastq file? Describe different ways to perform that.

Using `grep` and `wc`. Grep is used to search for patterns using regular expression in a file. With `'^@'` we search for every line that starts (`^`) with `@`. It outputs every line that matches our pattern. We pipe this in `wc -l`, as we used before, to count the number of lines.

```bash
grep -e '^@' SRR30833097.fastq | wc -l
```

The AWK language is useful for pattern scanning and text proccessing. It proccesses files line by line just like `grep`. In this case we don't look for a certain identifier, but we use the fact that every sequence occupies 4 lines to our advantage.

`NR` (number of records) is a builtin awk-variable. It records lines processed so far. When reading line by line from file, this is essentially the current line number.

With `NR % 4 == 1` we calculate if the current line is the first line of 4. This is true for lines 1, 5, 9... which are lines that contain `@` and represent one sequence each.

The filtered lines are returned, which we can again pipe into `wc -l`.

```bash
awk 'NR % 4 == 1' SRR30833097.fastq | wc -l
```

Here's a Python script that works the same way as the grep command.

```Python
# lines.py

import sys
from time import time

if len(sys.argv) != 2:
	raise Exception("Usage: python3 lines.py <fastq-file>")

fname = sys.argv[1]

print("Processing file...")
start_time = time()

with open(fname, "r") as file:
    n = 0
    for line in file:
        if line[0] == "@":
            n += 1

end_time = time()
final_time = end_time - start_time
print(f"number of sequences in file: {n}")
print(f"time: {final_time}")
```

The program outputs:

```
$ python3 lines.py SRR30833097.fastq
Processing file...
number of sequences in file: 4537365
time: 6.220864295959473 sec
```

We can modify it to use the tactic we used with AWK.

```Python
# nlines.py

# --//--
with open(fname, "r") as file:
    data = file.readlines()
    n = len(data) / 4
# --//--
```

```
$ python3 nlines.py SRR30833097.fastq
Processing file...
number of sequences in file: 4537365.0
time: 4.813834190368652 sec
```

Since Python is slow by its nature, a C program will do this faster.

```c
// clines.c
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MAX_LINE_LEN 100

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: ./clines.out <fastq-file>");
    exit(1);
  }

  FILE *file = fopen(argv[1], "r");

  if (file == NULL) {
    perror("Error opening file.");
    exit(1);
  }

  int n = 0;               // line counter
  char line[MAX_LINE_LEN]; // line buffer

  // for measuring time more precisely than time_t
  struct timeval start_time;
  struct timeval end_time;
  double final_time = 0;

  printf("Processing file...\n");
  gettimeofday(&start_time, NULL);

  while (fgets(line, sizeof(line), file) != NULL) {
    if (line[0] == '@') {
      n++;
    }
  }

  gettimeofday(&end_time, NULL);
  // first term is difference in seconds
  // second term is diference in microseconds converted to seconds
  final_time = (double)(end_time.tv_sec - start_time.tv_sec) +
               (double)(end_time.tv_usec - start_time.tv_usec) / 1000000.0;

  printf("sequences in file: %d\n", n);
  printf("time: %f sec\n", final_time);

  return 0;
}
```

```
$ gcc clines.c -o clines.out && ./clines SRR30833097.fastq
Processing file...
sequences in file: 4537365
time: 1.321521 sec
```

And the other variation with line calculations.

```c
  // nclines.c
  // --//--
  if (argc != 2) {
    printf("Usage: ./nclines.out <fastq-file>");
    exit(1);
  }
  // --//--
  while (fgets(line, sizeof(line), file) != NULL) {
    n++;
  }
  n = n / 4;
  // snip
```

```
$ gcc nclines.c -o nclines.out && ./nclines SRR30833097.fastq
Processing file...
sequences in file: 4537365
time: 1.296595 sec
```

As you can tell, the C program's are fastest among all, but with it we sacrifice a lot of simplicity that comes with a BASH or AWK one-liner.

> The scripts are in `scripts` directory.
### c) What about the quality of your reads?

Modules that failed: per base sequence content, per base GC content, adapter content. The overall quality can be better, which is what we adress in the next section.

### d) [Describe your fastqc and/or multiqc and interpret the results](https://mugenomicscore.missouri.edu/PDF/FastQC_Manual.pdf)

My fastqc report was explained above in the [[faks/ANS/notebook/README#Checking quality of sequences|Checking quality of sequences]] section.
