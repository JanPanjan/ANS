>[!TODO]
>- [ ] document tools' versions
>- [ ] multiqc report od vseh

---

# General instructions

Instructions for filling in the worksheets:

- Worksheets include questions (some are for refreshing the knowledge, whereas some will probably require a little bit of literature review) and practical work, where you carefully insert all codes and as much as possible comments.
    
- Each step of the bioinformatics analysis should be well documented and explained why it was performed and why each parameter/option of the command was used.
    
- Bioinformatics tools which you will use have several options. I encourage you to explore why they are used for. Make sure that you add this to the report as well.
    
- Comments regarding the worksheets, how can be improved, what should be added, etc. are welcome.

---

# Setting up your environment

**Virtual environments** are great, because they let you have separate environments for separate proejcts. This is advantageous, since one project could rely on a certain package version 3, while some other may require version 4.

## Conda

An advantage that **Conda** provides is not only for managing Python libraries, but also command line tools. This can make the tool's instalation process uniform and more generalized for users that don't work on the same systems.

I recommend installing Conda with these instructions: [docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
Essentially the difference between Miniconda and Anaconda is that with Miniconda you have to install many tools manually. Install whichever you like.

### Bioconda

To use certain bioinformatics tools, we need to use the **Bioconda** channel. No installation is needed, only this 3 commands that alter your `.condarc` configuration:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## Virtual environment

Finally, we can create an environment for this project using:

```bash
conda create --name ans
```

Answer the prompt with yes to create an environment and then **activate** the environment with `conda activate ans`. Activation of the environment is local to a terminal session, so you will have to use this command again if you open another terminal instance.
To deactivate it simply run `conda deactivate`.

# Finding our data

During this course we were tasked to work on the RNA-seq dataset linked to the article titled *"Wolbachia pipientis modulates germline stem cells and gene expression associated with ubiquitination and histone lysine trimethylation to rescue fertility defects in Drosophila"*.

## BioProject accession number

First thing to do is to find this article on **NCBI** or **PubMed**. This can be done with a quick Google search: https://academic.oup.com/genetics/article/229/3/iyae220/7934994#510949052. 

In the article, find the section *Data availability*. There will be an accession number for a **BioProject** (a BioProject is a collection of biological data related to a large-scale research effort).

## SRA accession number

We can search for sequencing data related to this BioProject through the command line using **NCBI Entrez Direct tool**. First we install it in our environment with `conda install bioconda::entrez-direct`. Then we make a query and save the output in `csv` format.

- `esearch` queries the database and returns all accession numbers that match that query
- `efetch` fetches the data linked to those accession numbers

```bash
esearch -db sra -query PRJNA1166928 | efetch -format runinfo > runinfo.csv
```

- `-db` specifies an NCBI database to search
- `-query` specifies the search query
- `|` pipes the output of the previous command as input for this command
- `-format` specifies the output format

If we take a quick look at this file, it has many different columns, while we are only interested in the SRA (sequence read archive) accession numbers. We can extract the first column with `cut` on the output of `tail` that will skip the first line (the header):

```bash
tail -n +2 runinfo.csv | cut -d ',' -f1 > SraAccList.txt
```

- `-n +2` tells `tail` to start displaying from the second line 
- The `-d` flag is used to specify a delimiter and the `-f1` tells `cut` to extract the first field.

Each of us students had to choose one accession number. I choose **SRR30833097**.  Save this into a file with `echo "SRR30833097" > OurAcc.txt`.

## Questions regarding research and sequencing dataset

### a) What is the aim of the study, described in the scientific article?

To study how *Wolbachia pipientis* influences **gene expression** (particulary the gene expression of genes associated with ubiquitin and histone lysine trimethylation) when rescuing the fertility defect on female hypomorph **bam gene** (*bag of marbles* - a gene involved in germ cell differentiation).

### b) Provide some info about the dataset you will work with and which sequencing technology was used.

![[Pasted image 20250421105830.png]]

---

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

## Prefetch

Sequencing data is obtained from the SRA database with the SRA Toolkit. It can be installed with `conda install -c bioconda sra-tools`. 

```bash
prefetch --option-file OurAcc.txt
```

**Prefetch** is used to obtain *Runs* (sequence files in compressed SRA format). The `--output-file` flag is used to use a file with a list of accession numbers as input. The command creates a directory named after the given accession number, where the downloaded files reside.

The prefetched runs can be converted into FastQ format using `fasterq-dump`, that takes the created directory as input:

```bash
fasterq-dump --skip-technical SRR30833097/
```

`--skip-technical` returns only biological reads.

Since we have only single-end sequences, it should output a single `.fastq` file in the current directory. You can check that the line count matches 18149460 with:

```bash
wc -l SRR30833097.fastq
```

# Generating a quality report

## a) Run FastQC over your dataset. Explain the results

### Command line FastQC

To get a **full report** on all sequences in our dataset, we can use **FastQC tool**. Install it with `conda install -c bioconda fastqc` and run it on the `.fastq` file. 

```bash
fastqc SRR30833097.fastq
```

The quality report is the generated `.html` file. To view the rendered file, you can open it with your browser (e.g. `opera *.html`). 

### GUI FastQC

Alternativelly, you can open a **graphical user interface** of FastQC tool by executing only `fastqc`. This will open a new window where you can open your `.fastq` file and generate and view the report.

![[Pasted image 20250418203153.png]]

Click on `File > Open` and select your file. When it's finished, you should see the report. Click on `File > Save report` and save the report in the working directory as `pre-filtertering_fastqc.html`.

### Explaining the results

FastQC shows a summary of modules that were run on our data. On the left there is a **green tick** if the module seems entirely normal, an **orange exclamation mark** if results are slightly abnormal and a **red cross** if they are very unsusual.

It's important to analyse these results in detail, because these three symbols might not show the whole context.

>[!warning]
>mislm da mas pngs od grafov v neki mapci v root, dej jih namesto teh screenshotov
#### Basic Statistics

![[Pasted image 20250421121853.png]]

Nothing that hasn't been mentioned before, except the `%GC` number, which tells us the overall **GC (guanosine + citosine) content** across all sequences.

It is the simplest known HI (homology independent) metric used as a genome signiature. It can be easily calculated from sequence data alone. It displays a huge variation across genomes and is reasonably constant within a given genome

Based on the table below, it seems GC content in our sequences is valid.

![[Pasted image 20250422113052.png|400]]
*source: www.researchgate.net/figure/GC-content-of-Human-Mouse-Drosophila-melanogaster-Caenorhabditis-elegans_tbl1_5485258?__cf_chl_tk=CQJZOaB7o7m7DdJQ09Z_pGqtyr92uh63l1kqDyk4nvI-1745314185-1.0.1.1-Js8Lw9JrJD2ALRLlK.cLmPitbjmIqgkt5bpJ_4_iZuc*

> It is also worth noting that the manuals says the Basic Statistics module never raises a warning or an error.
#### Per base sequence quality

![[Pasted image 20250421121906.png]]

It shows a range of quality values across all bases at each position. **Green area** are good quality calls, **orange area** are calls of reasonable quality and **red area** are calls of poor quality. The quality of calls will degrade as the run progresses, so the drop in quality at the end, as seen above, is not uncommon. 

This plot can alert us to whether there were any problems occuring during sequencing. Our mean quality (blue line) is consistently in the green area, which means our per base sequence quality of good.

#### Per sequence quality scores

![[Pasted image 20250421121926.png]]

This plot shows us the average quality score per read on the x-axis and the number of sequences with that average on the y-axis. Most of our reads have the highest quality scores.

####  Per base sequence content

![[Pasted image 20250421121941.png]]

This module presents **the proportion of each base position for which each of the 4 DNA bases had been called**. For RNA-seq the beginning is usually all over the place, because of adapters added during library preparation, although in our case there isn't a lot of noise in the beginning, but rather at the end, where the `%A` starts growing rapidly.

This module will fail when the percentage difference between AT and GC is greater that 20% in any position.

#### Per sequence GC content

![[Pasted image 20250421121950.png]]

This plot shows the GC distribution across all sequences (red) compared to a theoretical normal distribution (blue). As I mentioned, our organisms GC content is around 40%. A shifted theoretical distribution indicates this bias, but this is not the reason this module fails.

The module will show a failure when our GC distribution doesn't follow the theoretical one. This could indicate **contamination** with another organism within the library (broad peaks) or presence of **over-expressed sequences** (sharp peaks).

#### Per base N content

![[Pasted image 20250421121959.png]]

If a sequencer is unable to make a base call with sufficient confidence then it will substitute an N rather than a conventional base. This plot show that no N's were substituted in our sequences.

#### Sequence Length Distribution

![[Pasted image 20250421122009.png]]

This module generates a graph showing the distribution of read sizes. Generally all reads will be the same length (prior to trimming and preprocessing).

#### Sequence Duplication Levels

![[Pasted image 20250421122022.png]]

In a diverse library most sequences will occur only once in the final set. A high level of duplication is likely to indicate some kind of enrichment bias (e.g. PCR over amplification).

There are quite a few duplicate sequences in our dataset, but I am not sure how to interpret this.

#### Overrepresented sequences

![[Pasted image 20250421122035.png]]

**Overrepresented sequences** are sequences that are found more frequently than statistical models predict it should occur by chance in a random distribution of sequences of that length.

This model lists all of the sequences which make up more than 0.1 % of the total, but it tracks only the **first 200,000 sequences**, therefore some overrepresented sequences might be missed. Minimum read length is **25bp** and reads longer than **75bp** are **truncated to 50bp**.

The module looks for matches in a database of common contaminants and reports hits if they're found (*possible source* column). If a single sequence is overrepresented, this could mean that it is higly biologically significant or that the library is contaminated. A hit **may** indicate some form of contamination (e.g. adapters). 

There are not hits detected, which could mean that they are less common primers or that some other artifact is present.

#### Adapter content

![[Pasted image 20250421122129.png]]

This module shows the mean percentage of adapter content or polyA or polyG tail. There was a very small percentage of adapters detected in our reads, which is almost neglegible but will be dealt with in the next section.

A much higher percent of polyA tails was detected. This is not untypical, since our reads came from 3' ends of mRNA molecules and the polyadenilation is a known chemical modification of pre-mRNA molecules.

## b) Exchange FastQC results with your colleagues and run MultiQC to get a joined report for all datasets

# Checking quality of sequences

## Questions

### a) How is a fastq file composed?

FastQ files are text files that combine **FASTA formatted sequences** with their **quality scores** and is the standard format for storing the output of high-throughput sequencing instruments.

Since FastQ files bundle quality scores of sequences, we can take a quick look at the first stored sequence with `head -n 4 SRR30833097.fastq`. This will output

```
@SRR30833097.1 NB500947:1144:HM5GMBGXH:1:11101:4455:1091 length=86
GGCGGTCGAGTGCCTCACAGTGTATCAAGGGTNGGCCACGNTCCTNACTAATNGNGGCTNNTTGCGCCATCGTCTCANGCAATGTT
+SRR30833097.1 NB500947:1144:HM5GMBGXH:1:11101:4455:1091 length=86
AAAAAEEEEEEEEEEEEEEEEEEEE<E//EEE#EEEEEE<#EEEE#E/AEEE#/#EEEE##EEE/EEAEEEEAE6EE#/EEEAAEA
```

Each entry starts with `@` and is followed by a sequence identifier and some other information about the sequence. The line directly below it shows the raw sequence. The `+` line again can contain information about the sequence and below are the quality values for each nucleotide.

### b) How can I count the number of reads in a fastq file? Describe different ways to perform that.

#### BASH

Using `grep` and `wc`. Grep is used to search for patterns using regular expression in a file. With `'^@'` we search for every line that starts (`^`) with `@`. It outputs every line that matches our pattern. We pipe this in `wc -l`, as we used before, to count the number of lines. 

```bash
grep -e '^@' SRR30833097.fastq | wc -l
```

#### AWK

The AWK language is useful for pattern scanning and text proccessing. It proccesses files line by line just like `grep`. In this case we don't look for a certain identifier, but we use the fact that every sequence occupies 4 lines to our advantage.

`NR` (number of records) is a builtin awk-variable. It records lines processed so far. When reading line by line from file, this is essentially the current line number.

With `NR % 4 == 1` we calculate if the current line is the first line of 4. This is true for lines 1, 5, 9... which are lines that contain `@` and represent one sequence each.

The filtered lines are returned, which we can again pipe into `wc -l`.

```bash
awk 'NR % 4 == 1' SRR30833097.fastq | wc -l
```

#### Python

I made two scripts, one that uses pattern matching like the BASH command and one that uses calculating the mod of line numbers. They are a little more verbose than the other two options, but still pretty straightforward. 

I also measured time needed to process the file to determine which approach is faster.

##### Pattern matching

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

##### Line count calculation 

These are the only altered lines.

```Python
# nlines.py
# snip

with open(fname, "r") as file:
    data = file.readlines()
    n = len(data) / 4
	
# snip
```

```
$ python3 nlines.py SRR30833097.fastq 
Processing file...
number of sequences in file: 4537365.0
time: 4.813834190368652 sec
```

The second approach was faster by about 1.4 seconds, which is a difference, but it's not huge, since we are already using a file with 18 million lines. 
The first one may be more robust since it doesn't make assumptions about line content, while the second one *assumes* we have lines that can be evenly divided by 4.

#### C instead of Python?

Since Python is slow by its nature, maybe a simple C program will do this faster. 

##### Pattern matching

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

When we compile and run this program, we get this:

```
$ gcc clines.c -o clines.out && ./clines SRR30833097.fastq 
Processing file...
sequences in file: 4537365
time: 1.321521 sec
```

Now this is a great improvement from the Python scripts. Will the second approach with calculating the mod of lines be faster?

##### Line count calculation 

These are the only altered lines.

```c
  // nclines.c
  // snip
  if (argc != 2) {
    printf("Usage: ./nclines.out <fastq-file>");
    exit(1);
  }
  // snip
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

It's faster by a fraction. Again, not a big difference, but faster nevertheless, since we skip the comparisons for every line and do only one operation more after we read the whole file.

### c) What about the quality of your reads?

Modules that failed:

- per base sequence content
- per base GC content
- adapter content

The overall quality can be better, which is what we adress in the next section.

### d) [Describe your fastqc and/or multiqc and interpret the results](https://mugenomicscore.missouri.edu/PDF/FastQC_Manual.pdf)

My fastqc report was explained above in the [[practical-notebook#Checking quality of sequences|Checking quality of sequences]] section.

Multiqc report...

# Filtering and trimming based on quality parameters  

This is the **preprocessing** step, where we try and make the quality of our reads better using [fastp](https://github.com/OpenGene/fastp), an ultra-fast all-in-one FASTQ preprocessor. based on the previously generated FastQC report.

First install fastp with conda.

```bash
conda install -c bioconda fastp
```

The program filters our reads and generates a `fastq` file and a `html` report. We can try running the tool without any flags and check the quality of processed reads.

## Minimal trimming

```bash
fastp
-i SRR30833097.fastq
-o fastp-out-minimal.fq
--html fastp-report-minimal.html
--dont_eval_duplication # save time and memory
```

- `-i` specifies the input fastq file
- `-o` specifies name of the output fastq file
- `--html` specifies name for the output fastp report
- `--dont_eval_duplication` to save time and memory by not evaluating duplication rate of reads

The program outputs some text in the terminal. If we inspect it, we can inspect quality of reads before and after filtering.

```
...
No adapter detected for read1

Read1 before filtering:
total reads: 4537365
total bases: 390213390
Q20 bases: 359416837(92.1078%)
Q30 bases: 349585493(89.5883%)

Read1 after filtering:
total reads: 4445266
total bases: 382292876
Q20 bases: 355260023(92.9288%)
Q30 bases: 345871060(90.4728%)

Filtering result:
reads passed filter: 4445266
reads failed due to low quality: 90121
reads failed due to too many N: 1978
reads failed due to too short: 0
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
...
```

- no adapters were detected, which isn't very good, since FastQC detected some, but it will be dealt with in the next steps
- after filtering quality went up for almost 1%
- around 90k reads had too low quality
- around 2k reads were too long

So about 92k reads were removed and our quality increased only by almost 1%. Let's run fastqc again on these filtered reads, to see if they pass.

> If you missed the program's output in the terminal, you can open the generated HTML in the browser.

Runing fastqc on newly generated fastq file as before:

![[Pasted image 20250426172831.png]]

We can see that the polyA tail situation hasn't improved. This is why we need to add some flags to our `fastp` program, so it can do a better job at filtering low quality reads.

## Trim polyX tails

We can try with `--trim_poly_x`, which should (according to the documentation) trim the tails of our reads, where `x` represents any base.

```bash
fastp
-i SRR30833097.fastq
-o fastp-out-trimpolyx.fq
--trim_poly_x
--html fastp-report-trimpolyx.html
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4446923
total bases: 370914387
Q20 bases: 344936270(92.9962%)
Q30 bases: 335932934(90.5689%)

Filtering result:
reads passed filter: 4446923
reads failed due to low quality: 86186
reads failed due to too many N: 1944
reads failed due to too short: 2312
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 716650
bases trimmed in polyX tail: 12111332
...
```

Quality improved for almost 0.1%, which is not great progress. Many reads were detected with polX tail at the 3' end.

![[Pasted image 20250426174128.png]]

A little better, but not good enough.

## Trim polyG and polyX tails

PolyX tail trimming is disabled by default, but is similar to polyG tail trimming. PolyG tails can happen with Illumina NextSeq (which is the machine that was used to sequence our reads), since `G` means no signal in the Illumina two-color systems. If both trimming options are enabled, fastp will first do polyG trimming and then polyX.

Since we don't have a great amount (or any, its hard to see) of reads with a polyG tail, I don't think this will improve the quality greatly, but it doesn't make it worse, so I think it's fine to add to our command.

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out-trimpolygx.fq 
--trim_poly_g 
--trim_poly_x 
--html fastp-report-trimpolygx.html 
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4443750
total bases: 370455963
Q20 bases: 344591801(93.0183%)
Q30 bases: 335614022(90.5948%)

Filtering result:
reads passed filter: 4443750
reads failed due to low quality: 85378
reads failed due to too many N: 1943
reads failed due to too short: 6294
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 669709
bases trimmed in polyX tail: 11402404
...
```

There were less reads with polyX tails detected, probably because the polyG trimming was done first.

![[Pasted image 20250426174210.png]]

## Adjust minimum polyX tail length

We still have to add some options to our command. We can try adjusting minimum polyX tail length used for detecting. Default is 10. We can try with 5.

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out-trimpolygx5.fq 
--trim_poly_g 
--trim_poly_x 
--poly_x_min_len 5 
--html fastp-report-trimpolygx5.html 
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4445653
total bases: 370983902
Q20 bases: 345103461(93.0238%)
Q30 bases: 336143486(90.6086%)

Filtering result:
reads passed filter: 4445653
reads failed due to low quality: 84287
reads failed due to too many N: 1948
reads failed due to too short: 5477
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 897389
bases trimmed in polyX tail: 10973224
...
```

The quality increased by about 0.1% which is something, but not much. 

- Less reads failed due to being too short.
- A lot more reads were detected with polyX tail, which seems logical since we lowered the min tail size, but less were trimmed, which seems paradoxical.

The FastQC report shows that truly less reads were trimmed and consequently the polyA tail percentage is higher. The following picture shows two overlapping images.

![[Pasted image 20250426181750.png]]

---

We can also try to make the minimum length bigger, but I doubt it will benefit the quality. Let's make it 15 and check the results.

```
...
Read1 after filtering:
total reads: 4442425
total bases: 373018357
Q20 bases: 346884901(92.9941%)
Q30 bases: 337799434(90.5584%)

Filtering result:
reads passed filter: 4442425
reads failed due to low quality: 86689
reads failed due to too many N: 1957
reads failed due to too short: 6294
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 401716
bases trimmed in polyX tail: 8684854
...
```

As expected, even less reads were trimmed and the quality didn't improve.

## Use a sliding window

We have to look for alternative options for tail trimming. We can use the `--cut_tail` flag. This will create a sliding window from 3' to 5' end and drop bases in the window, if their mean quality drops below some threshold, otherwise it will stop trimming.

The window size is set using `--cut_window_size` and is 4 by default. The threshold is set using `--cut_tail_mean_quality` and is 20 by default, must be between 1 and 36.

Let's see what it does without changing these parameters.

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out-trimpolygxwindow.fq 
--trim_poly_g 
--trim_poly_x 
--cut_tail 
--html fastp-report-trimpolygxwindow.html 
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4458454
total bases: 359655463
Q20 bases: 338369588(94.0816%)
Q30 bases: 329735098(91.6808%)

Filtering result:
reads passed filter: 4458454
reads failed due to low quality: 62919
reads failed due to too many N: 1833
reads failed due to too short: 14159
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 969006
bases trimmed in polyX tail: 17650157
...
```

The quality increased more than 1%, which means we are on the right track. A lot more reads were trimmed, but less failed due to low quality. This could mean that it first cuts low quality bases in the window.

> The documentation has a "features" section where each feature is numbered. This option is before polyG and X tail trimming, which could mean that it does this first.

What does FastQC show us?

![[Pasted image 20250426183302.png]]

The line dropped again, but it's still just a little bit over the modules threshold for failure, which is 10%. We can play around with the parameters to get a better quality.

---

First, let's try to make the window a bigger size. It says that it stops if the mean quality is more than the threshold, which could in theory mean that it stops immeadiately on the first 4 bases. Let's set the window size to 10 with `--cut_window_size 10`.

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out-trimpolygxwindow20.fq 
--trim_poly_g 
--trim_poly_x 
--cut_tail 
--cut_window_size 10 
--html fastp-report-trimpolygxwindow20.html 
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4463315
total bases: 358814939
Q20 bases: 338128014(94.2347%)
Q30 bases: 329597651(91.8573%)

Filtering result:
reads passed filter: 4463315
reads failed due to low quality: 52540
reads failed due to too many N: 1749
reads failed due to too short: 19761
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 981280
bases trimmed in polyX tail: 17829600
...
```

The quality increased again, but not by much, even though a lot of reads failed due to low quality. Our reads are not very long, so we shouldn't increase this too much, but let's try with size 20.

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out-trimpolygxwindow20.fq 
--trim_poly_g 
--trim_poly_x 
--cut_tail 
--cut_window_size 20 
--html fastp-report-trimpolygxwindow20.html 
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4469047
total bases: 362654545
Q20 bases: 340736791(93.9563%)
Q30 bases: 332102803(91.5755%)

Filtering result:
reads passed filter: 4469047
reads failed due to low quality: 41519
reads failed due to too many N: 1675
reads failed due to too short: 25124
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 877220
bases trimmed in polyX tail: 15206635
...
```

The quality dropped again, since less tails were trimmed. Since we are scanning a bigger window, good quality bases can be included in the mean quality calculation and consequently low quality tails are saved from trimming. 

## Increase the threshold in the sliding window

Let's keep window size at 10 and increase the mean quality threshold from default 20 to max 36. This could potentially eliminate too many bases, but it would trim all reads with low quality tails.

>[!warning]
>Although it says that max is 36 on the github page, the command will fail with this value and print:
>```bash
> ERROR: the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.
>```
>So we will set it to 30.

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out-trimpolygxwindow10meanquality30.fq 
--trim_poly_g 
--trim_poly_x 
--cut_tail 
--cut_window_size 10 
--cut_mean_quality 30 
--html fastp-report-trimpolygxwindow10meanquality30.html 
--dont_eval_duplication
```

```
...
Read1 after filtering:
total reads: 4473488
total bases: 335974371
Q20 bases: 322232306(95.9098%)
Q30 bases: 315263902(93.8357%)

Filtering result:
reads passed filter: 4473488
reads failed due to low quality: 13029
reads failed due to too many N: 1065
reads failed due to too short: 49783
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
reads with polyX in 3' end: 1345574
bases trimmed in polyX tail: 19454923
...
```

Now this is real progress. The quality increased by about 2%. It trimmed about 1.5M more reads. A lot more reads failed due to being too short, which I don't know how to interpret.

The small adapter content present also went down, which is great.

![[Pasted image 20250426191615.png]]

>[!Note] *per base sequence content* and *per sequence GC content*
>These two modules still fail, but if we think logically, the first one can be dependent on the other. 
>Since our reads have skewed distribution of bases, the second one fails and with it the first one fails too, but this isn't critical in my opinion.

## Final command

Now we can remove all generated reports and filtered copies that we won0t need anymore. You're free to keep them, but I think our working directory will be more organized if we keep only one copy.

Remove all files whose names are matched with the pattern `fastp*` 

```bash
rm fastp*
```

The final command used to filter our fastq data is 

```bash
fastp 
-i SRR30833097.fastq 
-o fastp-out.fq 
--trim_poly_g 
--trim_poly_x 
--cut_tail 
--cut_window_size 10 
--cut_mean_quality 30 
--html fastp-report.html 
--dont_eval_duplication
```

# Downloading the reference genome and an annotation file

[[downloading reference genomes and annotation files]]

# Alignment using HISAT2

[[alignement using hisat2]]

# Quantifying the expression of transcripts using RNA- seq data with Salmon

[[quantifying the expression of transcripts using RNA-seq data with salmon]]