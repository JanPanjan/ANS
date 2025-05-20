# General instructions

Instructions for filling in the worksheets:

- Worksheets include questions (some are for refreshing the knowledge, whereas some will probably require a little bit of literature review) and practical work, where you carefully insert all codes and as much as possible comments.

- Each step of the bioinformatics analysis should be well documented and explained why it was performed and why each parameter/option of the command was used.

- Bioinformatics tools which you will use have several options. I encourage you to explore why they are used for. Make sure that you add this to the report as well.

- Comments regarding the worksheets, how can be improved, what should be added, etc. are welcome.

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

> [!bug] perl: error while loading shared libraries: libcrypt.so.1: cannot open shared object file: No such file or directory
>
> Modern linux distributions have moved away from `libcrypt.so.1`  to `libxcrypt`. If you get this error, the most straightforward fix is to install the missing package, e.g.
> 
> ```bash
> sudo pacman -S libxcrypt-compat
> ```
>
> Try running the prefetch command again. If the error persists, try installing through conda-forge with `conda install sra-tools` since it has top priority.

The prefetched runs can be converted into FastQ format using `fasterq-dump`, that takes the created directory as input:

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

My fastqc report was explained above in the [[practical-notebook#Checking quality of sequences|Checking quality of sequences]] section.

# Filtering and trimming based on quality parameters

This is the **preprocessing** step, where we try and make the quality of our reads better using [fastp](https://github.com/OpenGene/fastp), an ultra-fast all-in-one FASTQ preprocessor. based on the previously generated FastQC report.

First install it with conda.

```bash
conda install -c bioconda fastp
```

Create a directory for FASTP results.

```bash
mkdir fastp_results
cd fast_results
```

The program filters our reads and generates a `fastq` file and a `html` report. Try running the tool without any flags and check the quality of processed reads.

> Go to the end of this section to view the final command used.

## Minimal trimming

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-minimal.fq \
--html fastp-report-minimal.html
--dont_eval_duplication \
```

- `-i` specifies the input fastq file
- `-o` specifies name of the output fastq file
- `--html` specifies name for the output fastp report
- `--dont_eval_duplication` to save time and memory by not evaluating duplication rate of reads

The program outputs some text in the terminal. 

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

So about 92k reads were removed and our quality increased only by almost 1%. Run fastqc again on these filtered reads, to see if they pass.

> If you missed the program's output in the terminal, you can open the generated HTML in the browser.

Runing fastqc on newly generated fastq file as before:

```bash
fastqc fastp-out-minimal.fq --out-dir ../fastqc_reports
```

![[Pasted image 20250426172831.png]]

We see that the polyA tail situation hasn't improved. This is why we  add some flags to our `fastp` program, so it can do a better job at filtering low quality reads.

## Trim polyX tails

Try with `--trim_poly_x`, which should (according to the documentation) trim the tails of our reads, where `x` represents any base.

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolyx.fq \
--trim_poly_x \
--html fastp-report-trimpolyx.html \
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
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolygx.fq \
--trim_poly_g \
--trim_poly_x \
--html fastp-report-trimpolygx.html
--dont_eval_duplication \
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

We still have to add some options to our command. Try adjusting minimum polyX tail length used for detecting. Default is 10. We can try with 5.

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolygx5.fq \
--trim_poly_g \
--trim_poly_x \
--poly_x_min_len 5
--html fastp-report-trimpolygx5.html \
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

Try to make the minimum length bigger, but I doubt it will benefit the quality. Let's make it 15 and check the results.

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

We have to look for alternative options for tail trimming. Use the `--cut_tail` flag. This will create a sliding window from 3' to 5' end and drop bases in the window, if their mean quality drops below some threshold, otherwise it will stop trimming.

The window size is set using `--cut_window_size` and is 4 by default. The threshold is set using `--cut_tail_mean_quality` and is 20 by default, must be between 1 and 36.

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolygxwindow.fq \
--trim_poly_g \
--trim_poly_x \
--cut_tail \
--html fastp-report-trimpolygxwindow.html \
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

Let's make the window a bigger size. It says that it stops if the mean quality is more than the threshold, which could in theory mean that it stops immeadiately on the first 4 bases. Set the window size to 10 with `--cut_window_size 10`.

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolygxwindow20.fq \
--trim_poly_g \
--trim_poly_x \
--cut_tail \
--cut_window_size 10 \
--html fastp-report-trimpolygxwindow20.html
--dont_eval_duplication \
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

The quality increased again, but not by much, even though a lot of reads failed due to low quality. Our reads are not very long, so we shouldn't increase this too much, but try with size 20.

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolygxwindow20.fq \
--trim_poly_g \
--trim_poly_x \
--cut_tail \
--cut_window_size 20 \
--html fastp-report-trimpolygxwindow20.html
--dont_eval_duplication \
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

Keep window size at 10 and increase the mean quality threshold from default 20 to max 36. This could potentially eliminate too many bases, but it would trim all reads with low quality tails.

>[!warning]
>Although it says that max is 36 on the github page, the command will fail with this value and print:
>```bash
> ERROR: the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.
>```
>So we will set it to 30.

```bash
fastp \
-i ../SRR30833097.fastq \
-o fastp-out-trimpolygxwindow10meanquality30.fq \
--trim_poly_g \
--trim_poly_x \
--cut_tail \
--cut_window_size 10 \
--cut_mean_quality 30 \
--html fastp-report-trimpolygxwindow10meanquality30.html \
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

Remove all generated reports and filtered copies that you won't need anymore. Our working directory will be more organized if we keep only one copy.

```bash
rm fastp* # careful with this if you renamed something to have this prefix during this whole project
```

The final command used to filter the fastq data is

```bash
fastp \
-i ../SRR30833097.fastq \
-o filtered-reads.fq \
--trim_poly_g \
--trim_poly_x \
--cut_tail \
--cut_window_size 10 \
--cut_mean_quality 30 \
--html filtered-reads.html \
--dont_eval_duplication
```

Run FastQC again on filtered reads.

```bash
fastqc filtered-reads.fq --outdir ../fastqc_reports
```

Rename the report.

```bash
cd .. # make sure this puts you in root directory
cd fastqc_reports
mv filtered-reads_fastqc.html post-filtering_fastqc.html
rm *.zip
cd ..
```

# Downloading the reference genome and annotation file

The reference genome and its annotation can be downloaded using NCBI command line tools. Install them with

```bash
conda install -c conda-forge ncbi-datasets-cli
```

To find its accession number, got to [NCBI](https://www.ncbi.nlm.nih.gov/) and search *drosophila melanogaster*. Click on the genome section:

![[Pasted image 20250517181538.png|400]]

Copy and echo the RefSeq/GenBank accession number into a file:

```bash
echo "GCF_000001215.4" > accession_numbers/RefGenAcc.txt
```

To download the reference genome execute the following command:

```bash
datasets download genome accession $(cat accession_numbers/RefGenAcc.txt) \
--include genome,rna,protein,cds,gff3,gtf \
--filename genome.zip
```

- `genome` : since we want to download a genome
- `accession` : is clear
- `--include` : to specify what data files to include in the download
- `--filename` : desired name of the zip file

Unzip the downloaded file, rename it to something else and remove the zip:

> Make sure `unzip` is installed on your system.

```bash
unzip genome.zip # unzips into "ncbi_dataset/" directory 
mv ncbi_dataset ref_genome
rm genome.zip md5sum.txt README.md
```

Remove unneccessary files and move reference genome files into `ref_genome`. 

```bash
mv ref_genome/data/GCF_000001215.4/* ref_genome/
rm -rf ref_genome/data/
```

Make sure to run `ls -l ref_genome/` after. It should look like:

```
total 671356
-rw------- 1 X X  71855323 May 20 13:25 cds_from_genomic.fna
-rw------- 1 X X 145657746 May 20 13:25 GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
-rw------- 1 X X 164988613 May 20 13:26 genomic.gff
-rw------- 1 X X 185199259 May 20 13:26 genomic.gtf
-rw------- 1 X X  22970738 May 20 13:26 protein.faa
-rw------- 1 X X  96785880 May 20 13:26 rna.fna
```

> `X` will be your username.

Rename the genome assembly to something more readable.

```bash
mv ref_genome/GCF* ref_genome/GCF_000001215.4.fna
```

We downloaded coding DNA sequences, a RefSeq genome assembly, protein sequences, RNA sequences and annonation files (`gff` and `gtf`).

> RefSeq genome assemblies are NCBI-derived and mantained copies of GenBank assemblies (GCA), that include annotation, while the GCA might not.

## Questions regarding genome annotation files

### a) Describe the GTF

GTF is short for Gene Transfer Format. It's a tab-separated format. It's based on the general feature format (GFF) but contains additional gene information. It holds information about gene structure. Version 2 of GFF is identical to GTF.

There are 9 attributes with each row being a feature. Empty attributes should be donoted with a dot.

| column    | explained                                                                                                                                           |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| seqname   | name of the chromosome/scaffold. Can be given with a name or an ID (acc. number) and should match sequence names in the corresponding sequence file |
| source    | name of the program that generated this feature or the data source (database or project name).                                                      |
| feature   | name of the feature (e.g. gene, transcript, exon...).                                                                                               |
| start     | start position of the feature (base pair).                                                                                                          |
| end       | end position of the feature (base pair).                                                                                                            |
| score     | floating point value (no scores in our gtf).                                                                                                        |
| strand    | `+` defines a forward and `-` a reverse strand.                                                                                                     |
| frame     | one of 0, 1 or 2. Indicates the first base of the feature is the first base of a codon, 1 that the second base is the first base of a codon...      |
| attribute | a semicolon-separated list of tag-value pairs giving additional information about a feature.                                                        |

Example from our GTF:

| seqname     | source | feature    | start  | end    | score | strand | frame | attribute                                                  |
| ----------- | ------ | ---------- | ------ | ------ | ----- | ------ | ----- | ---------------------------------------------------------- |
| NC_004354.4 | RefSeq | gene       | 122493 | 122706 | .     | +      | .     | gene_id "Dmel_CR40469"; transcript_id ""; cyt_map "1A1-... |
| NC_004354.4 | RefSeq | transcript | 122493 | 122706 | .     | +      | .     | gene_id "Dmel_CR40469"; transcript_id "NR_003723.2"; db... |
| NC_004354.4 | RefSeq | exon       | 122493 | 122706 | .     | +      | .     | gene_id "Dmel_CR40469"; transcript_id "NR_003723.2"; db... |

### b) Examine GTF files. Which information can be found in these files?

In our GTF file everything except the score and frame could be found.

### c) How many genes are present?

Use grep to skip lines with `#` (comments) and then use AWK to check if third column matches "gene":

```bash
cat genomic.gtf | grep -v "#" | awk '$3=="gene"' | wc -l
```

`-v` of grep selects non-matching lines and `$3` extracts the third column.

### d) Provide me the commands and results for counting the number of sequences in the various fasta files (DNA, RNA, protein).

Commands and programs are found in the section [[practical-notebook#b) How can I count the number of reads in a fastq file? Describe different ways to perform that.]]

### e) Describe differences between different genomic fasta files.

Protein fasta cointains protein sequences (amino acid residues).

```
$ head -n 2 protein.faa
>NP_001007096.1 uncharacterized protein Dmel_CG42637, isoform C [Drosophila melanogaster]
MTRWPFNLLLLLSVAVRDCSNHRTVLTVGYLTALTGDLKTRQGLAISGALTMALDEVNKDPNLLPNVYLDLRWNDTKGDT
```

RNA fasta contains transcripts (mRNA).

```
$ head -n 2 rna.fna
>NM_001007095.3 Drosophila melanogaster uncharacterized protein, transcript variant C (CG42637), mRNA
TCACATATTCAAAATCGGGTAGGTAGTCGCGACGGAAAACGGGAAACGCGGACGAATCGCGGAGCCAGAGAAGCGGTAAA
```

CDS contains nucleotide coding sequences (exons).

```
$ head -n 4 cds_from_genomic.fna
>lcl|NC_004354.4_cds_NP_001096854.1_1 [gene=CG17636] [locus_tag=Dmel_CG17636] [db_xref=FLYBASE:FBpp0111834,GeneID:5740847] [protein=uncharacterized protein, isoform A] [protein_id=NP_001096854.1] [location=complement(join(124464..125409,125495..126259,126626..126630))] [gbkey=CDS]
ATGTCGTGCTGCAAGAAGTACGCCGTCTGCTGGATTATCCTGGTGGTGACCGCATTGGGTGTGACCTTGGGTCTGGTTTT
```

# Alignment using HISAT2

HISAT2 is an alignment program for mapping NGS reads against a single reference genome. It outputs alignments in the SAM format. It will build an index over our genome.

Install it with:

```bash
conda install -c bioconda hisat2
```

The command we'll use is `hisat2-build`, which will index our reference genome. It takes 2 parameters: 

- `reference_in` : comma-separated list of files with ref sequences.
- `hisat2_index_base` : writes output files to this *base*name.

Move into `ref_genome` and run:

```bash
hisat2-build \
GCF_000001215.4.fna \
genome \
-p 8
```

`-p` is the number of threads we wan't hisat2-build to spawn (more threads = more speed, but depends on your CPU, so be mindful).

The indexes are created in `genome.X.ht2` files. After creating an index, use the `hisat2` command to create pairwise alignments with our filtered reads. Move back to root and run:

```bash
hisat2 \
-x ref_genome/genome \
-U fastp_results/filtered-reads.fq \
-S hisat_alignment/HisatAlignment.sam \
-p 8
```

- `-x` : prefix of files with indexes.
- `-U` : files with unpaired reads.
- `-S` : file for SAM output.

# Quantifying the expression of transcripts using RNA- seq data with Salmon

[[quantifying the expression of transcripts using RNA-seq data with salmon]]
