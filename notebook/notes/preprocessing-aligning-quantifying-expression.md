---
id: preprocessing-aligning-quantifying-expression
aliases: []
tags: []
---

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

> [!WARNING]
> Although it says that max is 36 on the github page, the command will fail with this value and print:
> ```bash
>  ERROR: the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.
> ```
> So we will set it to 30.

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

> [!NOTE] *per base sequence content* and *per sequence GC content*
> These two modules still fail, but if we think logically, the first one can be dependent on the other.
> Since our reads have skewed distribution of bases, the second one fails and with it the first one fails too, but this isn't critical in my opinion.

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

Commands and programs are found in the section [[faks/ANS/notebook/README#b) How can I count the number of reads in a fastq file? Describe different ways to perform that.]]

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

Make a directory for the index and move into it.

```bash
mkdir -p ref_genome/hisat_index # run this in root
cd ref_genome/hisat_index
```

The command we'll use is `hisat2-build`, which will index our reference genome. It takes 2 parameters:

- `reference_in` : comma-separated list of files with ref sequences.
- `hisat2_index_base` : writes output files to this *base*name.

```bash
hisat2-build \
../GCF_000001215.4.fna \
genome \
-p 16
```

`-p` is the number of threads we wan't hisat2-build to spawn (more threads = more speed, but depends on your CPU, so be careful).

The indexes are created in `genome.X.ht2` files. After creating an index, use the `hisat2` command to create pairwise alignments with our filtered reads. Move back to root and run:

```bash
hisat2 \
-x ref_genome/hisat_index/genome \
-U fastp_results/filtered-reads.fq \
-S hisat_alignment/HisatAlignment.sam \
--summary-file hisat_alignment/alignment_summary.txt \
-p 8
```

- `-x` : prefix of files with indexes.
- `-U` : files with unpaired reads.
- `-S` : file for SAM output.

This is the alignment summary:

```
4473488 reads; of these:
  4473488 (100.00%) were unpaired; of these:
    376692 (8.42%) aligned 0 times
    3460208 (77.35%) aligned exactly 1 time
    636588 (14.23%) aligned >1 times
91.58% overall alignment rate
```

77.35 % of reads (3'460'208) were mapped uniquelly, 14.23 % (636'588) were mapped to multiple loci and 8.42 % (376'692) were not mapped.

# Quantifying the expression of transcripts using RNA-seq data with Salmon

[[quantifying the expression of transcripts using RNA-seq data with salmon]]
Salmon is a tool for transcript quantification from RNA-seq data. It typically works in two steps - **indexing** and **quantification**, if ran in mapping mode or only quantification if ran in alignment mode (with precomputed alignments).

We want to quantify how many of our transcripts are actually coding sequences (i.e. for proteints). We'll do this by using the downloaded genome **coding nucleotide sequences (cds)**. Since we haven't aligned any reads to them, we'll use the mapping route.

Install Salmon it with:

```bash
conda install -c bioconda salmon=1.10.3
```

Create a new directory.

```bash
mkdir quantification
```

## 1. Indexing

First, create an index for the transcriptome.

```bash
salmon index \
-t ref_genome/cds_from_genomic.fna \
-i quantification/cds_index \
-p 16
```

- `-t` : FASTA file with target transcripts
- `-i` : directory for the index (makes it if it doesn't exist)
- `-p` : number of threads to spawn

## 2. Quantification

Second, quantify the reads against this index.

```bash
salmon quant \
-i quantification/cds_index \
-l A \
-r fastp_results/filtered-reads.fq \
-o quantification \
-p 16
```

- `-i` : input index
- `-l` : library type
- `-r` : input **unpaired (single-end)** alignments
- `-o` : output directory
- `-p` : number of threads to spawn

> Since version 0.7.0, Salmon can automatically infer the type of the seqencing library, hence we use A (auto).

## Output

The output is a quantification file (`.sf`), plain text, tab-separated file with a single header. It's named `quant.sf`.
