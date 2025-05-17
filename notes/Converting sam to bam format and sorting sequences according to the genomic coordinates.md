Samtools package conflicts with hisat2, so create a new environment.

## collapse this

```
LibMambaUnsatisfiableError: Encountered problems while solving:
  - package python-3.13.2-hf623796_100_cp313 requires openssl >=3.0.15,<4.0a0, but none of the providers can be installed

Could not solve for environment specs
The following packages are incompatible
├─ hisat2 =* * is installable with the potential options
│  ├─ hisat2 [2.0.0beta|2.0.1beta|...|2.1.0] would require
│  │  └─ python [=2.7 *|>=2.7,<2.8.0a0 *], which can be installed;
│  ├─ hisat2 [2.0.0beta|2.0.1beta|...|2.0.5] would require
│  │  └─ python =3.4 *, which conflicts with any installable versions previously reported;
│  ├─ hisat2 [2.0.0beta|2.0.1beta|...|2.1.0] would require
│  │  └─ python [=3.5 *|>=3.5,<3.6.0a0 *], which can be installed;
│  ├─ hisat2 [2.0.5|2.1.0] would require
│  │  └─ perl =5.22.0 *, which conflicts with any installable versions previously reported;
│  ├─ hisat2 2.1.0 would require
│  │  └─ python >=3.6,<3.7.0a0 *, which can be installed;
│  ├─ hisat2 2.1.0 would require
│  │  └─ python >=3.7,<3.8.0a0 *, which can be installed;
│  ├─ hisat2 2.2.0 would require
│  │  └─ python_abi =2.7 *_cp27mu, which conflicts with any installable versions previously reported;
│  ├─ hisat2 2.2.0 would require
│  │  └─ python_abi =3.6 *_cp36m, which conflicts with any installable versions previously reported;
│  ├─ hisat2 2.2.0 would require
│  │  └─ python_abi =3.7 *_cp37m, which conflicts with any installable versions previously reported;
│  ├─ hisat2 2.2.1 would require
│  │  └─ python >3.5 * with the potential options
│  │     ├─ python [3.10.0|3.10.1|...|3.9.9], which can be installed;
│  │     ├─ python [3.13.0|3.13.1|3.13.2] would require
│  │     │  └─ openssl >=3.0.15,<4.0a0 *, which can be installed;
│  │     ├─ python [3.5.1|3.5.2|...|3.5.6], which can be installed;
│  │     ├─ python [3.6.0|3.6.1|...|3.6.9], which can be installed;
│  │     └─ python [3.7.0|3.7.1|...|3.7.9], which can be installed;
│  ├─ hisat2 2.2.1 would require
│  │  └─ libgcc >=13 *, which conflicts with any installable versions previously reported;
│  ├─ hisat2 2.2.1 would require
│  │  └─ libgcc-ng >=12 *, which conflicts with any installable versions previously reported;
│  ├─ hisat2 2.2.1 would require
│  │  └─ libgcc >=12 *, which conflicts with any installable versions previously reported;
│  └─ hisat2 2.2.1 would require
│     └─ python_abi =3.8 *_cp38, which conflicts with any installable versions previously reported;
├─ pin on python 3.13.* =* * is not installable because it requires
│  └─ python =3.13 *, which conflicts with any installable versions previously reported;
└─ samtools =* * is not installable because there are no viable options
   ├─ samtools [0.1.16|0.1.17|...|1.7] would require
   │  └─ openssl >=1.1.0,<=1.1.1 *, which conflicts with any installable versions previously reported;
   ├─ samtools [0.1.12|0.1.13|...|1.9] would require
   │  └─ ncurses [=5.9 *|==5.9 *|>=5.9,<5.10.0a0 *], which conflicts with any installable versions previously reported;
   ├─ samtools [0.1.18|0.1.19|...|1.6] would require
   │  └─ libgcc-ng >=12 *, which conflicts with any installable versions previously reported;
   ├─ samtools [0.1.18|0.1.19|...|1.6] would require
   │  └─ libgcc >=13 *, which conflicts with any installable versions previously reported;
   ├─ samtools 1.10 would require
   │  └─ htslib >=1.10.2,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools 1.10 would require
   │  └─ htslib >=1.10,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools 1.11 would require
   │  └─ htslib >=1.11,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools 1.12 would require
   │  └─ htslib >=1.12,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools 1.13 would require
   │  └─ htslib >=1.13,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools [1.14|1.15] would require
   │  └─ htslib >=1.14,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools [1.15|1.15.1] would require
   │  └─ htslib >=1.15,<1.22.0a0 *, which conflicts with any installable versions previously reported;
   ├─ samtools 1.21 would require
   │  └─ libgcc >=12 *, which conflicts with any installable versions previously reported;
   ├─ samtools [1.3|1.3.1|1.6] would require
   │  └─ openssl >=1.1.1l,<1.1.2a *, which conflicts with any installable versions previously reported;
   ├─ samtools [1.5|1.6|1.7] would require
   │  └─ zlib [=1.2.8 *|==1.2.8 *], which conflicts with any installable versions previously reported;
   ├─ samtools 1.9 would require
   │  └─ htslib >=1.9,<1.10.0a0 *, which conflicts with any installable versions previously reported;
   └─ samtools 1.9 would require
      └─ libdeflate >=1.0,<1.1.0a0 *, which conflicts with any installable versions previously reported.

Pins seem to be involved in the conflict. Currently pinned specs:
 - python=3.13
```

## new env

```bash
conda create -n sam
conda activate sam
```

```bash
conda install -c bioconda samtools
conda install -c bioconda tablet
```

After the alignment is produced, reads are in random order with respect to their position in the
reference genome. To examine the alignment with IGV or Tablet bam should be sorted and
indexes of bam and fasta files should be created.

`cd` into  `ref_genome/data/GCF_000001215.4` where the SAM file is located. If it's not, you can run from your root directory:

```
find . -type f -name "*.sam"
```

**Samtools** is  a  set of utilities that manipulate alignments in the BAM format. It imports from and
exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to
retrieve reads in any regions swiftly.

## Convert sam to bam

Use the `-o` flag to output the `sam`  file into `bam` format with command:

> [!bug] Samtools in bioconda
> neki nedela s samtools če installaš z bioconda, z apt dela.

```bash
samtools view -o HisatAlignment.bam HisatAlignment.sam
```

## Sort the bam file

Sort the bam file according to the location of reads on the reference genome.

```bash
samtools sort -o HisatAlignment.sorted.bam HisatAlignment.bam
```

## Create a new index

Tablet requires an index of the bam file to visualize the alignment. Create an index with:

```bash
samtools index HisatAlignment.sorted.bam
```

Create a new index for the reference genome also.

```bash
samtools faidx GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
```

The index is created as a `fai` file.

```bash
ls *.fai
```

> [!tip] `samtools faidx <ref.fasta> [region1 [...]]`
>
> Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence.
> If no region is specified, faidx will index the file and create `<ref.fasta>.fai` on the disk.

## Information about the SAM file

It can be accessed with:

```bash
samtools flagstat *.sam # this works since we have a single '.sam' file
```

It outputs this:

```
6703266 + 0 in total (QC-passed reads + QC-failed reads)
4473123 + 0 primary
2230143 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
6301123 + 0 mapped (94.00% : N/A)
4070980 + 0 primary mapped (91.01% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

> [!tip] `samtools flagstat in.sam|in.bam|in.cram`
>
> Does a full pass through the input file to calculate and print statistics to stdout.
>
> Provides  counts  for each of 13 categories based primarily on bit flags in the FLAG field. Each category in the output is broken down into
> QC pass and QC fail, which is presented as "#PASS + #FAIL" followed by a description of the category.

For example, description for flag 163 can be obtained by:

```bash
samtools flags 163
```

Read about the SAM file [here](https://www.metagenomics.wiki/tools/samtools/bam-sam-file-format).

> flags are explained too.

## Check the alignment with tablet or igv

I'll be using IGV as the viewing tool. The web-app is available [here](https://igv.org/app/).
Read about IGV [here](https://igv.org/doc/desktop/#FileFormats/DataTracks/#bam) and [here](https://igv.org/doc/webapp/#UserGuide/#loading-the-reference-genome).

### Include a screenshot of the sorted.bam file to confirm that the reads were sorted based on the reference genome position.

Not sure if this is correct.

![[~/Pictures/Screenshots/Screenshot From 2025-05-15 14-55-35]]

### Which column in the sam or bam file contains the leftmost mapping position?

### How to retrieve the unmapped reads from a bam file?