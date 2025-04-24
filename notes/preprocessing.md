Trimming reads (adapters). Tools: `cutadapt`, `fastp` (can automatically detect adapter sequences). They can also remove reads that are of too low quality.

```bash
$ conda install -bioconda fastp
```

# Run fastp

```bash
$ fastp -i SRR30833097.fastq -o fastp-out.fq
```

##### `-i SRR30833097.fastq`

Input fastq file.

##### `-o fastp-out.fq`

Name of the output fastq file.

This command outputs:

```
Detecting adapter sequence for read1...
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

Duplication rate (may be overestimated since this is SE data): 32.9681%

JSON report: fastp.json
HTML report: fastp.html

fastp -i SRR30833097.fastq -o fastp-out.fq 
fastp v0.22.0, time used: 17 seconds
```

SInce our fastqc report showed a lot of polyA tails, we can use the flag `--trim_poly_x`, which should enable trimming on the 3' end.

```bash
$ fastp -i SRR30833097.fastq -o fastp-out.fq --trim_poly_x
```

Interestingly this doesn't do almost anything for our trailing polyA tails. The minimum tail length it's searching for is 10, but lowering it to 5 made it worse. I kept the min tail length and trimmed 3' end tails based on their qualities.  The `--cut_tail` creates a sliding window (size 4 by default) that removes bases if their mean quality drops below the given threshold. 20 increased the quality a little bit, but 30 gave us a passing grade in fastqc.

```bash
$ fastp -i SRR30833097.fastq -o fastp-out.fq --trim_poly_x --poly_x_min_len 10 --cut_tail --cut_tail_mean_quality 30 --html fastp-report
```

Run `$ fastqc` and open the newly generated file `fastp-out.fq`. Save the report as `post-filtering_fastqc.html`.

![[Pasted image 20250424144014.png]]

# 