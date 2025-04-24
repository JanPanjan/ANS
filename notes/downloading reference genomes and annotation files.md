`GTF` short for genome transfer files, that contain genome annotation.
Find the reference genome [[navodila]] like here.
The accession number is GCF_000001215.4.

Save this in a file:

```bash
$ echo "GCF_000001215.4" > RefGenAcc.txt
```

Install NCBI toolset:

```bash
$ conda install conda-forge::ncbi-datasets-cli
```

With this we can install our reference genome by passing our accession number:

```bash
$ datasets download genome \
> accession GCF_000001215.4 \
> --include genome,rna,protein,cds,gff3,gtf \
> --filename genome.zip
```

Unzip and rename and remove the zip of the downloaded genome:

```bash
$ unzip genome.zip # creates a directory ncbi_dataset
$ mv ncbi_dataset ref_genome
$ rm genome.zip
```