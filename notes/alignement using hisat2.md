HISAT2 is a splice aware mapper (STAR is an alternative). With it we will create an index for our genome.

https://daehwankimlab.github.io/hisat2/manual/

Install the tool with:

```bash
$ conda install -c bioconda hisat2
```

Move into `ref_genome/data/G*` and run the next command to build an index for our genome:

```bash
$ hisat2-build GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna genome
```

Specify directory of indexes in an environment variable:

```bash
$ HISAT_INDEXES='full-path-to-where-you-ran-the-previous-directory'
# for me, this is /home/XXX/YYY/ZZZ/WORKING-DIRECTORY/ref_genome/data/GCF_000001215.4/
```

>[!bug] Environmental variables
>This did not work. I think the variable has to be in the `.bashrc` (or `.zshrc`)!

After that we can do sequence alignment within the same directory:

```bash
$ hisat2 -x genome -U ../../../fastp-out.fq -S HisatAlignment.sam
```