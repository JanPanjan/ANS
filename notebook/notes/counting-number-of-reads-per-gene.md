# Converting SAM to BAM format and sorting sequences according to the genomic coordinates

After the alignment is produced, reads are in **random order** with respect to their **position in the reference genome**. To examine the alignment with IGV or Tablet, BAM should be **sorted** and **indexes** of BAM and FASTA files should be created.

First, format the SAM file into BAM.

```bash
cd hisat_alignment # move into this directory first
```

```bash
samtools view \
-o HisatAlignment.bam \
HisatAlignment.sam \
--threads 16
```

- `-o` : output BAM file
- `--threads` : threads to spawn

Then sort the alignments in the bam file based on their leftmost coordinates.

```bash
samtools sort \
-o HisatAlignment.sorted.bam \
HisatAlignment.bam \
--threads 16
```

Create an index for this BAM file.

```bash
samtools index HisatAlignment.sorted.bam -@ 16 # output is .bai file
```

Create a new directory for this samtools index.

```bash
cd .. # move back to root
mkdir -p ref_genome/samtools_index
cd ref_genome
```

Create the index.

```bash
samtools faidx GCF*.fna
mv GCF*.fna.fai samtools_index # move the index
```

---

<< [[preprocessing, aligning, quantifying expression]] | [[]] >>