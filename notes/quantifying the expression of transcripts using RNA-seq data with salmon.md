https://combine-lab.github.io/salmon/

Install salmon using conda.

```bash
$ conda install -c conda-forge salmon
```

Create the index in the `/ref_genome/data/GCF_000001215.4` directory.

```bash
$ mkdir salmon-index # create dir for salmon index
$ salmon index -t cds_from_genomic.fna -i salmon-index/
```

Now quantify transcript abundances.

```bash
$ mkdir salmon-quant # create dir for salmon quantification
$ salmon quant -i salmon-index/ -l A -r ../../../fastp-out.fq -o salmon-quant
```

I'm guessing the `salmon-quant/quant.sf` file is the quantification output.