# starmapping

```bash
module load bioinfo/STAR-2.6.0c
module load bioinfo/RSEM-1.3.3
module load bioinfo/snakemake-4.8.0
```

```bash
snakemake --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```bash

If star index is not available
```bash
snakemake -s starindex.smk --jobs 2 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```bash
