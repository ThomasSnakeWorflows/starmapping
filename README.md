# starmapping

```bash
module load bioinfo/STAR-2.5.1b
conda activate star
```

```bash
snakemake --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```bash
