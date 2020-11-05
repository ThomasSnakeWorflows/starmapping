
import sys
import re
import yaml

from os.path import join
import pandas as pd
# validate is only available from sbakemake 5.1
#from snakemake.utils import validate


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


def get_fastq(wildcards):
    if hasattr(wildcards, "orient"):
        r = samples.loc[wildcards.sample, wildcards.orient]
        return {'R': r}
    else:
        r1, r2 = samples.loc[wildcards.sample, ['R1', 'R2']].dropna()
        return {'R1': r1, 'R2' : r2}


def get_fastq_files(wildcards):
    fastqs = get_fastq(wildcards)
    files = []
    for reads in fastqs.values():
        files.extend(reads.split(','))
    return files


configfile: "config.yaml"

modes = config['modes']
print(modes)

samples = pd.read_table(config["samples"], comment='#', dtype=str).set_index("sample", drop=False)
#validate(samples, "samples.schema.yaml")

workdir: config['wdir']
message("The current working directory is " + config['wdir'])

if "pe" in modes or "persem" in modes:
    orientations = ['R1', 'R2']
else:
    orientations = ['R1']


wildcard_constraints:
    sample="|".join(samples.index),


rule all:
    input:
        expand("{mode}.done", mode=modes)

rule se:
    input:
        expand("{sample}/se/{orient}/Log.final.out",
               sample=samples.index, orient=orientations),
        expand("{sample}/se/{orient}/Aligned.sortedByCoord.out.bam.bai",
               sample=samples.index, orient=orientations)
    output:
        touch("se.done")

rule pe:
    input:
        expand("{sample}/pe/Log.final.out",
               sample=samples.index),
        expand("{sample}/pe/Aligned.sortedByCoord.out.bam.bai",
               sample=samples.index)
    output:
        touch("pe.done")

rule persem:
    input:
        expand("{sample}/persem/Quant.genes.results",
               sample=samples.index)
    output:
        touch("persem.done")

rule index:
    input:
        "{prefix}/Aligned.sortedByCoord.out.bam"
    output:
        "{prefix}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule starpe:
    input:
        get_fastq_files
    output:
        "{sample}/pe/Log.final.out",
        "{sample}/pe/Aligned.sortedByCoord.out.bam"
    params:
        R1 =  lambda wildcards: get_fastq(wildcards)['R1'],
        R2 =  lambda wildcards: get_fastq(wildcards)['R2'],
        starrefdir = config['starrefdir']
    log:
        "logs/star/{sample}.log"
    threads: 8
    shell:
        "STAR "
        " --runThreadN {threads}"
        " --genomeDir {params.starrefdir}"
        " --outSAMtype BAM SortedByCoordinate"
        " --readFilesIn {params.R1} {params.R2}"
        " --readFilesCommand zcat "
        " --outFileNamePrefix {wildcards.sample}/pe/"
        " --outReadsUnmapped Fastx"
        " --outSJfilterOverhangMin 15 15 15 15"
        " --alignSJoverhangMin 15"
        " --alignSJDBoverhangMin 15"
        " --outFilterMultimapNmax 20"
        " --outFilterScoreMin 1"
        " --outFilterMatchNmin 1"
        " --outFilterMismatchNmax 2"
        " --chimSegmentMin 15"
        " --chimScoreMin 15"
        " --chimScoreSeparation 10"
        " --chimJunctionOverhangMin 15"

rule starrsem:
    input:
        get_fastq_files,
        gtf = config['ref']['gtf']
    output:
        "{sample}/persem/Log.final.out",
        "{sample}/persem/Aligned.sortedByCoord.out.bam",
        "{sample}/persem/Aligned.toTranscriptome.out.bam"
    params:
        R1 =  lambda wildcards: get_fastq(wildcards)['R1'],
        R2 =  lambda wildcards: get_fastq(wildcards)['R2'],
        starrefdir = config['starrefdir']
    log:
        "logs/star/persem/{sample}.log"
    threads: 4
    shell:
        "STAR --genomeDir {params.starrefdir} STARgenomeDir "
        " --readFilesIn {params.R1} {params.R2}"
        " --readFilesCommand zcat --outFilterType BySJout"
        " --outSAMunmapped Within --sjdbGTFfile {input.gtf}"
        " --outSAMattrIHstart 0"
        " --outFilterIntronMotifs RemoveNoncanonical" # for cufflinks
        " --alignSoftClipAtReferenceEnds Yes" # for cufflinks
        " --outSAMstrandField intronMotif " # for unstranded RNAseq data
        " --outSAMattributes NH HI AS NM MD"
        " --outFilterMultimapNmax 20"
        " --outFilterMismatchNmax 999"
        " --alignIntronMin 20"
        " --alignIntronMax 1000000"
        " --alignMatesGapMax 1000000"
        " --alignSJoverhangMin 8  --alignSJDBoverhangMin 1 "
        " --runThreadN {threads}"
        " --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"
        " --outWigType bedGraph  --outWigStrand Unstranded"
        " --outFileNamePrefix {wildcards.sample}/persem/ 2> {log}"


rule rsem:
    input:
        "{sample}/persem/Aligned.toTranscriptome.out.bam"
    output:
         "{sample}/persem/Quant.genes.results"
    params:
        rsem_prefix = config['rsemprefix']
    log:
        "logs/rsem/{sample}.log"
    threads: 4
    shell:
        "rsem-calculate-expression -p {threads} --paired-end --bam --estimate-rspd "
        "--calc-ci --no-bam-output --seed 12345 "
        "-p {threads} {input} {params.rsem_prefix} "
        "{wildcards.sample}/Quant 2> {log} "


rule starse:
    input:
        get_fastq_files
    output:
        "{sample}/se/{orient}/Log.final.out",
        "{sample}/se/{orient}/Aligned.sortedByCoord.out.bam"
    params:
        R =  lambda wildcards: get_fastq(wildcards)['R'],
        starrefdir = config['starrefdir']
    log:
        "logs/star/{sample}_{orient}.log"
    threads: 8
    shell:
        "STAR "
        " --runThreadN {threads}"
        " --genomeDir {params.starrefdir}"
        " --outSAMtype BAM SortedByCoordinate"
        " --readFilesIn {params.R}"
        " --readFilesCommand zcat"
        " --outFileNamePrefix {wildcards.sample}/se/{wildcards.orient}/"
        " --outReadsUnmapped Fastx"
        " --outSJfilterOverhangMin 15 15 15 15"
        " --alignSJoverhangMin 15"
        " --alignSJDBoverhangMin 15"
        " --seedSearchStartLmax 30"
        " --outFilterMultimapNmax 20"
        " --outFilterScoreMin 1"
        " --outFilterMatchNmin 1"
        " --outFilterMismatchNmax 2"
        " --chimSegmentMin 15"
        " --chimScoreMin 15"
        " --chimScoreSeparation 10"
        " --chimJunctionOverhangMin 15"
