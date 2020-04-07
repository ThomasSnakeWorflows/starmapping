
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
        return samples.loc[wildcards.sample, wildcards.orient]
    else:
        return samples.loc[wildcards.sample, ['R1', 'R2']].dropna()


configfile: "config.yaml"

modes = config['modes']

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
#validate(samples, "samples.schema.yaml")

workdir: config['wdir']
message("The current working directory is " + config['wdir'])

if "pe" in modes:
    orientations = ['R1', 'R2']
else:
    orientations = ['R1']

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

rule index:
    input:
        "{prefix}/Aligned.sortedByCoord.out.bam"
    output:
        "{prefix}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule starpe:
    input:
        get_fastq
    output:
        "{sample}/pe/Log.final.out",
        "{sample}/pe/Aligned.sortedByCoord.out.bam"
    params:
        starrefdir = config['starrefdir']
    log:
        "logs/star/{sample}.log"
    threads: 8
    shell:
        "STAR "
        " --runThreadN {threads}"
        " --genomeDir {params.starrefdir}"
        " --outSAMtype BAM SortedByCoordinate"
        " --readFilesIn {input}"
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

rule starse:
    input:
        get_fastq
    output:
        "{sample}/se/{orient}/Log.final.out",
        "{sample}/se/{orient}/Aligned.sortedByCoord.out.bam"
    params:
        starrefdir = config['starrefdir']
    log:
        "logs/star/{sample}_{orient}.log"
    threads: 8
    shell:
        "STAR "
        " --runThreadN {threads}"
        " --genomeDir {params.starrefdir}"
        " --outSAMtype BAM SortedByCoordinate"
        " --readFilesIn {input}"
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
