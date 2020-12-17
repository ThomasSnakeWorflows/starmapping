import sys
import re
import os


def get_threads(rule, default):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "threads" in cluster_config[rule]:
        return cluster_config[rule]["threads"]
    if "default" in cluster_config and "threads" in cluster_config["default"]:
        return cluster_config["default"]["threads"]
    return default


def get_all_output(wildcards):
    outputs = []
    if "starrefdir" in config:
        outputs.append(config['starrefdir']+"/SAindex")
    if "rsemprefix" in config:
        outputs.append(config['rsemprefix']+"/RSEMref.seq")
    return outputs


configfile: "config.yaml"

reference_fasta = os.path.abspath(config['ref']['fasta'])
reference_gtf = os.path.abspath(config['ref']['gtf'])

rule prepare:
    input:
        get_all_output

rule star_prepare:
    input:
        fa=reference_fasta,
        gtf=reference_gtf
    output:
        config['starrefdir']+"/SAindex"
    params:
        starrefdir=config['starrefdir']
    log:
        "logs/star_prepare.log"
    threads: get_threads("star_prepare", 4)
    shell:
        "STAR"
        " --runThreadN {threads}"
        " --runMode genomeGenerate"
        " --genomeDir '{params.starrefdir}'"
        " --genomeFastaFiles {input.fa}"
        " --sjdbGTFfile {input.gtf}"
        " --sjdbOverhang 100"
        " --outFileNamePrefix {params.starrefdir}"
        " --outStd Log"
        " 2> {log}"

if "rsemprefix" in config:
    rule rsem_prepare:
        input:
            fa=reference_fasta,
            gtf=reference_gtf
        output:
            config['rsemprefix']+"/RSEMref.seq"
        params:
            refsemrefprefix=config['rsemprefix']
        log:
            "logs/rsem_prepare.log"
        shell:
            "rsem-prepare-reference --gtf {input.gtf} {input.fa} {params.refsemrefprefix}"
