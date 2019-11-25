# This workflow performs GATK haplotypecaller on RNASeq data
# based on GATK best practice (for RNASeq), 
# followed by demuxlet on a scRNASeq library composing multiple
# samples to be de-mutiplexed. 
#
# Input/output/parameters should be specified in config.yaml
#
#

from snakemake.utils import validate, min_version
from snakemake.shell import shell

min_version("5.1.2")

configfile: "config.yaml"
#validate(config, schema="config.schema.yaml")

threads_max = config["params"]["threads"]

SAMPLES = config["sample"].split(",")

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

rule all:
    input:
        [ ".temp" if config["mode"] == "demuxlet" else "merged.vcf"]

rule markdup:
    input:
        "{sample}.bam"
    output:
        temp("{sample}.mkdup.bam")
    threads: threads_max
    benchmark:
        "benchmark_{sample}.tsv"
    shell:
        "gatk MarkDuplicatesSpark --input {input} --output {output}"


rule split_n:
    input:
        "{sample}.mkdup.bam"
    output:
        temp("{sample}.mkdup.splitN.bam")
    threads: 
        threads_max
    params:
     	fa=config["ref"]["fa"]
    shell:
        "gatk SplitNCigarReads --input {input} --output {output} --reference {params.fa}"
        
rule base_recalibrator:
    input:
        "{sample}.mkdup.splitN.bam"
    output:
        temp("{sample}.table")
    threads: threads_max
    params:
        fa=config["ref"]["fa"],
        dbsnp=config["known"]["dbsnp"],
        snp1k=config["known"]["snp1k"],
        knowindel=config["known"]["knowindel"],
        indel1k=config["known"]["indel1k"]
    shell:
          "gatk BaseRecalibrator -I {input} -O {output} -R {params.fa} "
          "--use-original-qualities --known-sites {params.dbsnp} "
          "--known-sites {params.snp1k} --known-sites {params.indel1k} "
          "--known-sites {params.knowindel}"

rule bqsr:
    input:
        table="{sample}.table",
        bam="{sample}.mkdup.splitN.bam"
    output:
        temp("{sample}.mkdup.splitN.bqsr.bam")
    threads: threads_max
    params:
        fa=config["ref"]["fa"]
    shell:
        "gatk ApplyBQSR --use-original-qualities -I {input.bam} "
        "-O {output} -R {params.fa} "
        "--bqsr-recal-file {input.table}"

rule haplotype_caller:
    input:
        "{sample}.mkdup.splitN.bqsr.bam"
    output:
        "{sample}.haplotypecaller.raw.g.vcf.gz"
    threads: threads_max
    params:
        fa=config["ref"]["fa"]
    shell:
        "gatk --java-options \"-Xmx8G\" HaplotypeCaller "
        "-R {params.fa} -I {input} -O {output} "
        "--dont-use-soft-clipped-bases true "
        "-stand-call-conf 20.0 "

rule filter_vcf:
    input:
        "{sample}.haplotypecaller.raw.g.vcf.gz"
    output:
        "{sample}.haplotypecaller.filtered.g.vcf.gz"
    threads: threads_max
    params:
        fa=config["ref"]["fa"],
        FS=config["vcf_filter"]["FS"],
        QD=config["vcf_filter"]["QD"]
    shell:
        "gatk VariantFiltration -V {input} "
        "-R {params.fa} "
        "--filter-name \"FS\" --filter-expression \"FS > {params.FS}\" "
        "--filter-name \"QD\" --filter-expression \"QD < {params.QD}\" "
        "-O {output}"

rule merge_vcfs:
    input:
        expand("{sample}.haplotypecaller.filtered.g.vcf.gz", sample=SAMPLES)
    output:
        "merged.vcf"
    threads: threads_max
    run:
        if len(SAMPLES)>1:
            shell("bcftools merge --threads {threads} -Ov -o {output} {input}")
        else:
            shell("gunzip -c {input} > {output}")

rule demuxlet:
    input:
        bam=config["demuxlet"]["bam"],
        vcf="merged.vcf",
    output:
        temp(".temp")
        #"{}.best".format(config["demuxlet"]["output_prefix"])
    params:
        prefix="{}_cov{}".format(config["demuxlet"]["output_prefix"], config["demuxlet"]["min_coverage"]),
        cov=config["demuxlet"]["min_coverage"],
        mq=config["demuxlet"]["min_mapq"]
    threads: threads_max
    shell:
        "demuxlet --alpha 0 --alpha 0.5 "
        "--sam {input.bam} --vcf {input.vcf} "
        "--field GT --out {params.prefix} "
        "--min-MQ {params.mq} --min-uniq {params.cov} | "
        "touch .temp"



