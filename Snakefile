# cutAnd snakemake pipeline for Cut&Run or Cut&Tag data
import glob
import os

# import configuration
configfile: "src/config_run.yml"

# map samples to fastqs
def get_samples(dir):
    '''Matches samples to their fastq files.'''
    samples = config["SAMPLES"]
    sdict = {s: glob.glob(f"{dir}/*{s}*") for s in samples}
    return sdict

# construct samples matching their fastqs dict
sampdict = get_samples("data/raw")
readpair=[[k+"_R1", k+"_R2"] for k in sampdict.keys()]
reads=[v for s in readpair for v in s]


def get_control(sample):
    '''Gets the control file ("IgG") for the sample'''
    t = sample.split("_")[0]
    c = glob.glob(f"data/mrkdup/{t}_IgG*.bam")
    if len(c)==0:
        control=""
    elif os.path.basename(c[0])==sample:
        control=""
    elif len(c)==1:
        control="-c "+str(c[0]) 
    else:
        print(f"Control file couldn't be found for {sample}")
        control=""
    return control

# define target output
rule all:
    input:
         "data/multiqc/multiqc_report.html",
         expand("data/fastqc/{read}.html", read=reads),
         expand([
                "data/macs2/{sample}_treat_pileup.bdg",
                "data/macs2/{sample}_control_lambda.bdg",
                "data/macs2/{sample}_FE.sort.bw"
                ],
                 sample=sampdict.keys()),

rule fastqc:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        html="data/fastqc/{read}.html",
	zip="data/fastqc/{read}_fastqc.zip" 
    log:
        "data/logs/fastqc_{read}.log"
    threads: 4
    wrapper:
        "0.49.0/bio/fastqc"
        

# align samples to genome
rule bowtie2:
    input:
        lambda wildcards: sampdict[wildcards.sample]
    output:
        "data/aligned/{sample}.bam"
    log:
        err="data/logs/bowtie2_{sample}.err"
    conda:
        "envs/align.yml"
    threads: 8
    shell:
        "bowtie2 --local --very-sensitive-local "
        "--no-unal --no-mixed --threads {threads} "
        "--no-discordant --phred33 "
        "-I 10 -X 700 -x {config[GENOME_IDX]} "
        "-1 {input[0]} -2 {input[1]} 2>{log.err} | samtools view -Sbh - > {output}"

# sort
rule sort:
    input:
        rules.bowtie2.output 
    output:
        temp("data/sort/{sample}.sorted.bam")
    conda:
        "envs/sam.yml"
    log: 
        "data/logs/sort_{sample}.log"
    threads: 8 
    shell:
        "samtools sort -T {wildcards.sample} -@ {threads} -o {output} {input} > {log} 2>&1"

# mark pcr duplicates
rule mrkdup:
    input:
        rules.sort.output
    output:
        "data/mrkdup/{sample}.sorted.mrkdup.bam"
    conda:
        "envs/sam.yml"
    log:
        "data/logs/mrkdup_{sample}.log"
    threads: 4
    shell:
        "sambamba markdup -t {threads} {input} {output} > {log} 2>&1"

# generate index file
rule index:
    input:
       rules.mrkdup.output 
    output:
        "data/mrkdup/{sample}.sorted.mrkdup.bam.bai"
    conda:
        "envs/sam.yml"
    log: 
        "data/logs/index_{sample}.log"
    threads: 4
    shell:
        "sambamba index -t 2 {input} > {log} 2>&1"


# keep all dups
rule macs2:
    input:
        rules.mrkdup.output 
    output:
        "data/macs2/{sample}_treat_pileup.bdg", "data/macs2/{sample}_control_lambda.bdg"
    conda:
        "envs/macs2.yml"
    params:
        outdir="data/macs2",
        control=lambda wildcards: get_control(wildcards.sample)
    log:
        "data/logs/macs2_{sample}.log"
    shell:
        """
        # macs2 call broad peaks no dups 
        macs2 callpeak \
                    -t {input} \
                    {params.control} \
                    -n {wildcards.sample} \
                    --outdir {params.outdir} \
                    -g {config[MACS2][G]} \
                    -q 0.01 \
                    --bdg \
                    --SPMR \
                    -f BAMPE \
                    --keep-dup all \
                    --broad \
                    --nomodel \
        """ 

rule get_fold_enrichment:
    input:
        "data/macs2/{sample}_treat_pileup.bdg","data/macs2/{sample}_control_lambda.bdg"
    output:
        "data/macs2/{sample}_FE.bdg"
    conda:
        "envs/macs2.yml"
    params:
        ""
    log:
        "data/logs/macs2_bdgcmp_{sample}.log"    
    shell:
        # calculate fold enrichment track
        "macs2 bdgcmp -t {input[0]} -c {input[1]} -o {output} -m FE"

rule sort_bg:
    input:
        rules.get_fold_enrichment.output  
    output:
        "data/macs2/{sample}_FE.sort.bdg"
    shell:
        "LC_COLLATE=C; sort -k1,1 -k2,2n {input} > {output}"

rule bdg_to_bw:
    input:
        rules.sort_bg.output
    output:
        "data/macs2/{sample}_FE.sort.bw"        
    conda:
        "envs/bgtobw.yml"
    log:
        "data/logs/bgtobw_{sample}.log"
    shell:
        "bedGraphToBigWig {input} {config[CHRSIZES]} {output}"

rule multiqc:
    input:
        expand("data/macs2/{sample}_peaks.xls", sample=sampdict.keys())
    output:
        "data/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    log:
        "data/logs/multiqc.log"
    shell:
        "multiqc --force -o data/multiqc data/ > {log} 2>&1"
    

