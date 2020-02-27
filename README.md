# cutAnd
Snakemake pipeline for cut and tag QC and peak calling

This pipeline will do the inital QC steps for cut&Tag or cut&Run data. 

# SETUP 

Clone this repository into your project directory:

```
git clone git@github.com:maxsonBraunLab/cutAnd.git my-project

# create a directory to hold your data
cd my-project 

# create directory for your fastq files
mkdir -p data/raw

# link your fastqs to here
ln -s /path/to/fastq/files/* data/raw

```

The configuration happens in config.yml 

 - Specify the bowtie index for the genome you are aligning to. 
 - List a control/background sample for peak calling, include in the sample list.
 - A list of sample names, the names must be a continuous portion of the existing filename, as the snakefile uses glob.glob(samplename) to fetch fastq files that match the given sample name. At the cost of flexibility, this is fairly convenient. Maybe it will change to a sample sheet format in the future. 

# Execution

Setup snakemake profile to run on compute cluster:

SLURM: follow [these instructions](https://github.com/Snakemake-Profiles/slurm)

Example of how to run pipeline

```
snakemake --use-conda --profile slurm -j 60 --latency-wait 60
```

To change runtime parameters for indvidual rules, you can provide a cluster configuration file via the `--cluster-config` flag, for example:

```
snakemake --use-conda --profile slurm --cluster-config src/cluster.yml -j 60 --latency-wait 60
```

Here is an example cluster.yml that sets defaults, and rule specific resources:

```
cat src/cluster.yml
__default__:
    memory: "20G"
    threads: 2
    partition: "exacloud"
    time: 14400 #time in seconds 
bowtie2:
    threads: 8
    time: 86400
```


The final QC document will be generated under `data/multiqc/multiqc_report.html`, were alignment statistics, duplication rates, and fastqc results can be compared across all samples. Additionally, a table of peak counts for each sample is output in data/macs2/peak_counts.tsv

